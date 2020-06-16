import pandas as pd
import numpy as np
import warnings
from scipy import stats
from scipy.stats import multinomial
# import statsmodels.api as sm

try:
    import parasail
except ImportError:
    print('Importing palmotif without alignment support from parasail.')

__all__ = ['compute_motif',
           'compute_pal_motif',
           'uniprot_frequency',
           'compute_relative_motif',
           'consensus_seq']

aa_alphabet = [aa for aa in 'ARNDCQEGHILKMFPSTWYVBZX-']

uniprot_frequency = {'A': 8.25,
                     'R': 5.53,
                     'N': 4.06,
                     'D': 5.45,
                     'C': 1.37,
                     'Q': 3.93,
                     'E': 6.75,
                     'G': 7.07,
                     'H': 2.27,
                     'I': 5.96,
                     'L': 9.66,
                     'K': 5.84,
                     'M': 2.42,
                     'F': 3.86,
                     'P': 4.70,
                     'S': 6.56,
                     'T': 5.34,
                     'W': 1.08,
                     'Y': 2.92,
                     'V': 6.87}

def compute_motif(seqs, reference_freqs=None, weights=None, align_first=False, gap_reduce=None, alphabet=None):
    """Compute heights for a sequence logo based on the frequency of AAs/symbols at each position of the
    aligned sequences. The output matrix/DataFrame contains the relative entropy of each AA/symbol in bits. It
    is computed relative to uniform frequencies (default) or some fixed set of frequencies (not position specific).

    The output is often used as symbol heights in a logo plot. Heights would reflect the fraction of total
    entropy contributed by that AA within each column of the alignment.

    Sequences must be pre-aligned using a tool like ClustalW, MUSCLE or other tool for multiple-sequence alignment.

    Parameters
    ----------
    seqs : list or pd.Series
        Alignment of AA sequences
    reference_freqs : dict or None
        Reference frequencies for each AA
    weights : np.array
        Weights for each sequence
    gap_reduce : None or float threshold
        Identifies columns with threshold fraction of gaps,
        removes the sequences that have AAs there and
        removes the column from the alignment.
    alphabet : str or list
        Limit frequencies to only a subset of AAs
        Default will consider all amino-acids.

    Returns
    -------
    motif : pd.DataFrame
        Heights for each symbol (index, rows) for each column.
        Heights reflect the fraction of total entropy contributed
        by that AA within each column of the alignment."""

    if alphabet is None:
        alphabet = aa_alphabet

    if weights is None:
        weights = np.ones(len(seqs))
    
    if not isinstance(seqs, pd.Series):
        align = pd.Series(seqs)
    else:
        align = seqs

    if not gap_reduce is None:
        align = reduce_gaps(align, thresh=gap_reduce)

    L = len(align.iloc[0])
    nAA = len(alphabet)

    freq = _get_frequencies(align, alphabet, weights)

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')

        if not reference_freqs is None:
            reference_freqs = np.array([reference_freqs.get(aa, 0) for aa in alphabet])
            reference_freqs = reference_freqs / np.sum(reference_freqs)
            reference_freqs = np.tile(reference_freqs[:, None], (1, L))
            heights = freq * np.log2(freq / reference_freqs)
            heights[np.isnan(heights)] = 0
        else:
            tot_entropy = np.log2(nAA) - np.nansum(-freq * np.log2(freq), axis=0, keepdims=True)
            heights = freq * tot_entropy

    return pd.DataFrame(heights, columns=list(range(L)), index=alphabet)

def reduce_gaps(align, thresh=0.7):
    """Identifies columns with thresh fraction of gaps,
    removes the sequences that have AAs there and
    removes the column from the alignment

    Parameters
    ----------
    align : pd.Series
        Contains one row per sequence, with each sequence assumed to
        have the same length and be aligned
    thresh : fraction
        Positions with more than thresh fraction of sequences having
        a gap character will be removed from the alignment and the 1 - thresh
        sequences that don't have a gap at that position will be removed
        from the aligmment

    Returns
    -------
    out : pd.Series
        Similar to input align, but with possibly fewer sequences (rows)
        and shorted sequences."""
    removePos = []
    for pos in range(len(align.iloc[0])):
        if np.mean(align.map(lambda seq: seq[pos] == '-')) > thresh:
            removePos.append(pos)
    for pos in removePos:
        align = align.loc[align.map(lambda seq: seq[pos] == '-')]
    return align.map(lambda seq: ''.join([aa for pos, aa in enumerate(seq) if not pos in removePos]))

def pairwise_alignment_frequencies(centroid, seqs, gopen=3, gextend=3, matrix=parasail.blosum62, alphabet=None):
    """Aligns each seq in seqs to centroid using global Needleman-Wunsch and then counts
    the AA/symbols in each seq at each alignment/centroid position. Alignment positions
    containing a gap/"-" in the centroid sequence are not counted.

    Parameters
    ----------
    centroid : str
        Sequence to align against. It does not need to be the centroid, though that
        would be a sensible choice.
    seqs : iterable
        List of sequences of varying length to be globally aligned to centroid seq
    gopen, gextend, matrix : alignment parameters for parasail

    Returns
    -------
    counts : np.ndarray [position x symbol x seq]
        A binary matrix indicating the presence of a symbol after alignment
        for each position of the centroid, each symbol and each sequence."""
    if alphabet is None:
        alphabet = aa_alphabet

    centroid_query = parasail.profile_create_stats_sat(centroid, matrix=matrix)
    
    seq_algn = np.zeros((len(centroid), len(alphabet), len(seqs)))
    for seqi, s in enumerate(seqs):
        # a = parasail.nw_trace(centroid, s, open=3, extend=3, matrix=parasail.blosum62)
        a = parasail.nw_trace_scan_profile_sat(centroid_query, s, open=gopen, extend=gextend)
        
        pos = 0
        for ref_aa, q_aa in zip(a.traceback.ref, a.traceback.query):
            if not q_aa == '-':
                seq_algn[pos, alphabet.index(ref_aa), seqi] = 1
                pos += 1
    return seq_algn

def _get_frequencies(seqs, alphabet, weights, add_one=False):
    L = len(seqs[0])
    nAA = len(alphabet)
    
    freq = np.zeros((nAA, L))
    for coli in range(L):
        for si, s in enumerate(seqs):
            freq[alphabet.index(s[coli]), coli] += weights[si]
    if add_one:
        freq = (freq + 1) / (freq + 1).sum(axis=0, keepdims=True)
    else:
        freq = freq / freq.sum(axis=0, keepdims=True)
    return freq

def compute_relative_motif(seqs, refs, alphabet=None):
    """Use statistic related to the Kullback-Liebler divergence to indicate how surprising
    it is to see the AAs in seqs vs. the refs.
    All seqs and refs must have the same length (i.e. be "aligned")
    
    For discussion about relative logos:
    https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2489-3
    https://academic.oup.com/nar/article/40/W1/W281/1076274

    Parameters
    ----------
    seqs : iterable
        List of AA sequences that are all similar to each other and the centroid
    refs : iterable
        List of all other sequences in the experiment as a reference.

    Returns
    -------
    A : pd.DataFrame [AA alphabet x position]
        A matrix of KL scores with a higher divergence
        being more surprising."""

    """
    p_i is reference
    q_i is observed

    A = q * np.log2(q / p)

    """
    """Adding 1 to the refs to avoid inf's, but this should be studied more carefully"""
    if alphabet is None:
        alphabet = aa_alphabet

    p = _get_frequencies(refs, alphabet, np.ones(len(refs)), add_one=True)
    q = _get_frequencies(seqs, alphabet, np.ones(len(seqs)))
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        A = q * np.log2(q/p)
    A[np.isnan(A)] = 0
    A = pd.DataFrame(A, index=alphabet)
    return A

def compute_pal_motif(centroid, seqs, refs=None,
                      gopen=3, gextend=3, matrix=None,
                      ref_freqs=None, alphabet=None, bootstrap_samples=0, alpha=0.05,
                      pseudo_count=1):
    """Compute pairwise alignments between the centroid and all sequences in seqs and refs. The motif
    will have the same length as the centroid with a statistic for each symbol at each position
    indicating the divergence of the observed distribution from the reference. Summing the statistics
    across the alphabet gives the Kullback-Liebler divergence for each position.

    Also returns the position specific log-likliehood of the sequences,
    given the reference uing a multinomial distribution on the alphabet.

    A 95% confidence interval on the output can be optionally computed with a
    non-parametric bootstrap. Recommend 10K samples for good estimation.

    Parameters
    ----------
    centroid : str
        Amino-acid sequence that is also among the seqs
    seqs : list
        List of AA sequences that are all similar to each other and the centroid
    refs : list
        List of all other sequences in the experiment as a reference.
    gopen : int
        Gap open penalty for parasail
    gextend : int
        Gap extend penalty for parasail
    alphabet : str
        String of characters whose frequency should be considered/counted
    matrix : substitution matrix
        Matrix from parasail for the alignment
        None for default parasail.blosum62
    ref_freqs : dict or None
        Reference frequencies for each AA
    pseudo_count : float
        A pseudocount added to each alphabet symbol's count in the reference set,
        to "shrink" the reference probabilities toward uniform and prevent infinite
        outputs.
    bootstrap_samples : bool
        With values > 0, compute percentile bootstrap confidence intervals for
        the position-specific KL and loglikelihood statistics.

    Returns
    -------
    A : pd.DataFrame [AA alphabet x position]
        A matrix of KL divergence scores that can be used directly with the plotPALLogo function
    out : pd.DataFrame [position x (KL, loglikelihood columns)]
        Per position KL divergence from, and log-likelihood of observing the sequences, given the reference.
    """
    if matrix is None:
        matrix = parasail.blosum62
    
    if alphabet is None:
        alphabet = aa_alphabet

    """Returned DataFrame is [positions x alphabet] opposite of output"""
    raw_seq_algn = pairwise_alignment_frequencies(centroid, seqs, gopen=gopen, gextend=gextend, matrix=matrix, alphabet=alphabet)
    seq_algn = pd.DataFrame(raw_seq_algn.sum(axis=2), index=list(centroid), columns=alphabet)
    L = seq_algn.shape[0]

    if refs is None:
        if ref_freqs is None:
            p = np.ones(len(aa_alphabet))
        else:
            p = np.array([ref_freqs.get(aa, 0) for aa in aa_alphabet])
        p = np.tile(p[None, :], (L, 1))
        p = p / np.sum(p, axis=1, keepdims=True)
    else:
        raw_ref_algn = pairwise_alignment_frequencies(centroid, refs, gopen=gopen, gextend=gextend, matrix=matrix)
        ref_algn = pd.DataFrame(raw_ref_algn.sum(axis=2), index=list(centroid), columns=alphabet)
        """Adding 1 to refs to avoid inf's, but this should be studied more carefully
        https://en.wikipedia.org/wiki/Additive_smoothing"""
        p = (ref_algn.values + pseudo_count) / (ref_algn.values + pseudo_count).sum(axis=1, keepdims=True)

    """
    Similar to the Kullback-Liebler divergence between
    the sequences and the reference. Sum the symbols at
    each position to get KL divergence for that position.
    Sum across positions to get overal KL divergence.
    p_i is reference
    q_i is observed

    A = q * np.log2(q / p)
    """
    q = (seq_algn.values) / (seq_algn.values).sum(axis=1, keepdims=True)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        A = q * np.log2(q/p)
        kl = _kl_stat(np.sum(raw_seq_algn, axis=2), p)
    A[np.isnan(A)] = 0
    A = pd.DataFrame(A, index=seq_algn.index, columns=seq_algn.columns)
        
    kl_samps = np.zeros((len(centroid), bootstrap_samples))

    loglik = _multinomial_loglik(np.sum(raw_seq_algn, axis=2), p)
    loglik_samps = np.zeros((len(centroid), bootstrap_samples))
    
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        for sampi in range(bootstrap_samples):
            rind = np.random.randint(raw_seq_algn.shape[2], size=raw_seq_algn.shape[2])
            loglik_samps[:, sampi] = _multinomial_loglik(np.sum(raw_seq_algn[:, :, rind], axis=2), p)
            kl_samps[:, sampi] = _kl_stat(np.sum(raw_seq_algn[:, :, rind], axis=2), p)
    if bootstrap_samples:
        kl_samps.sort(axis=1)
        loglik_samps.sort(axis=1)
        nvals = np.round((bootstrap_samples - 1) * np.array([alpha, 1-alpha])).astype(int)
        out = dict(kl=kl,
                    kl_ll=kl_samps[:, nvals[0]],
                    kl_ul=kl_samps[:, nvals[1]],
                    loglik=loglik,
                    loglik_ll=loglik_samps[:, nvals[0]],
                    loglik_ul=loglik_samps[:, nvals[1]])
    else:
        out = dict(kl=kl, loglik=loglik)
    out = pd.DataFrame(out, index=list(centroid))
    return A.T, out

def _kl_stat(counts, ref_pr):
    q = counts / np.sum(counts, axis=1, keepdims=True)
    return np.nansum(q * np.log2(q / ref_pr), axis=1)

def _multinomial_loglik(counts, ref_pr):
    """A and seq_algn have alphabet in columns (before A gets transposed)"""
    loglik = np.zeros(counts.shape[0])
    for posi in range(counts.shape[0]):
        loglik[posi] = multinomial.logpmf(counts[posi, :], np.sum(counts[posi, :]), ref_pr[posi, :])
    return loglik

def consensus_seq(motif):
    """With a motif pd.DataFrame, the symbol with the highest value at each
    position is the "consensus" sequence."""
    return ''.join(motif.idxmax(axis=0).values)

def multinomial_confint(counts, alpha=0.05):
    """Confidence interval on observed counts from a multinomial distribution.
    Code for the Goodman (1965) method is lifted from statsmodels 0.12
    which was not currently released.

    License for statsmodels package is below. This function will be replaced by
    a dependency on statsmodels 0.12 upon release.

    Copyright (C) 2006, Jonathan E. Taylor
    All rights reserved.

    Copyright (c) 2006-2008 Scipy Developers.
    All rights reserved.

    Copyright (c) 2009-2018 statsmodels Developers.
    All rights reserved.


    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:

      a. Redistributions of source code must retain the above copyright notice,
         this list of conditions and the following disclaimer.
      b. Redistributions in binary form must reproduce the above copyright
         notice, this list of conditions and the following disclaimer in the
         documentation and/or other materials provided with the distribution.
      c. Neither the name of statsmodels nor the names of its contributors
         may be used to endorse or promote products derived from this software
         without specific prior written permission.


    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
    AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
    ARE DISCLAIMED. IN NO EVENT SHALL STATSMODELS OR CONTRIBUTORS BE LIABLE FOR
    ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
    DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
    SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
    CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
    LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
    OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
    DAMAGE."""

    n = counts.sum()
    k = len(counts)
    proportions = counts / n
    chi2 = stats.chi2.ppf(1 - alpha / k, 1)
    delta = chi2 ** 2 + (4 * n * proportions * chi2 * (1 - proportions))
    region = ((2 * n * proportions + chi2 +
               np.array([- np.sqrt(delta), np.sqrt(delta)])) /
              (2 * (chi2 + n))).T
    return region