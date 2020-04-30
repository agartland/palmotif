import pandas as pd
import numpy as np
import warnings

try:
    import parasail
except ImportError:
    print('Importing seqlogo without alignment support from parasail.')

from .muscle import muscle_align

__all__ = ['compute_motif',
           'compute_pal_motif',
           'uniprot_frequency',
           'compute_relative_motif']

aa_alphabet = [aa for aa in '-ABCDEFGHIKLMNPQRSTVWXYZ']

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
    """Compute heights for a sequence logo

    Parameters
    ----------
    seqs : list or pd.Series
        Alignment of AA sequences
    reference_freqs : dict or None
        Reference frequencies for each AA
    weights : np.array
        Weights for each sequence
    align_first : bool
        Attempt a multiple sequence alignment using MUSCLE.
        (must be installed and on the path)
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
    
    if align_first:
        align = muscle_align(seqs)
    else:
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
    removes the column from the alignment"""

    removePos = []
    for pos in range(len(align.iloc[0])):
        if np.mean(align.map(lambda seq: seq[pos] == '-')) > thresh:
            removePos.append(pos)
    for pos in removePos:
        align = align.loc[align.map(lambda seq: seq[pos] == '-')]
    return align.map(lambda seq: ''.join([aa for pos, aa in enumerate(seq) if not pos in removePos]))

def pairwise_alignment_frequencies(centroid, seqs, gopen=3, gextend=3, matrix=parasail.blosum62):
    alphabet = sorted([aa for aa in skbio.sequence.Protein.alphabet if not aa in '*.'])
    centroid_query = parasail.profile_create_stats_sat(centroid, matrix=matrix)
    
    seq_algn = np.zeros((len(centroid), len(alphabet)))
    for s in seqs:
        # a = parasail.nw_trace(centroid, s, open=3, extend=3, matrix=parasail.blosum62)
        a = parasail.nw_trace_scan_profile_sat(centroid_query, s, open=gopen, extend=gextend)
        
        pos = 0
        for ref_aa, q_aa in zip(a.traceback.ref, a.traceback.query):
            if not q_aa == '-':
                seq_algn[pos, alphabet.index(ref_aa)] = seq_algn[pos, alphabet.index(ref_aa)] + 1
                pos += 1
    return pd.DataFrame(seq_algn, index=list(centroid), columns=alphabet)

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

def compute_relative_motif(seqs, refs):
    """Use log-OR scores indicating how likely it was to see the AA
    in the seqs vs. the refs. Seqs and refs must have equal length.

    Parameters
    ----------
    seqs : list
        List of AA sequences that are all similar to each other and the centroid
    refs : list
        List of all other sequences in the experiment as a reference.

    Returns
    -------
    A : pd.DataFrame [AA alphabet x position]
        A matrix of log-OR scores"""

    alphabet = aa_alphabet

    """
    p_i is reference
    q_i is observed

    A = q * np.log2(q / p) - c(N)

    """
    """Adding 1 to avoid inf's, but this should be studied more carefully"""
    p = _get_frequencies(refs, alphabet, np.ones(len(refs)), add_one=True)
    q = _get_frequencies(seqs, alphabet, np.ones(len(seqs)))
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        A = q * np.log2(q/p)
    A[np.isnan(A)] = 0
    A = pd.DataFrame(A, index=alphabet)
    return A

def compute_pal_motif(centroid, seqs, refs, gopen=3, gextend=3, matrix=None):
    """Compute pairwise alignments between the centroid and all sequences in seqs and refs. The motif
    will have the same length as the centroid with log-OR scores indicating how likely it was to see the AA
    in the seqs vs. the refs.

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
    matrix : substitution matrix
        Matrix from parasail for the alignment

    Returns
    -------
    A : pd.DataFrame [AA alphabet x position]
        A matrix of log-OR scores that can be used directly with the plotPALLogo function"""
    if matrix is None:
        matrix = parasail.blosum62

    seq_algn = pairwise_alignment_frequencies(centroid, seqs, gopen=gopen, gextend=gextend, matrix=matrix)
    ref_algn = pairwise_alignment_frequencies(centroid, refs, gopen=gopen, gextend=gextend, matrix=matrix)

    """
    p_i is reference
    q_i is observed

    A = q * np.log2(q / p) - c(N)

    """
    """Adding 1 to avoid inf's, but this should be studied more carefully"""
    p = (ref_algn.values + 1) / (ref_algn.values + 1).sum(axis=1, keepdims=True)
    q = (seq_algn.values) / (seq_algn.values).sum(axis=1, keepdims=True)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        A = q * np.log2(q/p)
    A[np.isnan(A)] = 0
    A = pd.DataFrame(A, index=ref_algn.index, columns=ref_algn.columns)
    #pdf = pd.DataFrame(p, index=ref_algn.index, columns=ref_algn.columns)
    #qdf = pd.DataFrame(q, index=ref_algn.index, columns=ref_algn.columns)
    return A.T