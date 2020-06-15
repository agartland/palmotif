"""
python -m unittest palmotif/tests/test_compute.py
"""
import unittest
import numpy as np
import pandas as pd
from os.path import join as opj

from palmotif import *

from palmotif.tests.sequences import seqs1, seqs2


class TestCompute(unittest.TestCase):

    def test_compute_motif(self):
        motif = compute_motif(seqs1, reference_freqs=None, weights=None, align_first=False, gap_reduce=None, alphabet=None)
        self.assertTrue(True)

    def test_reference_freqs_motif(self):
        motif = compute_motif(seqs1, reference_freqs=uniprot_frequency, weights=None, align_first=False, gap_reduce=None, alphabet=None)
        self.assertTrue(True)

    def test_compute_motif(self):
        motif = compute_relative_motif(seqs1, seqs2)
        self.assertTrue(True)

    def test_pal_no_ref(self):
        mixed_seqs = seqs1 + [s[:9] for s in seqs1]
        motif = compute_pal_motif(mixed_seqs[0], mixed_seqs, refs=None, gopen=3, gextend=3, matrix=None, ref_freqs=None)

    def test_pal_no_ref(self):
        mixed_seqs1 = seqs1 + [s[:9] for s in seqs1]
        mixed_seqs2 = seqs2 + [s[:9] for s in seqs2]
        motif = compute_pal_motif(mixed_seqs1[1], mixed_seqs1, refs=None, gopen=3, gextend=3, matrix=None, ref_freqs=uniprot_frequency)

    def test_pal_yes_ref(self):
        mixed_seqs1 = seqs1 + [s[:9] for s in seqs1]
        mixed_seqs2 = seqs2 + [s[:9] for s in seqs2]
        motif = compute_pal_motif(mixed_seqs1[2], mixed_seqs1, refs=mixed_seqs2, gopen=3, gextend=3, matrix=None)
    def test_pal_loglik(self):
        mixed_seqs1 = seqs1 + [s[:9] for s in seqs1]
        mixed_seqs2 = seqs2 + [s[:9] for s in seqs2]
        motif, ll = compute_pal_motif(mixed_seqs1[2], mixed_seqs1, refs=mixed_seqs2, gopen=3, gextend=3, matrix=None, return_loglikelihood=True)

    def test_pal_yes_ref_with_liklihood(self):
        mixed_seqs1 = ["A","A","A","A"]
        mixed_seqs2 = ["T","T","T","T"]
        motif,ll = compute_pal_motif(mixed_seqs1[2], mixed_seqs1, refs=mixed_seqs2, gopen=3, gextend=3, matrix=None, return_loglikelihood=True)
        assert ll[0] < 0.000000000000000000000000001 

    def test_pal_yes_ref_with_liklihood2(self):
        """Assert that number is very small in case where seqs match ref when smooth = True, and 0 if smooth = false """
        mixed_seqs1 = ["A","A","A"]
        mixed_seqs2 = ["A","A","A"]
        _,ll = compute_pal_motif(mixed_seqs1[2], mixed_seqs1, refs=mixed_seqs2, gopen=3, gextend=3, matrix=None, return_loglikelihood=True, smooth = True)
        assert ll < 0.0001#,#[np.log(1**3)])#atol=1)
        _,ll = compute_pal_motif(mixed_seqs1[2], mixed_seqs1, refs=mixed_seqs2, gopen=3, gextend=3, matrix=None, return_loglikelihood=True, smooth = False)
        assert ll[0] == np.log(1**3)

    def test_pal_yes_ref_with_liklihood3(self):
        """simple tests that log liklihood is expected in simple bernoulli trial i.e. (.5)^4 if smooth = False and very close if smooth is True"""
        mixed_seqs1 = ["A","A","A","A"]
        mixed_seqs2 = ["T","T","A","A"]
        _,ll = compute_pal_motif(mixed_seqs1[2], mixed_seqs1, refs=mixed_seqs2, gopen=3, gextend=3, matrix=None, return_loglikelihood=True, smooth = True)
        assert np.isclose(ll, [np.log(.5**4)], rtol=0.00001)
        _,ll = compute_pal_motif(mixed_seqs1[2], mixed_seqs1, refs=mixed_seqs2, gopen=3, gextend=3, matrix=None, return_loglikelihood=True, smooth = False)
        assert ll[0] == np.log(.5**4)

    def test_pal_yes_ref_with_liklihood4(self):
        """
        This is an explicit test that laplace smoothing is working such that the log liklihood is not -inf
        even thought T never appears in the refernce set.
        """
        mixed_seqs1 = ["T","A","A"]
        mixed_seqs2 = ["A","A","A"]
        _,ll = compute_pal_motif(mixed_seqs1[2], mixed_seqs1, refs=mixed_seqs2, gopen=3, gextend=3, matrix=None, return_loglikelihood=True, smooth = True)
        assert ll[0]> np.NINF
  
    def test_pal_yes_ref_with_liklihood_no_smoothing(self):
        """
        With no smoothing, get -inf when amino in seqs not in ref seqs
        """
        mixed_seqs1 = ["T","A","A"]
        mixed_seqs2 = ["A","A","A"]
        _,ll = compute_pal_motif(mixed_seqs1[2], mixed_seqs1, refs=mixed_seqs2, gopen=3, gextend=3, matrix=None, return_loglikelihood=True, smooth = False)
        assert ll[0] == np.NINF

if __name__ == '__main__':
    unittest.main()
