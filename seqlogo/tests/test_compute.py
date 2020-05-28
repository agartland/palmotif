"""
python -m unittest seqlogo/tests/test_compute.py
"""
import unittest
import numpy as np
import pandas as pd
from os.path import join as opj

from seqlogo import *

from seqlogo.tests.sequences import seqs1, seqs2


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


   

if __name__ == '__main__':
    unittest.main()
