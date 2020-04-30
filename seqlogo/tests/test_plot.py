"""
python -m unittest seqlogo/tests/test_plot.py
"""
import unittest
import numpy as np
import pandas as pd
from os.path import join as opj

from seqlogo import *

from seqlogo.tests.sequences import seqs1, seqs2


class TestSeqlogo(unittest.TestCase):
    """TODO: generate these test files and then run tests against them"""
    
    def test_svg_plot(self):
        fn = opj('seqlogo', 'tests', 'test.svg')
        a = pd.DataFrame(np.array([[0.1, 0.2, 0, 0.5],[0, 0.3, 0.8, 0.1], [0.5, 0.1, 0.4, 0.]]).T, columns=[0, 1, 2], index=['A', 'C', 'T', 'G'])
        svg_logo(a, fn, color_scheme='nucleotide')
        self.assertTrue(True)
    
    def test_svg_negative(self):
        fn = opj('seqlogo', 'tests', 'negative.svg')
        a = pd.DataFrame(np.array([[0.1, -0.2, 0, 0.5],[0, 0.3, -0.8, 0.1]]).T, columns=[0, 1], index=['C', 'R', 'W', 'S'])
        svg_logo(a, fn)
        self.assertTrue(True)
    
    def test_svg_alphabet(self):
        fn = opj('seqlogo', 'tests', 'alphabet.svg')
        aas = [aa for aa in 'ABCDEFGHIKLMNPQRSTVWXYZ']

        np.random.seed(110820)
        a = pd.DataFrame(np.log10(10*np.random.rand(len(aas), 10)), columns=['P%d' % (i + 1) for i in range(10)], index=aas)
        svg_logo(a, fn)
        self.assertTrue(True)

    def test_svg_taylor(self):
        fn = opj('seqlogo', 'tests', 'taylor.svg')
        aas = [aa for aa in 'ABCDEFGHIKLMNPQRSTVWXYZ']

        np.random.seed(110820)
        a = pd.DataFrame(np.log10(10*np.random.rand(len(aas), 10)), columns=['P%d' % (i + 1) for i in range(10)], index=aas)
        svg_logo(a, fn, color_scheme='taylor')
        self.assertTrue(True)

    def test_compute_motif(self):
        motif = compute_motif(seqs1, reference_freqs=None, weights=None, align_first=False, gap_reduce=None, alphabet=None)
        fn = opj('seqlogo', 'tests', 'compute.svg')
        svg_logo(motif, fn)
        self.assertTrue(True)

    def test_reference_freqs_motif(self):
        motif = compute_motif(seqs1, reference_freqs=uniprot_frequency, weights=None, align_first=False, gap_reduce=None, alphabet=None)
        fn = opj('seqlogo', 'tests', 'ref_freq.svg')
        svg_logo(motif, fn)
        self.assertTrue(True)

    def test_compute_motif(self):
        motif = compute_relative_motif(seqs1, seqs2)
        fn = opj('seqlogo', 'tests', 'relative_motif.svg')
        svg_logo(motif, fn)
        self.assertTrue(True)
   

if __name__ == '__main__':
    unittest.main()
