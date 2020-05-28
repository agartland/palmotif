# palmotif
Python package for computing and plotting sequence logos, with output as SVG or matplotlib

Use PyPI package `seqlogo` for standard logos. This package allows for creating logos of sequences with different lengths.

Credit to developers of more sophisticated [LogoJS](https://github.com/weng-lab/logojs-package) for the SVG alphabet in this package.

Color schemes are derived from those in [weblogo](https://github.com/ostrokach/weblogo)

# Examples
```python
from palmotif import compute_motif, svg_logo
motif = compute_motif(seqs)
svg_logo(motif, 'test.svg', color_scheme='taylor')
```
 - [Protein or DNA sequence logos](https://raw.githubusercontent.com/agartland/palmotif/master/palmotif/tests/test.svg)
 - [Allows negative values](https://raw.githubusercontent.com/agartland/palmotif/master/palmotif/tests/negative.svg)
 - [Shapely color scheme](https://raw.githubusercontent.com/agartland/palmotif/master/palmotif/tests/alphabet.svg)
 - [Taylor color scheme](https://raw.githubusercontent.com/agartland/palmotif/master/palmotif/tests/taylor.svg)
