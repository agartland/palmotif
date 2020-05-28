"""Color schemes derived from those in weblogo:
https://github.com/ostrokach/weblogo/blob/master/weblogolib/colorscheme.py"""

#  Copyright (c) 2003-2005 The Regents of the University of California.
#  Copyright (c) 2005 Gavin E. Crooks

#  This software is distributed under the MIT Open Source License.
#  <http://www.opensource.org/licenses/mit-license.html>
#
#  Permission is hereby granted, free of charge, to any person obtaining a 
#  copy of this software and associated documentation files (the "Software"),
#  to deal in the Software without restriction, including without limitation
#  the rights to use, copy, modify, merge, publish, distribute, sublicense,
#  and/or sell copies of the Software, and to permit persons to whom the
#  Software is furnished to do so, subject to the following conditions:
#
#  The above copyright notice and this permission notice shall be included
#  in all copies or substantial portions of the Software.
#
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
#  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN 
#  THE SOFTWARE.

""" Popular color codings for nucleic and amino acids. 
Classes:
    ColorScheme -- A color scheme
    ColorGroup  
    
    
Generic
    monochrome
Nucleotides
    nucleotide
    base pairing
Amino Acid
    hydrophobicity
    chemistry
    charge
    taylor
Status : Beta - Needs documentation.
"""
# Good online references include bioruby and the JalView alignment editor.
# Clamp, M., Cuff, J., Searle, S. M. and Barton, G. J. (2004), 
# "The Jalview Java Alignment Editor," Bioinformatics, 12, 426-7
# http://www.jalview.org

import pandas as pd
import io

__all__ = ['aa_grouping',
           'color_group',
           'nucleotide',
           'base_pairing',
           'hydrophobicity',
           'chemistry',
           'charge',
           'taylor',
           'logojs',
           'shapely',
           'black']

class aa_grouping(object):
    """Builds lookup dicts of color and label assignments from the set of color_groups"""
    def __init__(self, groups, title='', description=''):
        self.groups = groups
        self.labelDict = {}
        self.colorDict = {}
        for g in groups:
            self.colorDict.update({aa:g.color for aa in g.aas})
            self.labelDict.update({aa:g.label for aa in g.aas})

    def label(self, aa):
        return self.labelDict.get(aa, '')

    def color(self, aa):
        return self.colorDict.get(aa, 'black')

class color_group(dict):
    """Builds a lookup dict for the amino-acids specified"""
    def __init__(self, aas, color, label=''):
        dict.__init__(self)
        self.color = color
        self.label = label
        self.aas = aas

black = aa_grouping([])

nucleotide = aa_grouping([
    color_group("G", "orange"),
    color_group("TU", "red"),
    color_group("C",  "blue"),
    color_group("A",  "green")]) 

base_pairing = aa_grouping([
    color_group("TAU",  "darkorange", "Weak (2 Watson-Crick hydrogen bonds)"),
    color_group("GC",    "blue", "Strong (3 Watson-Crick hydrogen bonds)")])

# From Crooks2004c-Proteins-SeqStr.pdf
hydrophobicity = aa_grouping([
    color_group( "RKDENQ",   "blue", "hydrophilic"),
    color_group( "SGHTAP",   "green", "neutral"  ),
    color_group( "YVMCLFIW", "black",  "hydrophobic") ])

# from makelogo
chemistry = aa_grouping([
  color_group( "GSTYC",  "green",   "polar"),
  color_group( "NQ",      "purple", "neutral"), 
  color_group( "KRH",     "blue",   "basic"),
  color_group( "DE",      "red",    "acidic"),
  color_group("PAWFLIMV", "black",  "hydrophobic") ])   

charge = aa_grouping([
    color_group("KRH", "blue", "Positive" ),
    color_group( "DE", "red", "Negative") ])


taylor = aa_grouping([
    color_group( 'A', '#CCFF00' ),
    color_group( 'C', '#FFFF00' ),
    color_group( 'D', '#FF0000'),
    color_group( 'E', '#FF0066' ),
    color_group( 'F', '#00FF66'),
    color_group( 'G', '#FF9900'),
    color_group( 'H', '#0066FF'),
    color_group( 'I', '#66FF00'),
    color_group( 'K', '#6600FF'),
    color_group( 'L', '#33FF00'),
    color_group( 'M', '#00FF00'),
    color_group( 'N', '#CC00FF'),
    color_group( 'P', '#FFCC00'),
    color_group( 'Q', '#FF00CC'),
    color_group( 'R', '#0000FF'),
    color_group( 'S', '#FF3300'),
    color_group( 'T', '#FF6600'),
    color_group( 'V', '#99FF00'),
    color_group( 'W', '#00CCFF'),
    color_group( 'Y', '#00FFCC')],
    title = "Taylor",
    description = "W. Taylor, Protein Engineering, Vol 10 , 743-746 (1997)")

_logojs_colors = {'A': '#000000',
                 'B': '#bb8800',
                 'C': '#008811',
                 'D': '#ff0000',
                 'E': '#ff0022',
                 'F': '#333333',
                 'G': '#007700',
                 'H': '#220099',
                 'I': '#111111',
                 'K': '#0000aa',
                 'L': '#002222',
                 'M': '#220022',
                 'N': '#009911',
                 'P': '#080808',
                 'Q': '#00aa00',
                 'R': '#0022aa',
                 'S': '#008f00',
                 'T': '#006600',
                 'V': '#222200',
                 'W': '#080808',
                 'Y': '#00a800',
                 'Z': '#aaaa00',
                 '-': '#000000'}
"""Amino acids colored according to chemical properties (acidic, basic, and non-polar are red, blue, and black, respectively"""
logojs = aa_grouping([color_group(aa, color) for aa, color in _logojs_colors.items()],
                     title='LogoJS',
                     description='https://github.com/weng-lab/logojs-package')


_shapely_colors = 'Amino Acids,Color Name,RGB Values,Hexadecimal\nD|E,bright red,"[230,10,10]",E60A0A\nC|M,yellow,"[230,230,0]",E6E600\nK|R,blue,"[20,90,255]",145AFF\nS|T,orange,"[250,150,0]",FA9600\nF|Y,mid blue,"[50,50,170]",3232AA\nN|Q,cyan,"[0,220,220]",00DCDC\nG,light grey,"[235,235,235]",EBEBEB\nL|V|I,green,"[15,130,15]",0F820F\nA,dark grey,"[200,200,200]",C8C8C8\nW,pink,"[180,90,180]",B45AB4\nH,pale blue,"[130,130,210]",8282D2\nP,flesh,"[220,150,130]",DC9682'
_colors_df = pd.read_csv(io.StringIO(_shapely_colors), delimiter=',')
shapely = aa_grouping([color_group(r['Amino Acids'].replace('|', ''), '#%s' % r['Hexadecimal']) for i,r in _colors_df.iterrows()],
                      title='Shapely',
                      description='http://acces.ens-lyon.fr/biotic/rastop/help/colour.htm#shapelycolours')
