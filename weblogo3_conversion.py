import numpy as np
import pandas as pd
import weblogo
from weblogo import ColorScheme, SymbolColor

import IPython.display

seqs = ['CACADLGAYPDKLIF',
         'CACDALLAYTDKLIF',
         'CACDAVGDTLDKLIF',
         'CACDDVTEVEGDKLIF',
         'CACDFISPSNWGIQSGRNTDKLIF',
         'CACDILLGDTADKLIF',
         'CACDIVLSGGLDTRQMFF',
         'CACDLLLRQSSTDKLIF',
         'CACDNLSETTDKLIF',
         'CACDPLGTDKLIF',
         'CACDPMGGSGGLSWDTRQMFF',
         'CACDPVLGDTRLTDKLIF',
         'CACDPVQGYSGQNRAYTDKLIF',
         'CACDSILGDTLYTDKLIF',
         'CACDSLTSHTGGFGPDKLIF',
         'CACDSTGDLSSWDTRQMFF',
         'CACDSVESRNVLGDPTTDKLIF',
         'CACDSVLSRDLGDSELIF',
         'CACDTAAGGYASSWDTRQMFF',
         'CACDTAPHGGRTWDTRQMFF',
         'CACDTGGYVNWDTRQMFF',
         'CACDTGRLLGDTADTRQMFF',
         'CACDTIRGFSSWDTRQMFF',
         'CACDTIVAPALDKLIF',
         'CACDTLFLGEDTPTDKLIF']

sys.path.append(opj(_git, 'palmotif'))
import palmotif

motif, s = palmotif.compute_pal_motif(seqs[0], seqs)

"""Weblogo3 does not support negative values!"""
# motif = motif * np.random.choice([-1, 1], size=motif.shape)

#fseqs = [s[:5] for s in seqs]
#logoseqs = weblogo.seq.SeqList(alist=fseqs, alphabet=weblogo.seq.Alphabet(''.join(motif.index)))
#logodata = weblogo.LogoData().from_seqs(logoseqs)

"""motif is a pd.DataFrame with logo positions across the columns,
a row index containing the alphabet symbols, and each value containing
an arbitrary float value for symbol height"""
alphabet = weblogo.seq.Alphabet(''.join(motif.index))
hydrophobicity = ColorScheme([SymbolColor("RKDENQ", "blue", "hydrophilic"),
                              SymbolColor("SGHTAP", "green", "neutral"),
                              SymbolColor("YVMCLFIW", "black", "hydrophobic"),
                              SymbolColor("-*ZX", "gray", "default")],
                             alphabet=alphabet)

logodata = weblogo.LogoData(length=motif.shape[1],
                            alphabet=alphabet,
                            entropy=motif.sum(axis=0),
                            weight=np.ones(motif.shape[1]),
                            counts=motif.values.T)

logooptions = weblogo.LogoOptions(  resolution=200,
                                    show_fineprint=False,
                                    yaxis_scale=1.1*np.max(motif.values),
                                    color_scheme=hydrophobicity,
                                    unit_name='nats',
                                    yaxis_label='bits')
logoformat = weblogo.LogoFormat(logodata, logooptions)
logobytes = weblogo.logo_formatter.png_formatter(logodata, logoformat)
IPython.display.Image(data=logobytes)