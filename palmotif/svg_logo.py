import numpy as np
import pandas as pd
import svgwrite
import matplotlib as mpl

from .svg_alphabet import add_letter
from . import aacolors

__all__ = ['svg_logo']

"""Note for converting SVG to PNG with imagemagick:
import subprocess

cmd = ['convert', '-density 200', 'test.svg', 'test.png']
subprocess.call(' '.join(cmd), shell=True)"""

_basic_html = """<html>
<svg width="{width}" height="{height}">
{svg}
</svg>
</html>"""

def svg_logo(motif, filename=None,
             color_scheme='chemistry',
             return_str=False, return_html=False,
             yscale=1, svg_height='500px', svg_width='500px',
             ylabel='',
             yrange=None):
    """Sequence logo of the motif using SVG.

    Parameters
    ----------
    motif : pd.DataFrame
        Matrix of scores for each AA symbol (index, rows) and each
        column in an alignment (columns). Values will be used to scale
        each symbol linearly.
    filename : path
        File for saving the SVG.
    color_scheme : str
        An aa_grouping instance from aacolors.py or custom define.
        Options: nucleotide, base_pairing, hydrophobicity, chemistry, charge, taylor, logojs, shapely
    return_str : bool
        Optionally return SVG as text instead of saving to filename."""
    if filename is None:
        return_str = True
    if return_html:
        return_str = True

    colors = getattr(aacolors, color_scheme)

    yscale *= 4
    margin = 10
    left_margin = 100
    bottom_margin = 50
    xpad = 4
    HW = 100
    fontsize='20pt'
    ylabel_fontsize='28pt'

    if yrange is None:
        """Scale height of 100px to absolute max value"""
        # mx_value = np.max(np.abs(motif.values))
        mx_value = np.max(np.sum(np.abs(motif.values), axis=0))
        hscale = HW / mx_value
        wscale = 1

        mx_pos = 0
        mn_neg = 0
        for j in range(motif.shape[1]):
            tmp = motif.values[:, j]
            tot = np.sum(tmp[tmp>0])
            if tot > mx_pos:
                mx_pos = tot
            tot = np.sum(tmp[tmp<0])
            if tot < mn_neg:
                mn_neg = tot
    else:
        mn_neg, mx_pos = yrange
        mx_value = mx_pos - mn_neg
        hscale = HW / mx_value
        wscale = 1
        
    
    yticklabels = mpl.ticker.MaxNLocator(nbins=5, steps=[1, 2, 2.5, 5, 10]).tick_values(mn_neg, mx_pos)
    yticks = [hscale * yt * yscale for yt in yticklabels]

    mx_pos = hscale * mx_pos
    mn_neg = hscale * mn_neg

    mx_pos = np.max([mx_pos, yticks[-1]])
    mn_neg = np.min([mn_neg, yticks[0]])

    yzero = margin + mx_pos
    xzero = left_margin

    height = mx_pos - mn_neg + margin + bottom_margin
    width = wscale * (HW + xpad) * motif.shape[1] + margin + left_margin + xpad

    xticks = [xzero + (i + 1) * wscale * (HW + xpad) - wscale * HW / 2 for i in range(motif.shape[1])]
    # xticklabels = ['%d' % (i + 1) for i in range(motif.shape[1])]
    xticklabels = motif.columns

    # dwg = svgwrite.Drawing(filename=filename, height=height, width=width)
    # dwg = svgwrite.Drawing(filename=filename, size=(width, height), viewBox='0 0 %f %f' % (width, height))
    dwg = svgwrite.Drawing(filename=filename, size=('100%', '100%'), viewBox='0 0 %f %f' % (width, height))

    letter_groups = {}
    for xi in range(motif.shape[1]):
        xshift = xzero + xi*(HW + xpad) + xpad
        scores = motif.iloc[:, xi]
        pos_scores = scores[scores > 0].sort_values()
        neg_scores = (-scores[scores < 0]).sort_values()
        
        posshift = 2
        for yi, (aa, score) in enumerate(pos_scores.items()):
            scaled_height = hscale * score * yscale
            translate = (xshift, yzero - posshift - scaled_height)
            transform = 'translate({xtrans} {ytrans}) scale({xscale} {yscale})'.format(xtrans=translate[0], ytrans=translate[1],
                                                                                        xscale=wscale, yscale=yscale * score/mx_value)
            letter_groups[(xi, aa)] = add_letter(dwg, aa, group_id='%d_%s' % (xi, aa), color=colors.color(aa), background='white', transform=transform)
            #box = dwg.add(dwg.rect(insert=(x*cm, y*cm), size=(X*cm, score*fontsize*cm), fill='green', opacity=score))
            #print(aa, x, y, score)
            posshift += scaled_height

        negshift = 2
        for yi, (aa, score) in enumerate(neg_scores.items()):
            scaled_height = hscale * score * yscale
            translate = (xshift, yzero + negshift)
            transform = 'translate({xtrans} {ytrans}) scale({xscale} {yscale})'.format(xtrans=translate[0], ytrans=translate[1],
                                                                                        xscale=wscale, yscale=yscale * score/mx_value)
            letter_groups[(xi, aa)] = add_letter(dwg, aa, group_id='%d_%s' % (xi, aa), color=colors.color(aa), background='white', transform=transform)
            #box = dwg.add(dwg.rect(insert=(x*cm, y*cm), size=(X*cm, score*fontsize*cm), fill='green', opacity=score))
            #print(aa, x, y, score)
            negshift += scaled_height
    
    axes = dwg.add(dwg.g(id='axes', stroke_width=1.5, stroke='#000000'))
    axes.add(dwg.path(d='M {xz} {yz} h {len}'.format(xz=xzero, yz=yzero, len=wscale * (HW + xpad) * motif.shape[1] + 2 * xpad)))
    axes.add(dwg.path(d='M {xz} {yz} v {len}'.format(xz=xzero, yz=yzero - np.min(yticks), len=-(yticks[-1] - yticks[0]))))
    for yt, ytl in zip(yticks, yticklabels):
        axes.add(dwg.path(d='M {xz} {yz} h {len}'.format(xz=xzero, yz=yzero - yt, len=-(margin / 2))))
        axes.add(dwg.text(ytl, (xzero - (margin/2) - xpad, yzero - yt),
                          fill='#000000',
                          font_size=fontsize,
                          font_family='sans-serif',
                          font_weight='normal',
                          text_anchor='end',
                          dominant_baseline='middle'))
    for xt, xtl in zip(xticks, xticklabels):
        axes.add(dwg.text(xtl, (xt, yzero - yticks[0] + 10),
                          fill='#000000',
                          font_size=fontsize,
                          font_family='sans-serif',
                          font_weight='normal',
                          text_anchor='middle',
                          dominant_baseline='hanging'))
    x, y = (0, yzero - yticks[len(yticks) // 2])
    axes.add(dwg.text(ylabel, (x, y),
                          fill='#000000',
                          font_size=ylabel_fontsize,
                          font_family='sans-serif',
                          font_weight='normal',
                          text_anchor='middle',
                          dominant_baseline='hanging',
                          transform=f"rotate(-90, {x}, {y})"))

    if not return_str:
        dwg.save()
        return filename
    else:
        svg = dwg.tostring()
        if return_html:
            return _basic_html.format(svg=svg, height=svg_height, width=svg_width)
        else:
            return svg

