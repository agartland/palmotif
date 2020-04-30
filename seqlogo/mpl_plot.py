import matplotlib.pyplot as plt
import matplotlib.patheffects
import matplotlib.transforms as mtrans
from matplotlib import transforms
import numpy as np
from . import aacolors

__all__ = ['mpl_logo']

class Scale(matplotlib.patheffects.RendererBase):
    def __init__(self, sx=1., sy=1.):
        self._sx = sx
        self._sy = sy

    def draw_path(self, renderer, gc, tpath, affine, rgbFace):
        affine = affine.identity().scale(self._sx, self._sy) + affine
        renderer.draw_path(gc, tpath, affine, rgbFace)

def _extend_bbox(bbox, ext):
    if ext.x0 < bbox[0, 0]:
        bbox[0, 0] = ext.x0
    if ext.y0 < bbox[0, 1]:
        bbox[0, 1] = ext.y0
    if ext.x1 > bbox[1, 0]:
        bbox[1, 0] = ext.x1
    if ext.y1 < bbox[1, 1]:
        bbox[1, 1] = ext.y1
    return bbox

def mpl_logo(x, y, motif, axh=None, fontsize=16, color_scheme='shapely'):
    """Sequence logo of the motif at data coordinates x,y using matplotlib.
    
    Not a great solution since logos are difficult to place in coordinate space and
    may be distorted by zooming after it is plotted.

    Inspiration:
    https://github.com/saketkc/motif-logos-matplotlib/blob/master/Motif%20Logos%20using%20matplotlib.ipynb

    Parameters
    ----------
    x,y : float
        Position for the bottom-left corner of the logo in the axes data coordinates
    motif : pd.DataFrame
        Matrix of scores for each AA symbol (index, rows) and each
        column in an alignment (columns). Values will be used to scale
        each symbol linearly.
    axh : matplotlib axes handle
        Will use plt.gca() if None
    fontsize : float
        Pointsize of font passed to axh.text
    color_scheme : str
        An aa_grouping instance from aacolors.py or custom define.
        Options: nucleotide, base_pairing, hydrophobicity, chemistry, charge, taylor, logojs, shapely

    Returns
    -------
    letters : list of lists of tuples (Text, height)
        List of lists of the Text objects in the order they are plotted
        (lower left to upper right). Useful for putting annotations relative
        to specific objects in the motif.
    bbox : [[x0, y0], [x1, y1]]
        Full extent of the logo in screen coodinates.
    bbox_data : [[x0, y0], [x1, y1]]
        Full extent of the logo in data coodinates."""
    colors = getattr(aacolors, color_scheme)
    if axh is None:
        axh = plt.gca()
  
    trans_offset = transforms.offset_copy(axh.transData, 
                                      fig=axh.figure, 
                                      x=0, y=0, 
                                      units='dots')
    bottom_offset = trans_offset
    neg_trans_offset = trans_offset
    first_pos = None
    last_neg = None
    bbox = np.array([[x, y], [x, y]])
    mxy = y 
    mny = y
    letters = [[] for i in range(motif.shape[1])]
    for xi in range(motif.shape[1]):
        scores = motif.iloc[:, xi]
        pos_scores = scores[scores > 0].sort_values()
        neg_scores = (-scores[scores < 0]).sort_values()
        yshift = 0
        neg_trans_offset = trans_offset
        for yi, (aa, score) in enumerate(pos_scores.items()):
            txt = axh.text(x, 
                          y, 
                          aa, 
                          transform=trans_offset,
                          fontsize=fontsize, 
                          color=colors.color(aa),
                          ha='left',
                          family='monospace')
            txt.set_path_effects([Scale(1.5, score)])
            axh.figure.canvas.draw()
            window_ext = txt.get_window_extent(txt._renderer)
            
            #bbox_data = axh.transData.inverted().transform(np.array([[window_ext.x0, window_ext.y0], [window_ext.x1, window_ext.y1]]))
            #bbox_data = window_ext.inverse_transformed(axh.transData)
            #plt.plot(bbox_data.corners()[:, 0], bbox_data.corners()[:, 1], '-k')

            bbox = _extend_bbox(bbox, window_ext)
            letters[xi].append((txt, window_ext.height * score))
            if first_pos is None:
                first_pos = window_ext

            if yshift == 0:
                bottom_offset = transforms.offset_copy(txt._transform, 
                                                      fig=axh.figure,
                                                      x=window_ext.width * 1.5,
                                                      units='dots')

            yshift += window_ext.height * score
            
            # print(xi, aa, window_ext.x0, window_ext.y0, '%1.1f' % score, '%1.0f' % window_ext.height, '%1.0f' % window_ext.width, '%1.1f' % yshift)
            trans_offset = transforms.offset_copy(txt._transform, 
                                                  fig=axh.figure,
                                                  y=window_ext.height * score,
                                                  units='dots')
        trans_offset = bottom_offset
        if yshift > mxy:
            mxy = yshift

    trans_offset = transforms.offset_copy(axh.transData, 
                                          fig=axh.figure, 
                                          x=0, y=0, 
                                          units='dots')    
    for xi in range(motif.shape[1]):
        scores = motif.iloc[:, xi]
        pos_scores = scores[scores > 0].sort_values()
        neg_scores = (-scores[scores < 0]).sort_values()
        yshift = 0
        for yi, (aa, score) in enumerate(neg_scores.items()):
            txt = axh.text(x, 
                          y, 
                          aa, 
                          transform=trans_offset,
                          fontsize=fontsize, 
                          color=colors.color(aa),
                          ha='left',
                          va='top',
                          family='monospace')
            txt.set_path_effects([Scale(1.5, -score)])
            axh.figure.canvas.draw()
            window_ext = txt.get_window_extent(txt._renderer)
            bbox = _extend_bbox(bbox, window_ext)
            letters[xi].append((txt, window_ext.height * score))
            last_neg = window_ext

            if yshift == 0:
                bottom_offset = transforms.offset_copy(txt._transform, 
                                                      fig=axh.figure,
                                                      x=window_ext.width * 1.5,
                                                      units='dots')
            yshift -= window_ext.height * score
            
            # print(xi, aa, '%1.1f' % score, '%1.0f' % window_ext.height, '%1.0f' % window_ext.width, '%1.1f' % yshift)
            trans_offset = transforms.offset_copy(txt._transform, 
                                                      fig=axh.figure,
                                                      y=-window_ext.height * score,
                                                      units='dots')
        trans_offset = bottom_offset
        if yshift < mny:
            mny = yshift
        

    if last_neg is None:
        last_neg = first_pos
    plt.show()
    
    #bbox = np.array([[first_pos.x0, last_neg.y0], [window_ext.x1, window_ext.y0 + window_ext.height * score]])
    bbox_data = axh.transData.inverted().transform(bbox)
    return letters, bbox, bbox_data

