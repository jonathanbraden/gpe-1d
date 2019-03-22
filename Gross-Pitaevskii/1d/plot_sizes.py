#!/usr/bin/env python
import numpy as np

fig_wid = 5.93 # change based on journal
fig_half_wid = 2.97
fig_single_wid = 0.75*fig_wid

# Basic sizes of stuff
tick_size = 2.
tick_pad = 2.5
tick_font = 9.
label_pad = 4.
label_font = 11.
box_pad = 5.
cb_wid = 1.5*label_font
cb_size = label_font + tick_size + 3.*tick_font + tick_pad + label_pad + cb_wid  # unscale to axis size here

pnt_to_in = 1./72.  # Converts points to inches
gr = 0.5*(1.+np.sqrt(5.)); gr_i = 1./gr  # Golden ratio and inverse

# Axis padding in inches
left_margin_in = pnt_to_in*(box_pad + label_font + label_pad + 3.*tick_font + tick_pad + tick_size)
top_margin_in = pnt_to_in*(box_pad + label_font + label_pad)
bottom_margin_in = pnt_to_in*(label_font + label_pad + tick_size + tick_pad + tick_font + box_pad)
right_margin_in = pnt_to_in*(box_pad + cb_size)  # Change this

def add_cb_above(_f,frac,pad,nc=1,nr=1,wspace=0.04,hspace=0.04):
    f_w, f_h = _f.get_size_inches()
    t_c = _f.subplotpars.top
    b_c = _f.subplotpars.bottom
    ax_h = f_h*(t_c-b_c)/(nr+(nr-1)*hspace)
    
    cb_v = (frac+pad)*ax_h
    cb_v_lab = pnt_to_in*(label_font+label_pad+tick_size+tick_pad+tick_font+box_pad)-top_margin_in
    f_h_ = f_h+cb_v+cb_v_lab

    l_c = _f.subplotpars.left; r_c = _f.subplotpars.right

    t_n = t_c*(f_h/f_h_); b_n = b_c*(f_h/f_h_)
    _f.set_size_inches([f_w,f_h_])
    _f.subplots_adjust(bottom=b_n,top=t_n)
    cax = _f.add_axes([l_c+0.125*(r_c-l_c),t_n + pad*(ax_h/f_h_)+pnt_to_in*(tick_pad+tick_font+tick_size)/f_h_,0.75*(r_c-l_c),frac*ax_h/f_h_])
    
    # Now, need to put the axis label on top, shrink to be some fraction of existing axes
    return cax

# Do to, adjust this so the axis is centered
def compute_fig_size(f_w,ax_w,nr=1,nc=1,wspace=0.04,hspace=0.02,cb_=False,cb_top=False,frac=0.1,pad=0.05):
    ax_h = ax_w*gr_i
    h_size = nc*ax_w*(1.+(nc-1)*wspace) + left_margin_in    # Put in wspace as well
    v_size = nr*ax_h*(1.+(nr-1)*hspace) + bottom_margin_in  # Put in hspace as well
    cb_v = 0.; cb_h = 0.
    cb_v_lab = 0.; cb_h_lab=0.
    if cb_:
        if cb_top:
            cb_v = (frac+pad)*ax_h
            cb_v_lab = pnt_to_in*(label_font+label_pad+tick_size+tick_pad+tick_font)
        else:
            cb_h = (frac+pad)*ax_w
    
    f_h = v_size + cb_v + cb_v_lab + top_margin_in  # Need to add space for colorbar labels
    
    l_frac = left_margin_in / f_w
    b_frac = bottom_margin_in / f_h
    r_frac = h_size / f_w
    t_frac = v_size / f_h # 1. - top_margin_in / f_h
    return (f_w,f_h), (l_frac,b_frac,r_frac,t_frac)

def adjust_figure_margins(f,ax_w,nc=1,nr=1,wspace=0.02,hspace=0.02):
    f_w = f.get_figwidth()

    ax_h = ax_w*gr_i
    h_size = nc*ax_w*(1.+(nc-1)*wspace) + left_margin_in    # Put in wspace as well
    v_size = nr*ax_h*(1.+(nr-1)*hspace) + bottom_margin_in  # put in hspace as well

    f_h = v_size + top_margin_in # Fix if I include a colorbar
    
    l_frac = left_margin_in / f_w
    b_frac = bottom_margin_in / f_h
    r_frac = h_size / f_w
    t_frac = 1. - top_margin_in / f_h

    f.set_size_inches([f_w,f_h])
    f.subplots_adjust(left=l_frac,bottom=b_frac,right=r_frac,top=t_frac)
    return

def adjust_figure_colorbar_top(f,frac,pad):
    """
    Readjust size of figure once we've added a horzontal colorbar on top
    """
    f_w,f_h = f.get_size_inches()
    t_sz_in = f_h - f.buplotpars.top*f_h; b_sz_in = f.subplotpars.bottom*f_h
    f_h = f_h*(1.+frac+pad)
    f.set_size_inches([f_w,f_h])
    # Now adjust subplot sizes
    return

def adjust_figure_single(f,cb_=False):
    adjust_figure(f,ax_w=0.8*fig_single_wid,fig_w=fig_single_wid)
    return

def adjust_figure_half(f):
    adjust_figure(f,ax_w=2.3,fig_w=fig_half_wid)
    return

def adjust_figure(f,ax_w,fig_w=None):
    """
    Adjust the margins on the given figure so that the width of the axis has the given size
    in inches, and the axis aspect ratio is the golden ratio.

    If the figure width isn't specified, the current figure width is used

    """
    if fig_w == None:
        fig_w = f.get_size_inches()[0]
        
    ax_h   = a_w * gr_i
    h_size = ax_w + left_margin_in
    v_size = ax_h + bottom_margin_in
    fig_h  = v_size + top_margin_in
    
    l_frac = left_margin_in / fig_w
    b_frac = bottom_margin_in / fig_h
    r_frac = h_size /fig_w
    t_frac = 1. - top_margin_in / fig_h

    f.set_size_inches(fig_w,fig_h)
    f.subplots_adjust(left=l_frac,right=r_frac,bottom=b_frac,top=t_frac)
    return
