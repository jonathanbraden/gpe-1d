#!/usr/bin/env python
import palettable as pal

# Formatting for Nonlinear BEC Paper
cosCMap = pal.colorbrewer.diverging.RdYlBu_11_r.mpl_colormap #'coolwarm'
rhoCMap = 'PuOr'

# Notation for Nonlinear BEC Paper
lPar = r'\lambda'
densTot = r'\varrho'
densRel = r'\varepsilon'
phaseTot = r'\vartheta'
phaseRel = r'\varphi'

gS = r'g_{\rm s}'
rhoBG = r'\bar{n}'
nuZero = r'\nu_0'
nuAmp = r'\delta'
nuFreq = r'\omega'

tDim = r'\hbar^{-1}%s%s' % (gS,rhoBG)
xDim = r'\hbar^{-1}\sqrt{%s%s m}' % (gS,rhoBG)
kDim = r'\frac{\hbar k_{\rm nyq}}{\sqrt{%s%s m}}' % (gS,rhoBG)

xUnit = r'$\hbar^{-1}\sqrt{g_{\rm s}\bar{n}m}\, x$'
tUnit = r'$\hbar^{-1}g_{\rm s}\bar{n}\,t$'
xBar = r'$\bar{x}$'
tBar = r'$\bar{t}$'

xLab = xUnit
tLab = tUnit
