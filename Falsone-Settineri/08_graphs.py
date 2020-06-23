# -*- coding: utf-8 -*-
#
# Graphic comparison between TBM and EBBM and the Falsone-Settineri approach.
# ------------------------------------------------------------------------------
# By       : Michael Heredia Pérez.
# Date     : June/2020.
# e-mail   : mherediap@unal.edu.co
# Universidad Nacional de Colombia sede Manizales.
# ------------------------------------------------------------------------------ 

# Libraries.
import numpy as np
import matplotlib.pyplot as plt 
import pandas as pd 

# Reading EBBM file.
EBBM = pd.read_excel('results_EBBM.xlsx', sheet_name=None)

xx_EBBM = EBBM['xvt']['xx']
vv_EBBM = EBBM['xvt']['vv']
tt_EBBM = EBBM['xvt']['tt']

# Reading EBBM file.
TBM = pd.read_excel('results_TBM.xlsx', sheet_name=None)

xx_TBM = TBM['xvt']['xx']
vv_TBM = TBM['xvt']['vv']
tt_TBM = TBM['xvt']['tt']

# Reading Falsone-Settineri file
FS = pd.read_excel('results_FS.xlsx', sheet_name=None)

xx_FS = FS['xvt']['xx']
vv_FS = FS['xvt']['vv']
tt_FS = FS['xvt']['tt']
MM_FS = FS['xvt']['MM']
VV_FS = FS['xvt']['VV']

# Graphs.
# ------------------------------------------------------------------------------


plt.figure()
plt.plot(xx_EBBM, vv_EBBM, 'b-', label = 'EBBM')
plt.plot( xx_TBM,  vv_TBM, 'g-', label = 'TBM')
plt.plot(  xx_FS,   vv_FS, 'm-', label = 'FS')
plt.grid()
plt.legend()
plt.title('Vertical displacement')
plt.xlabel(r'$x $ [m]')
plt.ylabel(r'$v(x)$ [m]')
plt.grid('on')
plt.show()

plt.figure()
plt.plot(xx_EBBM, tt_EBBM, 'b-', label = 'EBBM')
plt.plot( xx_TBM,  tt_TBM, 'g-', label = 'TBM')
plt.plot(  xx_FS,   tt_FS, 'm-', label = 'FS')
plt.grid()
plt.legend()
plt.title('Transversal section angle rotation')
plt.xlabel(r'$x $ [m]')
plt.ylabel(r'$\theta(x)$ [rad]')
plt.grid('on')
plt.show()