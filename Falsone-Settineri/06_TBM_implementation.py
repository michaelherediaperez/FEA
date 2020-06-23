# -*- coding: utf-8 -*-
#
# Implementation of the TBM with a 2-noded FE.
# ------------------------------------------------------------------------------
# By       : Michael Heredia Pérez.
# Based on : Diego Andrés Álvarez Marin - daalvarez@unal.edu.co
# Date     : June/2020.
# e-mail   : mherediap@unal.edu.co
# Universidad Nacional de Colombia sede Manizales.
# ------------------------------------------------------------------------------ 


# Libraries.
import numpy as np
import matplotlib.pyplot as plt 
import pandas as pd 

# Constants so the code is easier to understand.
Y = 0; TH = 1
NL1 = 0; NL2 = 1


# Data can be red from a external file.
# Beam information, material propperties and loads (Modify this info if needed).
# ------------------------------------------------------------------------------
P  = 20              # kN      concentrated load in the middle x = L/2.
L  = 2               # m       lenght of the beam.
b  = 0.1             # m       transversal section base.   
h  = 0.3             # m       transversal section height.
E  = 2.1e8           # kPa     Young modulus.
nu = 0.3             # ad      Poisson coefficient.
G  = E/(2*(1+nu))    # kPa     shear modulus.
A  = b*h             # m**2    transversal section area.
I  = b*h**3/12       # m**4    moment of inercia.

alpha = 5/6          # ad      area shear coefficient.

# Supports relationship.
rest = np.array([[0, Y],     # Vertical displ restriction at DOF 0.
                 [8, Y]])    # Vertical displ restriction at DOF 8.

# FE organization.
# ------------------------------------------------------------------------------
# Nodes information.
xnod = np.linspace(0, L, 5)         # Nodes posittion.
nnds = xnod.shape[0]                # Number of nodes. 
Le = np.diff(xnod)[0]               # Lenght of the FE.
nef  = nnds -1                      # Number of FEs.
ndof = 2*nnds                       # Number od DOFs.

# Conectivity matriz from local to global LaG and DOFs matrix. 
LaG = np.c_[np.arange(0, nnds-1), np.arange(1, nnds)]
dof = np.c_[np.arange(0, ndof-1, 2), np.arange(1, ndof, 2)]

# Arrays with the information for each FE along the beam.

# Known (c) and Unknown (d) DOFs movements.
c  = rest[:, 0]     
ac = rest[:, 1]                     # Known magnitud of movements.
d  = d = np.setdiff1d(range(ndof), c)


# Assembling K and f.
# ------------------------------------------------------------------------------
idx = nef*[None]
K   = np.zeros((ndof, ndof))       # Initial Stiffness matrix.
f   = np.zeros(ndof)               # Equivalent nodal force vector.

f[4] = - P       # P applied in the middle.


for e in np.arange(nef):
    idx[e] = np.r_[ dof[LaG[e, NL1], :], dof[LaG[e, NL2], :]  ]
    
    # Stiffness bending matrix, 1 point of GL.
    Kb_T = (E*I/L) * np.array([[0,  0, 0,  0],
                               [0,  1, 0, -1],
                               [0,  0, 0,  0],
                               [0, -1, 0,  1] ])

    # Stiffness shear matrix, 1 point of GL.
    Ks1_T = (G*A*alpha/Le) * np.array([[   1,    Le/2,    -1,    Le/2],
                                       [Le/2, Le**2/4, -Le/2, Le**2/4],
                                       [  -1,   -Le/2,     1,   -Le/2],
                                       [Le/2, Le**2/4, -Le/2, Le**2/4] ])

    Ke = Kb_T + Ks1_T
    # FE (e) nodal equivalent forces vector for normal loads (here must go).

    K[np.ix_(idx[e], idx[e])] += Ke
    #f[idx[e],:]     += fe

# Show configuration of K
plt.figure()
plt.spy(K)
plt.title('Black points are elements different to zero')
plt.show()

# Extract the submatrix and specify the known quantities
# ------------------------------------------------------------------------------
# f = equivalent nodal forces vector.
# q = equilibrium nodal forces vector.
# a = movements.

#| qd |   | Kcc Kcd || ac |   | fd |  # qc=0 (siempre)
#|    | = |         ||    | - |    |
#| qc |   | Kdc Kdd || ad |   | fc |
Kcc = K[np.ix_(c,c)];  Kcd = K[np.ix_(c,d)]; fd = f[c]
Kdc = K[np.ix_(d,c)];  Kdd = K[np.ix_(d,d)]; fc = f[d]

# Solving the equation system.
ad = np.linalg.solve(Kdd, fc - Kdc@ac) # Unknown movements.
qd = Kcc@ac + Kcd@ad - fd              # Unknown equilibrium forces.

# Movements vector (a) and forces one (q)
a = np.zeros(ndof); q = np.zeros(ndof) # Save memory.
a[c] = ac;          a[d] = ad          
q[c] = qd         # q[d] = qc = 0      

# Movements inside the FE.
# ------------------------------------------------------------------------------
nint = 10                           # Number of points for interpolation.
xi   = np.linspace(-1, 1, nint)     # Normalized space in xi.

xx = nef * [None]       # Geometry interpolation.
vv = nef * [None]       # Displacement interpolation.
tt = nef * [None]       # Rotation angle interpolation.

for e in np.arange(nef):

    # Shape functions matrix.
    N = np.array([[xi**3/4 - 3*xi/4 + 1/2         ],
                  [xi**3/4 - xi**2/4 - xi/4 + 1/4 ],
                  [-xi**3/4 + 3*xi/4 + 1/2        ],
                  [xi**3/4 + xi**2/4 - xi/4 - 1/4 ] ]) 
 
    
    # Derivative of the shape functions matrix.
    dN_dxi = np.array([[3*xi**2/4 - 3/4        ],
                       [3*xi**2/4 - xi/2 - 1/4 ],
                       [3/4 - 3*xi**2/4        ],
                       [3*xi**2/4 + xi/2 - 1/4 ] ])

    # Nodal movement vector for the FE a^{(e)}
    ae = a[idx[e]]

    # Geometry interpolation. 
    xx[e] = Le*xi / 2 + ( xnod[LaG[e,NL1]] + xnod[LaG[e,NL2]] ) /2
        
    # Displacement inside the FE.
    vv[e] = N.T @ ae

    # Rotation angle inside the FE.
    tt[e] = np.arctan((dN_dxi.T*2/Le) @ ae)


# Preset result.
# ------------------------------------------------------------------------------

# Flatening the information.
xx = [val for sublist in xx for val in sublist]
vv = [val for sublist in vv for val in sublist]
tt = [val for sublist in tt for val in sublist]

# Info is set into a pandas dataframe.
table_xvt = pd.DataFrame(data = np.c_[xx, vv, tt],
                         columns = ['xx', 'vv', 'tt'])

# The results are saved into a MS EXEL spreadsheet.
filename = 'results_TBM.xlsx'
writer   = pd.ExcelWriter(filename, engine='xlsxwriter')
table_xvt.to_excel(writer, sheet_name='xvt') # pylint: disable=abstract-class-instantiated
writer.save()
print('\nInformation results are saved into a MS EXCEL spreadsheet.')

# END :)