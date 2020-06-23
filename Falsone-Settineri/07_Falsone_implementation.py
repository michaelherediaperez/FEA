# -*- coding: utf-8 -*-
#
# Implementation of the Falsone approximation for the TBM with a 2-noded FE.
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
  
    # FE (e) stiffness matrix calculated in 03_K_f_falsone.py 
    # Of course it is symetric

    # As we are considering E, I, A, G constant through the beam, K is constant.

    k1_1 =  3*A*E*G*I*Le*alpha/(4*(A*G*alpha + 3*E*I))
    k1_2 =  3*A*E*G*I*Le*alpha/(4*(A*G*alpha + 3*E*I))
    k1_3 = -3*A*E*G*I*Le*alpha/(4*A*G*alpha + 12*E*I)
    k1_4 =  3*A*E*G*I*Le*alpha/(4*(A*G*alpha + 3*E*I))

    k2_2 =  E*I*Le*(4*A*G*alpha + 3*E*I)/(4*(A*G*alpha + 3*E*I))
    k2_3 = -3*A*E*G*I*Le*alpha/(4*A*G*alpha + 12*E*I) 
    k2_4 =  E*I*Le*(2*A*G*alpha - 3*E*I)/(4*(A*G*alpha + 3*E*I))

    k3_3 =  3*A*E*G*I*Le*alpha/(4*(A*G*alpha + 3*E*I))
    k3_4 = -3*A*E*G*I*Le*alpha/(4*A*G*alpha + 12*E*I)

    k4_4 = E*I*Le*(4*A*G*alpha + 3*E*I)/(4*(A*G*alpha + 3*E*I))

    Ke = np.array([ [k1_1,  k1_2,  k1_3,  k1_4],
                    [k1_2,  k2_2,  k2_3,  k2_4],
                    [k1_3,  k2_3,  k3_3,  k3_4],
                    [k1_4,  k2_4,  k3_4,  k4_4] ])

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

# This vector (a) is actually the TBM with real variables.

# Movements inside the FE.
# ------------------------------------------------------------------------------
nint = 10                           # Number of points for interpolation.
xi   = np.linspace(-1, 1, nint)     # Normalized space in xi.

xx   = nef * [None]       # Geometry interpolation.
vv   = nef * [None]       # Displacement interpolation.
vv_b = nef * [None]       # Displacement interpolation.
tt   = nef * [None]       # Rotation angle interpolation.
MM   = nef * [None]       # Bending moment interpolation.
VV   = nef * [None]       # Shear force interpolation.

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

    # Second derivative of the shape functions matrix.
    dN2_dxi2 = np.array([[ 3*xi/2       ],
                         [ 3*xi/2 - 1/2 ],
                         [ -3*xi/2      ],
                         [ 3*xi/2 + 1/2 ]])

    # Third derivative of the shape functions matrix.
    dN3_dxi3 = np.array([[3/2 *xi/xi], 
                         [3/2 *xi/xi], 
                         [-3/2*xi/xi],
                         [3/2 *xi/xi]])

    # Correlation factor (I - C).inv()
    factor = np.array(
             [[(A*G*alpha + 3*E*I/2)/(A*G*alpha + 3*E*I), 
               -E*I*(4*A*G*alpha + 3*E*I)/(2*A*G*alpha*(A*G*alpha + 3*E*I)), 
               3*E*I/(2*(A*G*alpha + 3*E*I)), 
               -E*I*(2*A*G*alpha - 3*E*I)/(2*A*G*alpha*(A*G*alpha + 3*E*I))   ], 
              [0, 1, 0, 0], 
              [3*E*I/(2*(A*G*alpha + 3*E*I)), 
               E*I*(2*A*G*alpha - 3*E*I)/(2*A*G*alpha*(A*G*alpha + 3*E*I)), 
               (A*G*alpha + 3*E*I/2)/(A*G*alpha + 3*E*I), 
               E*I*(4*A*G*alpha + 3*E*I)/(2*A*G*alpha*(A*G*alpha + 3*E*I))    ], 
              [0, 0, 0, 1] ])

    # Fictitious nodal movement vector for the FE \bar{a}^{(e)}
    ae_b = factor @ a[idx[e]]

    # Geometry interpolation. 
    xx[e] = Le*xi / 2 + ( xnod[LaG[e,NL1]] + xnod[LaG[e,NL2]] ) /2
        
    # Fictitious displacement inside the FE.
    vv_b[e] = N.T @ ae_b
    
    # Correlations for the other important variables.
    tt[e]   = dN_dxi.T @ ae_b
    MM[e]   = E*I * dN2_dxi2.T  @ ae_b
    VV[e]   = E*I * dN3_dxi3.T  @ ae_b
    vv[e]   = vv_b[e] - E*I/(alpha*A*G) * dN2_dxi2.T  @ ae_b


# Preset result.
# ------------------------------------------------------------------------------

# Flatening the information.
xx = [val for sublist in xx for val in sublist]
vv = [val for sublist in vv for val in sublist]
tt = [val for sublist in tt for val in sublist]
MM = [val for sublist in MM for val in sublist]
VV = [val for sublist in VV for val in sublist]


# Info is set into a pandas dataframe.
table_xvt = pd.DataFrame(data = np.c_[xx, vv, tt, MM, VV],
                         columns = ['xx', 'vv', 'tt', 'MM', 'VV'])

# The results are saved into a MS EXEL spreadsheet.
filename = 'results_FS.xlsx'
writer   = pd.ExcelWriter(filename, engine='xlsxwriter')
table_xvt.to_excel(writer, sheet_name='xvt') # pylint: disable=abstract-class-instantiated
writer.save()
print('\nInformation results are saved into a MS EXCEL spreadsheet.')

# Graphs.
# ------------------------------------------------------------------------------

plt.figure()
plt.plot(xx, vv, 'b-')
plt.grid()
plt.title('Vertical displacement')
plt.xlabel(r'$x $ [m]')
plt.ylabel(r'$v(x)$ [m]')
plt.grid('on')
plt.show()

plt.figure()
plt.plot(xx, tt, 'b-')
plt.grid()
plt.title('Transversal section angle rotation')
plt.xlabel(r'$x $ [m]')
plt.ylabel(r'$\theta(x)$ [rad]')
plt.grid('on')
plt.show()

plt.figure()
plt.plot(xx, MM, 'b-')
plt.grid()
plt.title('Bending moment (mathematical graph)')
plt.xlabel(r'$x $ [m]')
plt.ylabel(r'$M(x)$ [kN.m]')
plt.grid('on')
plt.show()

plt.figure()
plt.plot(xx, VV, 'b-')
plt.grid()
plt.title('Shear Force')
plt.xlabel(r'$x $ [m]')
plt.ylabel(r'$V(x)$ [kN]')
plt.grid('on')
plt.show()

# VERTICAL DISPLACEMENTS AND ANGLE ROTATIONS ARE TO BIG. 
# STILL LOOKING FOR THE ERROR.

# END :)