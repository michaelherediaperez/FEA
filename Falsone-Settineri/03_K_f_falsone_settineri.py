# -*- coding: utf-8 -*-
#
# Stiffness Matrix and equivalent nodal forces vector for n = 4.
# ------------------------------------------------------------------------------
# By       : Michael Heredia Pérez.
# Date     : June/2020.
# e-mail   : mherediap@unal.edu.co
# Universidad Nacional de Colombia sede Manizales.
# ------------------------------------------------------------------------------ 

# Libraries.
import sympy as sp 
sp.init_printing(pretty_print = True)

# Symbolic variables.
xi, E, I, alpha, A, G, x1, x2, Le, q= sp.symbols('xi E I alpha A G x1 x2 Le q')
xi1 = -1 ; xi2 = 1

# Hermite shape functions.
N1  = xi**3/4 - 3*xi/4 + 1/2
N1b = xi**3/4 - xi**2/4 - xi/4 + 1/4
N2  = -xi**3/4 + 3*xi/4 + 1/2
N2b = xi**3/4 + xi**2/4 - xi/4 - 1/4

# The jacobian is:
J = Le/2

# Hermite Shape function vector and fictious nodal displacement vector.
N  = sp.Matrix([[N1, N1b, N2, N2b]])

# The first and second derivatives of N.
dN_dxi   = sp.diff(N, xi)
dN2_dxi2 = sp.diff(dN_dxi, xi)
dN3_dxi3 = sp.diff(dN2_dxi2, xi)

print('dN_dxi = ');   print(dN_dxi)
print('dN2_dxi2 = '); print(dN2_dxi2)
print('dN3_dxi3 = '); print(dN3_dxi3)

# Empty vector. 
empty = sp.zeros(1, 4)

# Tha matrix C for the FE (e) with n = 4.
C = (E*I)/(alpha*A*G) * sp.Matrix([dN2_dxi2.subs(xi, xi1), 
                                   empty, 
                                   dN2_dxi2.subs(xi, xi2), 
                                   empty])

# Identity matrix of order n = 4.
In = sp.eye(4)

# Factor of inv(I - C)
factor = (In - C).inv(method = "LU")

# The stiffness matrix for the fictitious variable.
Kb = sp.simplify(sp.integrate(dN2_dxi2.T*E*I*dN2_dxi2*J, (xi, -1, 1)))
print('Kb = '); sp.pprint(Kb, num_columns=150); print()

# The TBM stiffness matrix after the adjustment.
K = sp.simplify(Kb * factor)
print('K = '); sp.pprint(K, num_columns=150); print()

# The equivalent nodal forces vector
f = sp.simplify(sp.integrate(N.T*q*J, (xi, -1, 1)))
print('f = '); sp.pprint(f, num_columns=150); print()

print(sp.simplify(factor))