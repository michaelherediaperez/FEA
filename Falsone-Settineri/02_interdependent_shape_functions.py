# -*- coding: utf-8 -*-
#
# Matris C and Interdependent shape functions for nodal movements 
# ------------------------------------------------------------------------------
# By       : Michael Heredia Pérez.
# Date     : June/2020.
# e-mail   : mherediap@unal.edu.co
# Universidad Nacional de Colombia sede Manizales.
# ------------------------------------------------------------------------------ 

# Libreries.
import sympy as sp 
import numpy as np
import matplotlib.pyplot as plt
sp.init_printing(pretty_print = True)

# Symbolic variables.
xi, E, I, alpha, A, G, L = sp.symbols('xi E I alpha A G L')
xi1 = -1 ; xi2 = 1

# Hermite shape functions.
N1  = xi**3/4 - 3*xi/4 + 1/2
N1b = xi**3/4 - xi**2/4 - xi/4 + 1/4
N2  = -xi**3/4 + 3*xi/4 + 1/2
N2b = xi**3/4 + xi**2/4 - xi/4 - 1/4

# Hermite Shape function vector.
N = sp.Matrix([[N1, N1b*L/2, N2, N2b*L/2]])

# The first and second derivatives of N.
dN_dxi = sp.diff(N, xi)
dN2_dxi2 = sp.diff(dN_dxi, xi)

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

# Interdependent shape functions for vertical displacement v(x) 
Nv = (N - (E*I)/(alpha*A*G)*dN2_dxi2) @ factor

# Interdependent shape functions for rotation angle theta(x) 
Nt = -dN_dxi @ factor

# Python results because pretties are awfuls .
print('*'*60)
print('C =  '); sp.pprint(C)       # Just this one is pretty.
print('C =  '); print(C)
print('Nt = '); print(Nt)
print('Nv = '); print(Nv)

# Making evaluable the shape functions.
f_Nv1 = sp.lambdify([xi, E, I, A, G, alpha, L], Nv[0])
f_Nv2 = sp.lambdify([xi, E, I, A, G, alpha, L], Nv[1])
f_Nv3 = sp.lambdify([xi, E, I, A, G, alpha, L], Nv[2])
f_Nv4 = sp.lambdify([xi, E, I, A, G, alpha, L], Nv[3])

f_Nt1 = sp.lambdify([xi, E, I, A, G, alpha, L], Nt[0])
f_Nt2 = sp.lambdify([xi, E, I, A, G, alpha, L], Nt[1])
f_Nt3 = sp.lambdify([xi, E, I, A, G, alpha, L], Nt[2])
f_Nt4 = sp.lambdify([xi, E, I, A, G, alpha, L], Nt[3])

xxi  = np.linspace(xi1, xi2, num = 100)

# Giving values to the constants for the graphs.
b = 0.1; h = 0.25
E = 200e9
A = b*h
G = 70e9
I = b*h**3/12
alpha = 5/6
L = 10

# Grapghic of the 2-noded TBM FE interdependent displacement shape functions.
plt.figure()
plt.plot(xxi, f_Nv1(xxi, E, I, A, G, alpha, L), label = r'$N_{v1}(\xi)$')
plt.plot(xxi, f_Nv2(xxi, E, I, A, G, alpha, L), label = r'$N_{v2}(\xi)$')
plt.plot(xxi, f_Nv3(xxi, E, I, A, G, alpha, L), label = r'$N_{v3}(\xi)$')
plt.plot(xxi, f_Nv4(xxi, E, I, A, G, alpha, L), label = r'$N_{v4}(\xi)$')
plt.xlabel(r'$\xi$')
plt.title('2-noded TBM FE interdependent displacemnt shape functions')
plt.grid()
plt.axis('equal')
plt.legend()
plt.show()

# Grapghic of the 2-noded TBM FE interdependent displacement shape functions.
plt.figure()
plt.plot(xxi, f_Nt1(xxi, E, I, A, G, alpha, L), label = r'$N_{\theta 1}(\xi)$')
plt.plot(xxi, f_Nt2(xxi, E, I, A, G, alpha, L), label = r'$N_{\theta 2}(\xi)$')
plt.plot(xxi, f_Nt3(xxi, E, I, A, G, alpha, L), label = r'$N_{\theta 3}(\xi)$')
plt.plot(xxi, f_Nt4(xxi, E, I, A, G, alpha, L), label = r'$N_{\theta 4}(\xi)$')
plt.xlabel(r'$\xi$')
plt.title('2-noded TBM FE interdependent angle rotation shape functions')
plt.grid()
plt.axis('equal')
plt.legend()
plt.show()

# We most choose the adecuated interdependent shape functions.
# End :)