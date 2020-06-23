# -*- coding: utf-8 -*-
#
# Hermite polynimials - 2-noded EBBM FE shape functions. 
# ------------------------------------------------------------------------------
# By       : Michael Heredia Pérez.
# Based on : Diego Andrés Álvarez - daalvare@unal.edu.co
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
xi, v1, v2, v1p, v2p = sp.symbols('xi v1 v2 v1p v2p')
xi1 = -1 ; xi2 = 1

# Lagrangians plynomials.
L1 = sp.interpolate([(-1, 1), (1, 0)], xi)
L2 = sp.interpolate([(-1, 0), (1, 1)], xi)

# Pi and Qi polynomials.
P1 = (1 - 2*(sp.diff(L1, xi)).subs(xi, xi1)*(xi-xi1))*L1**2
P2 = (1 - 2*(sp.diff(L2, xi)).subs(xi, xi2)*(xi-xi2))*L2**2
Q1 = (xi-xi1)*L1**2
Q2 = (xi-xi2)*L2**2

# Polynomial interplation.
v =  sp.simplify(P1*v1 + P2*v2 + Q1*v1p + Q2*v2p)

# Calculate de shape functions.
N1  = sp.expand(v.subs([(v1, 1), (v1p, 0), (v2, 0), (v2p, 0)]))
N1b = sp.expand(v.subs([(v1, 0), (v1p, 1), (v2, 0), (v2p, 0)]))
N2  = sp.expand(v.subs([(v1, 0), (v1p, 0), (v2, 1), (v2p, 0)]))
N2b = sp.expand(v.subs([(v1, 0), (v1p, 0), (v2, 0), (v2p, 1)]))

# Pretty results.
print('-'*60)
print('N_1(x)  ='); sp.pprint(N1)
print('Nb_1(x) ='); sp.pprint(N1b) 
print('N_2(x)  ='); sp.pprint(N2)
print('Nb_2(x) ='); sp.pprint(N2b) 

# LaTex format results.
print('-'*60)
print(r'$N_1(x) = $');       print(sp.latex(N1))
print(r'$\bar{N_1}(x) = $'); print(sp.latex(N1b)) 
print(r'$N_2(x) = $');       print(sp.latex(N2))
print(r'$\bar{N_2}(x) = $'); print(sp.latex(N2b)) 

# Python results.
print('-'*60)
print('N1(x)  = '); print(N1)
print('N1b(x) = '); print(N1b) 
print('N2(x)  = '); print(N2)
print('N2b(x) = '); print(N2b) 

# Making evaluable the shape functions.
f_N1  = sp.lambdify([xi], N1)
f_N1b = sp.lambdify([xi], N1b)
f_N2  = sp.lambdify([xi], N2)
f_N2b = sp.lambdify([xi], N2b)
xxi   = np.linspace(xi1, xi2, num = 100)

# Grapghic or the 2-noded EBBM FE shape functions.
plt.figure()
plt.plot(xxi, f_N1(xxi),  label = r'$N_1(\xi)$')
plt.plot(xxi, f_N1b(xxi), label = r'$\bar{N_1}(\xi)$')
plt.plot(xxi, f_N2(xxi),  label = r'$N_2(\xi)$')
plt.plot(xxi, f_N2b(xxi), label = r'$\bar{N_2}(\xi)$')
plt.xlabel(r'$\xi$')
plt.title('2-noded EBBM FE shape functions')
plt.grid()
plt.axis('equal')
plt.legend()
plt.show()

# End :)