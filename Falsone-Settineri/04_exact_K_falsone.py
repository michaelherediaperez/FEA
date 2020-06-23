# -*- coding: utf-8 -*-
#
# Exact Stiffness Matrix for n = 4 following the Falsone article.
# ------------------------------------------------------------------------------
# By       : Michael Heredia Pérez.
# Date     : June/2020.
# e-mail   : mherediap@unal.edu.co
# Universidad Nacional de Colombia sede Manizales.
# ------------------------------------------------------------------------------ 

# Libreries
import sympy as sp

# Symbolic variables.
x, L, V, M, t, v, vb = sp.symbols('x L V M t v vb')
alpha, A, G, E, I = sp.symbols('alpha A G E I')
C1, C2, C3, C4 = sp.symbols('C1, C2, C3, C4')  
q = 0

# Saving memory for the stiffness matrix.
K = sp.zeros(4)

# As border conditions we have the following differential equations.


V = sp.integrate(q, x) + C1  
M = sp.integrate(V, x) + C2
t = (sp.integrate(M, x) + C3)/EI
w = sp.integrate(t, x) + C4

for i in range(4):

    sol = sp.solve([w.subs(x, 0) - int((i == 0)),  # Condiciones de frontera
                    t.subs(x, 0) - int((i == 1)),
                    w.subs(x, L) - int((i == 2)),
                    t.subs(x, L) - int((i == 3))],
                   [C1, C2, C3, C4])

    constantes = [(C1, sol[C1]), (C2, sol[C2]), (C3, sol[C3]), (C4, sol[C4])]

    K[:, i] = [+ (V.subs(constantes)).subs(x, 0),   # Y1  se evaluan las
               - (M.subs(constantes)).subs(x, 0),   # M1  reacciones verticales
               - (V.subs(constantes)).subs(x, L),   # Y2  y los momentos en los
               + (M.subs(constantes)).subs(x, L)]  # M2  apoyos

# %%Se imprime la solución
print(f'K_EB = (EI/L**3) * \n{sp.pretty(K/(EI/L**3))}')

# Fin del programa