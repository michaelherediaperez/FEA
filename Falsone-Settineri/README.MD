# AN EULER-BERNOULLI-LIKE FINITE ELEMENT METHOD FOT TIMOSHENKO BEAMS
## Implementation of the G. Falsone and D. Settineri article.

Article can be found here:


Based on the formulation made by Falsone and Settineri, the procedure to simplify the problem of TBM to a EBBM-like FE is presentened in the following document:


## Stiffness matrix and Equivalent nodal forces vector

The computation of **K** and **f** for the 2-noded Timoshenko finite element proposed by the authors was made employing numerical integration over a Gauss-Legendre cuadrature domain. The results are shown here:

```
K =

⎡   3⋅A⋅E⋅G⋅I⋅Le⋅α        3⋅A⋅E⋅G⋅I⋅Le⋅α            -3⋅A⋅E⋅G⋅I⋅Le⋅α           3⋅A⋅E⋅G⋅I⋅Le⋅α       ⎤    
⎢ ─────────────────   ─────────────────      ────────────────      ─────────────────   ⎥ 
⎢ 4⋅(A⋅G⋅α + 3⋅E⋅I)     4⋅(A⋅G⋅α + 3⋅E⋅I)          4⋅A⋅G⋅α + 12⋅E⋅I         4⋅(A⋅G⋅α + 3⋅E⋅I)    ⎥
⎢                                                                                      ⎥
⎢   3⋅A⋅E⋅G⋅I⋅Le⋅α     E⋅I⋅Le⋅(4⋅A⋅G⋅α + 3⋅E⋅I)      -3⋅A⋅E⋅G⋅I⋅Le⋅α       E⋅I⋅Le⋅(2⋅A⋅G⋅α - 3⋅E⋅I)  ⎥
⎢ ─────────────────  ────────────────────   ────────────────     ────────────────────  ⎥
⎢ 4⋅(A⋅G⋅α + 3⋅E⋅I)     4⋅(A⋅G⋅α + 3⋅E⋅I)          4⋅A⋅G⋅α + 12⋅E⋅I      4⋅(A⋅G⋅α + 3⋅E⋅I)       ⎥
⎢                                                                                      ⎥
⎢ -3⋅A⋅E⋅G⋅I⋅Le⋅α        -3⋅A⋅E⋅G⋅I⋅Le⋅α            3⋅A⋅E⋅G⋅I⋅Le⋅α       -3⋅A⋅E⋅G⋅I⋅Le⋅α           ⎥
⎢ ────────────────   ────────────────       ─────────────────      ────────────────    ⎥
⎢ 4⋅A⋅G⋅α + 12⋅E⋅I       4⋅A⋅G⋅α + 12⋅E⋅I          4⋅(A⋅G⋅α + 3⋅E⋅I)      4⋅A⋅G⋅α + 12⋅E⋅I       ⎥
⎢                                                                                      
⎢   3⋅A⋅E⋅G⋅I⋅Le⋅α      E⋅I⋅Le⋅(2⋅A⋅G⋅α - 3⋅E⋅I)    -3⋅A⋅E⋅G⋅I⋅Le⋅α       E⋅I⋅Le⋅(4⋅A⋅G⋅α + 3⋅E⋅I)  ⎥
⎢ ─────────────────  ─────────────────────  ────────────────   ────────────────────   ⎥
⎣ 4⋅(A⋅G⋅α + 3⋅E⋅I)     4⋅(A⋅G⋅α + 3⋅E⋅I)         4⋅A⋅G⋅α + 12⋅E⋅I      4⋅(A⋅G⋅α + 3⋅E⋅I)       ⎦

```

```
f =

⎡       0.5⋅Le⋅q         ⎤
⎢                       ⎥
⎢0.166666666666667⋅Le⋅q  ⎥
⎢                       ⎥
⎢       0.5⋅Le⋅q         ⎥
⎢                       ⎥
⎣-0.166666666666667⋅Le⋅q ⎦

```

Both were calculated with this code:


## Implementation example

The beam shown in the figure was analized with the folllowing code:

<img src="figs/beam.png">

The properties of the beam can be found in the code. It is worth to say that it may be an error in the computation of the matrix **K**, because as we see in the graphics of the vertical displacement, rotation angle, bending moment and shear force, ALL OF THEM ARE MORE THAN WHAT THEY ARE SUPPOSED TO BE.

<img src="figs/vv.png">
<img src="figs/tt.png">
<img src="figs/MM.png">
<img src="figs/VV_s.png">


### NOTE

I am still looking forward for the error.