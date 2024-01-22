# Discretization of governing equations
This chapter provide concrete discretization of heat and mass transfer equations presented in technical document. Details on how discretization works are given in Chapter 2 and here most of details will be omitted.

## Notation
Some of common notation will be used in this document:

{V} - vector of variables or constants
[M] - matrix
$\psi \equiv \psi(x,y)$ - shape function in global coordinate system

$\tilde{\psi} \equiv \psi(\xi, \eta)$​ - shape function in local coordinate system​

$V^{(a)}$ - variable from iteration "a"

$V_{[i]}$ - variable from timestep "i"

$[M]_{ij}$ - matrix element at position (i, j)

$\{V\}_i$ - vector element at position "i"

## Pressure

Pressure equation from technical document is given in a following form:

$\frac{\partial \gamma_{air}}{\partial P_{air}} \frac{\partial P_{air}}{\partial t} = \frac{\partial}{\partial x}(\frac{k_{air} \cdot \rho_{air}}{\mu_{air}} \cdot \frac{\partial P_{air}}{\partial x}) + \frac{\partial}{\partial y}(\frac{k_{air} \cdot \rho_{air}}{\mu_{air}} \cdot \frac{\partial P_{air}}{\partial y})$

When discretizing this equation, it will be considered that $\frac{k_{air} \cdot \rho_{air}}{\mu_{air}}$ and $\frac{\partial \gamma_{air}}{\partial P_{air}}$ are constant. Because that is not true, it is important to note that above equation will be considered non-linear where solution will be found through iterations. Meaning that program will solve for the  set of $\Delta P_{air}$ and then will recalculate $\frac{k_{air} \cdot \rho_{air}}{\mu_{air}}$ and $\frac{\partial \gamma_{air}}{\partial P_{air}}$ for new solution.

Iterations will be repeated till tolerance is achieved.































