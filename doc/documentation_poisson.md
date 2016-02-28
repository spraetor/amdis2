# AMDiS - Adaptive MultiDimensional Simulations

## Poisson equation

As example for a stationary problem, we choose the Poisson equation

$$-\Delta u = f ~~~\mbox{in } \Omega \subset \mathbb{R}^{dim}$$

$$u = g ~~~ \mbox{on } \partial\Omega$$

with

$$f(x) = -\left( 400x^2 - 20 dow \right) e^{-10x^2}$$

$$g(x) = e^{-10x^2}.$$

$$dim$$ is the problem dimension and thus the dimension of the mesh the problem will be discretized on. $$dow$$ is the dimension of the world the problem lives in. So world coordinates are always real valued vectors of size $$dow$$. Note that the problem is defined in a dimension independent way. Furthermore, $$dow$$ has not to be equal to $$dim$$ as long as $$1 \leq dim \leq dow \leq 3$$ holds.

Although the implementation described below is dimension independent, we focus on the case $$dim = dow = 2$$ for the rest of this section. The analytical solution on $$\Omega = [0,1]\times [0,1]$$ is shown in the following Figure:

### Source Code

For this first example, we give the complete source code. But to avoid
loosing the overview, we sometimes interrupt the code to give some
explaining comment. The first three lines of the application code are:

``` c++
#include "AMDiS.h"
#include "GenericOperatorTerm.h"
using namespace AMDiS;
```

`AMDiS` to avoid potential naming conflicts with other libraries.

Now, the functions $$f$$ and $$g$$ will be defined by the functor-classes $$F$$
and $$G$$:

``` c++
class G : public FunctorBase
{
public:
  // define the return type of the functor
  typedef double value_type;

  // degree of a polynomial that approximates the function best
  int getDegree(int d0) const
  {
    return 4*d0; 
  }

  // the actual evaluation of the functor
  double operator()(const WorldVector<double>& x) const
  {
    return exp(-10.0 * (x * x));
  }
};
```

`G` is a sub class of the base class `FunctorBase` that represents 
arbitrary functors. Here, we want to define a mapping from
$$\mathbb{R}^{dow}$$, implemented by the class
`WorldVector<double>`, to $$\mathbb{R}$$, represented by the data
type `double`. The actual mapping is defined by overloading the
`operator()`. `x*x` stands for the scalar product of vector
`x` with itself.

For a functor it is assumed that it provies a `value_type` and a function `getDegree(...)`
that returns the polynomial degree used to choose an appropriate quadrature rule during assembling. The arguments
of the function `getDegree(...)` are the polynomial degrees of the objects passed to the `operator()`. 
For coordinates it is always 1, but a functor can takes values of a `DOFVectors` where now degree is related to the 
polynomial degree of the underlying basis functions.

The class`F` is defined in a similar way:

``` c++
class F : public FunctorBase
{
public:
  typedef double value_type;
  int getDegree(int d0) const { return 4*d0; }

  double operator()(const WorldVector<double>& x) const 
  {
    int dow = Global::getGeo(WORLD);
    double r2 = (x * x);
    double ux = exp(-10.0 * r2);
    return -(400.0 * r2 - 20.0 * dow) * ux;
  }
};
```

`F` gets the world dimension from the class `Global` by a
call of the static function `getGeo(WORLD)`. 

Now, we start with the main program:

``` c++
// ===== main program =====
int main(int argc, char* argv[])
{
  // ===== initialize AMDiS =====
  AMDiS::init(argc,argv);
```

The first command in each AMDiS-program should be `AMDiS::init(...)`. 
For arguments passed to the executable it is assumed that the first parameter is
the filename of a parameter-file (described below). The init-method reads this parameter-file
and in parallel initializes also MPI communicators and other necessary data-structures.

Now, a stationary problem with name `ellipt` is created and initialized:  

``` c++
  // ===== create and init the scalar problem ===== 
  ProblemStat ellipt("ellipt");
  ellipt.initialize(INIT_ALL);
```

The name argument of the problem is used to identify parameters in the
parameter file that belong to this problem. In this case, all
parameters with prefix `ellipt->` are associated to this
problem. The initialization argument `INIT_ALL` means that all
problem modules are created in a standard way. Those are: The finite
element space including the corresponding mesh, the required system
matrices and vectors, an iterative solver, an estimator, a marker, and
a file writer for the computed solution. The initialization of these
components can be controlled through the parameter-file.

The next steps are the creation of the adaptation loop and the
corresponding `AdaptInfo`:

``` c++
  // === create adapt info ===
  AdaptInfo adaptInfo("ellipt->adapt");

  // === create adapt ===
  AdaptStationary adapt("ellipt->adapt", ellipt, adaptInfo);
```

The `AdaptInfo` object contains information about the current
state of the adaptation loop as well as user given parameters
concerning the adaptation loop, like desired tolerances or maximal
iteration numbers. Using `adaptInfo`, the adaptation loop can be
inspected and controlled at runtime.  Now, a stationary adaptation
loop is created, which implements the standard ''assemble-solve-estimate-adapt'' loop. Arguments are the name, again
used as parameter prefix, the problem as implementation of an
iteration interface, and the `AdaptInfo` object. The adaptation
loop only knows when to perform which part of an iteration. The
implementation and execution of the single steps is delegated to an
iteration interface, here implemented by the stationary problem
`ellipt`.

The operators now are defined as follows:

``` c++
  // ===== create matrix operator =====
  Operator matrixOperator(ellipt.getFeSpace(), ellipt.getFeSpace());
  addSOT(matrixOperator, 1.0);
  ellipt.addMatrixOperator(matrixOperator, 0, 0);

  // ===== create rhs operator =====
  Operator rhsOperator(ellipt.getFeSpace());
  addZOT(rhsOperator, function_<F>(X()) );
  ellipt.addVectorOperator(rhsOperator, 0);
```

We define a matrix operator (left hand side operator) on the finite
element space of the problem and a vector operator (right hand side operator).
Therefore we write the differential equation above in a weak form: Find $$u\in V_g$$, s.t.

$$\langle \nabla u, \nabla v\rangle = \langle f, v\rangle\qquad\forall v\in V_0$$

Here $$u$$ denotes the trial function and $$v$$ the test function. The 
space $$V$$ denotes the finite-element space, where the subscript indices indicate the
Dirichlet-boundary conditions.

All operators in AMDiS are constructed by defining the coefficient-functions in front of the second-order terms
(terms containing two derivatives), first-order terms (just one derivative) and zero-order terms. A general
second-order term would look like

$$\langle A(x,u,\nabla u)\cdot \nabla u, \nabla v\rangle$$

with a matrix-valued function $$A$$. The special case of a scaled identity matrix can simply be written as 
a scalar-valued coefficient function. In the case of the Poisson-equation above, the coefficient-function is simply given by
$$A \equiv 1$$ and the coefficient-function on the right-hand side by $$f(x)$$

While $$A$$ can be implemented as `addSOT(operator, 1.0)`, or `addSOT(operator, constant(1.0))` 
the right-hand side term $$f$$ must be evaluated at the coordinates of the quadrature points. This can be achieved 
using the expression

``` c++
addZOT(operator, function_(F(), X()) );
// or
addZOT(operator, function_<F>(X()) );
```

Now, we define boundary conditions:

``` c++
  // ===== add boundary conditions =====
  ellipt.addDirichletBC(1, 0, 0, new G);
```

We have one Dirichlet boundary condition associated with identifier
$$1$$. All nodes belonging to this boundary are set to the value of
function `G` at the corresponding coordinates. In the macro file the Dirichlet boundary is marked
with identifier $$1$$, too. So the nodes can be uniquely determined. As
with adding operators to the operator matrix, we have to define the
operator, on which the boundary condition will be applied. Thus we
have to provide the matrix position indices after the boundary
identifier.

Finally we start the adaptation loop and afterwards write out the results:

``` c++
  // ===== start adaption loop =====
  adapt.adapt();

  // ===== write result =====
  ellipt.writeFiles(adaptInfo, true);
}
```

The second argument of `writeFiles` forces the file writer to
print out the results. In time dependent problems it can be useful to
write the results only every $$i$$-th timestep. To allow this behavior
the second argument has to be `false`.

== Parameter file ==
The name of the parameter file must be given as command line
argument. In the 2d case we call:

```
> ./ellipt init/ellipt.dat.2d
```

In the following, the content of file `init/ellipt.dat.2d` is
described:

```
dimension of world:                2

elliptMesh->macro file name:       ./macro/macro.stand.2d
elliptMesh->global refinements:    0
```

The dimension of the world is 2, the macro mesh with name
`elliptMesh` is defined in file `./macro/macro.stand.2d`. The mesh is not globally refined
before the adaptation loop. A value of $$n$$ for
`elliptMesh->global refinements` means $$n$$ bisections of every
macro element. Global refinements before the adaptation loop can be
useful to save computation time by starting adaptive computations with
a finer mesh.

```
ellipt->mesh:                      elliptMesh
ellipt->dim:                       2
ellipt->polynomial degree[0]:      1
ellipt->components:                1
```

Now, we construct the finite element space for the problem
`ellipt` (see Section ''''Source Code''''). We use the mesh
`elliptMesh`, set the problem dimension to 2, and choose Lagrange
basis functions of degree 1. The number of components, i.e.
variables, in the equation is set to 1, since we are about to define a
scalar PDE.

```
ellipt->solver:                    cg   
ellipt->solver->max iteration:     1000
ellipt->solver->tolerance:         1.e-8
ellipt->solver->left precon:       diag
ellipt->solver->right precon:      no
```

We use the ''conjugate gradient method'' as iterative solver. The
solving process stops after maximal $$1000$$ iterations or when a
tolerance of $$10^{-8}$$ is reached. Furthermore, we apply diagonal left
preconditioning, and no right preconditioning.

```
ellipt->estimator[0]:                 residual
ellipt->estimator[0]->error norm:     1
ellipt->estimator[0]->C0:             0.1
ellipt->estimator[0]->C1:             0.1
```

As error estimator we use the residual method. The used error norm is
the H1-norm (instead of the L2-norm: 2). Element residuals (C0) and
jump residuals (C1) both are weighted by factor $$0.1$$.

```
ellipt->marker[0]->strategy:          2     % 0: no 1: GR 2: MS 3: ES 4:GERS
ellipt->marker[0]->MSGamma:           0.5
```

After error estimation, elements are marked for refinement and
coarsening. Here, we use the maximum strategy with $$\gamma = 0.5$$.

```
ellipt->adapt[0]->tolerance:          1e-4
ellipt->adapt[0]->refine bisections:  2

ellipt->adapt->max iteration:      10
```

The adaptation loop stops, when an error tolerance of $$10^{-4}$$ is
reached, or after maximal $$10$$ iterations. An element that is marked
for refinement, is bisected twice within one iteration. Analog
elements that are marked for coarsening are coarsened twice per
iteration.

```
ellipt->output->filename:          output/ellipt
ellipt->output->ParaView format:   1
```

The result is written in ParaView-format to the file
`output/ellipt.vtu`.

\subsection{Macro file}
[[Bild::ellipt_macro.png|thumb|left|400px|Two dimensional macro mesh.]]
In the Figure one can see the macro mesh which is
described by the file `macro/macro.stand.2d`.  First, the
dimension of the mesh and of the world are defined:

```
DIM: 2
DIM_OF_WORLD: 2
```

Then the total number of elements and vertices are given:

```
number of elements: 4  
number of vertices: 5  
```

The next block describes the two dimensional coordinates of the five vertices:  

```
vertex coordinates:
 0.0 0.0  
 1.0 0.0
 1.0 1.0
 0.0 1.0
 0.5 0.5
```

The first two numbers are interpreted as the coordinates of vertex 0,
and so on.

Corresponding to these vertex indices now the four triangles are
given:

```
element vertices:
0 1 4 
1 2 4  
2 3 4 
3 0 4 
```

Element 0 consists in the vertices 0, 1 and 4. The numbering is done
anticlockwise starting with the vertices of the longest edge.

It follows the definition of boundary conditions:

```
element boundaries:
0 0 1 
0 0 1
0 0 1 
0 0 1 
```

The first number line means that element 0 has no boundaries at edge 0
and 1, and a boundary with identifier 1 at edge 2. The edge with
number $$i$$ is the edge opposite to vertex number $$i$$. The boundary
identifier 1 corresponds to the identifier 1 we defined within the
source code for the Dirichlet boundary. Since all elements of the
macro mesh have a Dirichlet boundary at edge 2, the line `0 0 1`
is repeated three times.

The next block defines element neighborships. `-1` means there is
no neighbor at the corresponding edge. A non-negative number
determines the index of the neighbor element.

```
element neighbours:
1 3 -1 
2 0 -1 
3 1 -1 
0 2 -1
```

This block is optional. If it isn't given in the macro file, element
neighborships are computed automatically.

== Output ==
Now, the program is started by the call
`./ellipt init/ellipt.dat.2d`. After 9 iterations the tolerance
is reached and the files `output/ellipt.mesh` and
`output/ellipt.dat` are written. In the following Figure the solution and the corresponding mesh is shown. The visualizations was done by the
VTK based tool ''''CrystalClear'''':

[[Bild:ellipt_dat.jpg|thumb|left|400px|Solution after 9 iterations.]]
[[Bild:ellipt_mesh.jpg|thumb|left|400px|corresponding mesh]]

