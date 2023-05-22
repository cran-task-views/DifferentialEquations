---
name: DifferentialEquations
topic: Differential Equations
maintainer: Thomas Petzoldt, Karline Soetaert
email: thomas.petzoldt@tu-dresden.de
version: 2023-05-22
source: https://github.com/cran-task-views/DifferentialEquations/
---


Differential equations (DE) are mathematical equations that describe how
a quantity changes as a function of one or several (independent)
variables, often time or space. Differential equations play an important
role in biology, chemistry, physics, engineering, economy and other
disciplines.

Differential equations can be separated into stochastic versus
deterministic DEs. Problems can be split into initial value problems
versus boundary value problems. One also distinguishes ordinary
differential equations from partial differential equations, differential
algebraic equations and delay differential equations. All these types of
DEs can be solved in R. DE problems can be classified to be either stiff
or nonstiff; the former type of problems are much more difficult to
solve.

The [dynamic models SIG](https://stat.ethz.ch/mailman/listinfo/r-sig-dynamic-models)
is a suitable mailing list for discussing the use of R for solving
differential equation and other dynamic models such as individual-based
or agent-based models.

This task view was created to provide an overview on the topic. If something
is missing, or if a new package should be mentioned here, please e-mail the
maintainers or submit an issue or pull request in the GitHub
repository linked above.

### Stochastic Differential Equations (SDEs)

In a stochastic differential equation, the unknown quantity is a
stochastic process.

-   The package `r pkg("sde", priority = "core")` provides
    functions for simulation and inference for stochastic differential
    equations. It is the accompanying package to the book by Iacus
    (2008).
-   The package `r pkg("pomp")` contains functions for
    statistical inference for partially observed Markov processes.
-   Packages `r pkg("adaptivetau")` and
    `r pkg("GillespieSSA")` implement Gillespie's "exact"
    stochastic simulation algorithm (direct method) and several
    approximate methods.
-   The package `r pkg("resde")` computes maximum likelihood parameter
    estimates for univariate reducible stochastic differential equation
    models.
-   The package `r pkg("Sim.DiffProc")` provides functions
    for simulation of Itô and Stratonovitch stochastic differential
    equations.
-   Package `r pkg("diffeqr")` can solve SDE problems using
    the **DifferentialEquations.jl** package from the Julia programming
    language.

### Ordinary Differential Equations (ODEs)

In an ODE, the unknown quantity is a function of a single independent
variable. Several packages offer to solve ODEs.

-   The "odesolve" package was the first to solve ordinary
    differential equations in R. It contained two integration methods.
    It has been replaced by the package
    `r pkg("deSolve", priority = "core")`.
-   The package `r pkg("deSolve")` contains several solvers
    for solving ODE, DAE, DDE and PDE. It can deal with stiff and
    nonstiff problems.
-   The package `r pkg("odeintr")` generates and compiles
    C++ ODE solvers on the fly using Rcpp and
    [Boost](http://www.boost.org/) [odeint](http://www.odeint.com/) .
-   The R package `r pkg("diffeqr")` provides a seamless
    interface to the **DifferentialEquations.jl** package from the Julia
    programming language. It has unique high performance methods for
    solving ODE, SDE, DDE, DAE and more. Models can be written in either
    R or Julia. It requires an installation of the Julia language.
-   Package `r pkg("pracma")` implements several adaptive
    Runge-Kutta solvers such as ode23, ode23s, ode45, or the
    Burlisch-Stoer algorithm to obtain numerical solutions to ODEs with
    higher accuracy.
-   Package `r pkg("rODE")` (inspired from the book of
    Gould, Tobochnik and Christian, 2016) aims to show physics, math and
    engineering students how ODE solvers can be made with R's S4
    classes.
-   Package `r pkg("sundialr")` provides a way to call the
    'CVODE' function from the 'SUNDIALS' C ODE solving library. The
    package requires the ODE to be written as an 'R' or 'Rcpp'
    function.
-   The package `r pkg("mrgsolve")` compiles ODEs on the fly
    and allows shorthand prescription dosing.
-   The package `r pkg("rxode2")` is similar to
    `r pkg("mrgsolve")`, but has the added value of being
    the backend of the nonlinear mixed effects modeling R package
    `r pkg("nlmixr2")`.

### Delay Differential Equations (DDEs)

In a DDE, the derivative at a certain time is a function of the variable
value at a previous time.

-   The `r pkg("dde")` package implements solvers for
    ordinary (ODE) and delay (DDE) differential equations, where the
    objective function is written in either R or C. Suitable only for
    non-stiff equations. Support is also included for iterating
    difference equations.
-   The package `r pkg("PBSddesolve")` (originally published
    as "ddesolve") includes a solver for non-stiff DDE problems.
-   Functions in the package `r pkg("deSolve")` can solve
    both stiff and non-stiff DDE problems.
-   Package `r pkg("diffeqr")` can solve DDE problems using
    the **DifferentialEquations.jl** package from the Julia programming
    language.

### Partial Differential Equations (PDEs)

PDEs are differential equations in which the unknown quantity is a
function of multiple independent variables. A common classification is
into elliptic (time-independent), hyperbolic (time-dependent and
wavelike), and parabolic (time-dependent and diffusive) equations. One
way to solve them is to rewrite the PDEs as a set of coupled ODEs, and
then use an efficient solver.

-   The R-package `r pkg("ReacTran")` provides functions for
    converting the PDEs into a set of ODEs. Its main target is in the
    field of ''reactive transport'' modelling, but it can be used to
    solve PDEs of the three main types. It provides functions for
    discretising PDEs on cartesian, polar, cylindrical and spherical
    grids.
-   The package `r pkg("deSolve")` contains dedicated
    solvers for 1-D, 2-D and 3-D time-varying ODE problems as generated
    from PDEs (e.g. by `r pkg("ReacTran")`).
-   The package `r pkg("rootSolve", priority = "core")`
    contains optimized solvers for 1-D, 2-D and 3-D algebraic problems
    generated from (time-invariant) PDEs. It can thus be used for
    solving elliptic equations.

Note that, to date, PDEs in R can only be solved using finite
differences. At some point, we hope that finite element and spectral
methods will become available.

### Differential Algebraic Equations (DAEs)

Differential algebraic equations comprise both differential and
algebraic terms. An important feature of a DAE is its differentiation
index; the higher this index, the more difficult to solve the DAE.

-   The package `r pkg("deSolve")` provides two solvers,
    that can handle DAEs up to index 3.
-   Package `r pkg("diffeqr")` can solve DAE problems using
    the **DifferentialEquations.jl** package from the Julia programming
    language.

### Boundary Value Problems (BVPs)

BVPs have solutions and/or derivative conditions specified at the
boundaries of the independent variable.

-   The package `r pkg("ReacTran")` can solve BVPs that
    belong to the class of reactive transport equations.
-   Package `r pkg("diffeqr")` can also solve BVPs using the
    **DifferentialEquations.jl** package from the Julia programming
    language.

### Population ODE modeling

-   The package `r pkg("nlmixr2")` fits ODE-based nonlinear
    mixed effects models using `r pkg("rxode2")`.

### Other

-   The `r pkg("simecol")` package provides an interactive
    environment to implement and simulate dynamic models. Next to DE
    models, it also provides functions for grid-oriented,
    individual-based, and particle diffusion models.
-   In the package `r pkg("FME")` are functions for inverse
    modelling (fitting to data), sensitivity analysis, identifiability
    and Monte Carlo Analysis of DE models.
-   `r pkg("mkin")` provides routines for fitting kinetic
    models with one or more state variables to chemical degradation
    data.
-   Package `r pkg("dMod")` provides functions to generate
    ODEs of reaction networks, parameter transformations, observation
    functions, residual functions, etc. It follows the paradigm that
    derivative information should be used for optimization whenever
    possible.
-   The package `r pkg("CollocInfer")` implements
    collocation-inference for continuous-time and discrete-time
    stochastic processes.
-   Root finding, equilibrium and steady-state analysis of ODEs can be
    done with the package `r pkg("rootSolve")`.
-   The `r pkg("PBSmodelling")` package adds GUI functions
    to models.
-   Package `r pkg("cOde")` supports the automatic creation
    of dynamically linked code for packages
    `r pkg("deSolve")` (or a built-in implementation of the
    sundials cvode solver) from inline C embedded in the R code.
-   Package `r pkg("rodeo")` is an object oriented system
    and code generator that creates and compiles efficient Fortran code
    for `r pkg("deSolve")` from models defined in
    stoichiometry matrix notation.
-   Package `r pkg("ecolMod")` contains the figures, data
    sets and examples from a book on ecological modelling (Soetaert and
    Herman, 2009).



### Links
-   Wikipedia: [Differential equation](http://en.wikipedia.org/wiki/Differential_equation)
-   R-Forge website: [deSolve](https://deSolve.R-Forge.R-project.org) (differential equation solvers)
-   Github website: [pomp](https://kingaa.github.io/pomp/) (partially observed Markov process)
-   Book: [Iacus, SM. 2008. Simulation and Inference for Stochastic Differential Equations: with R examples, Springer](http://www.springer.com/978-0-387-75838-1)
-   Book: [Soetaert, K. and P.M.J. Herman, 2009. A Practical Guide to Ecological Modelling, using R as a simulation Platform, Springer.](http://www.springer.com/life+sciences/ecology/book/978-1-4020-8623-6)
-   Book: [Stevens, H, 2009. A Primer of Ecology with R, Springer](http://www.springer.com/life+sci/ecology/book/978-0-387-89881-0) and the [2021 online edition](https://hankstevens.github.io/Primer-of-Ecology/) on Github.
-   Book: [Soetaert, K., Cash, J. and Mazzia, F. 2012. Solving Differential Equations in R, Springer.](http://www.springer.com/statistics/computanional+statistics/book/978-3-642-28069-6)
