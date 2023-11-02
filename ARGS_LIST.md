```
./MultirateIntegrators --help
This code implements several explicit stabilized methods, among which mRKC, RKC, ROCK2.
Standard methods as EE, IE, RK4, midpoint methods are implemented for comparison.
Different problems are hardcoded in the executable. Look into src/OdeProblems.cpp.
Run the code from the ./build folder as ./MultirateIntegrators OPTIONS, where OPTIONS is a combination of the below.
Results are stored in the ./results folder.

The following options are available:
--- General options:
    -test       : a number in 1-18 specifying the problem that we want to solve. Default: 1
                  The list of problems is given below.
    -ofile      : name of output file. Default: sol
    -refsol     : name of the reference solution (if available). Default: ""
                  If available, errors are computed at the end of the simulation.
    -bin        : Writes solution in binary format. Default: false
    -spec       : Writes solution in specific format for the current ode problem. Default: false
    -matlab     : Writes solution in matlab format. Default: false
    -ofreq      : Output frequency. Default: -1
                  - If 0 writes solution only at the end of simulation.
                    In general, used to generate a reference solution (combine it with '-bin true').
                  - If >0 writes solution every ofreq time steps,
                  - If -1 never writes the solution.
                    In general, used to generate a solution and just compare it against a reference solution.
--- Time integration options:
    -rk         : The name of the numerical integrator to use. Default: RKC1
                  The list of integrators is given below.
    -dt         : Time step size. Default: 1e-2
                  When running with an adaptive step size, this is the initial step size.
    -rfreq      : Frequency at which the spectral radii are re-estimated. Default: 10
    -verb       : Enables or disables verbosity. Default: true
    -dtadap     : Enables/disables time step adaptivity. Default: false
    -atol       : sets the absolute tolerance for error control. Default: 1e-2
    -rtol       : sets the relative tolerance. Default: rtol=atol
    -oec        : sets the error controller. Values are 1,2,3 for
                  I=Integral, PI=Proportional I, PPI= Predictive PI controllers
                  respectively. Default: 3
    -convtest   : If provided, performs a time convergence test, i.e. runs several simulations and checks errors. Default: false
                  If a reference solution is not provided, it is computed on the fly..
    -maxpow     : Minimal step size used for the convergence test is dt=tend/2^maxpow. Default: 6
    -minpow     : Maximal step size used for the convergence test is dt=tend/2^minpow. Default: 3
--- List of problems:
 This is the list of problems hardcoded in src/OdeProblems.cpp. You choose them via the -test option. To change a parameter go to src/OdeProblems.cpp, change it and recompile.
 The ones marked with (m) are solvable with both a standard and a multirate solver.
     1 (m) : The Dahlquist test equation.
     2     : A scalar non stiff non linear problem.
     3     : The Neuron cable problem.
     4     : The famous ODE Brusselator benchmark.
     5 (m) : The famous PDE Brusselator benchmark.
     6     : The Krogh10 benchmark.
     7     : A population dynamics problem.
     8     : The famous Van der Pol benchmark.
     9     : Solves a ionic model. By default its Courtemanche1998, but
             HodgkinHuxley1952, LuoRudy1991, and Fox2002 are available as well.
    10 (m) : A heat equation on a non uniform mesh.
    11 (m) : An infectious disease transimission problem.
    12 (m) : The Robertson chemical system benchmark.
    13     : The Oregonator benchmark.
    14 (m) : The CUSP benchmark.
    15 (m) : Four coupled reaction-diffusion problems.
    16 (m) : A radiation-diffusion problem.
    17 (m) : An integro-differential equation.
    18 (m) : The monodomain model in 2D.
--- List of numerical integrators:
 Those marked with (m) are multirate integrators and can be applied only to (m) problems.
 For implicit solvers, the Jacobians are either provided analytically or
 computed numerically (expensive) depending on the problem.
     EE      : Explicit Euler.
     IE      : Implicit Euler.
     EM      : Explicit Midpoint.
     IM      : Implicit Midpoint.
     RK4     : Runge-Kutta 4 (Classical RK).
     IERKC   : Splitting method: stiff term with IE, mildly stiff term with RKC.
     IE      : Implicit Euler.
     mRKC (m): Multirate RKC1 method.
     RKC1    : First order Runge-Kutta-Chebyshev method.
     RKC2    : Second order Runge-Kutta-Chebyshev method.
     RKL1    : First order Runge-Kutta-Legendre method.
     RKL2    : Second order Runge-Kutta-Legendre method.
     RKU1    : First order Runge-Kutta-Chebyshev method with second kind polynomials.
     RKU2    : Second order Runge-Kutta-Chebyshev method with second kind polynomials.
     ROCK2   : Second order stabilized method with optimal polynomials.
     DROKC2  : Damped ROCK2.
   
```