#include "Parameters.h"
#include "ClassicalOdeRungeKuttaIntegrators.h"
#include "StabilizedOdeRungeKuttaIntegrators.h"
#include "MultirateOdeProblems.h"
#include <GetPot>
#include <filesystem>

Parameters::Parameters()
{
    problem_size = -1; // must be redefined later when initializind odes, sdes. If not an error will occur
    refsol_path = string("");
}

Parameters::~Parameters()
{
}

bool Parameters::initOde(Ode *&ode)
{
    if (ntest == 1)
        ode = new MultirateDahlquistTestProblem();
    else if (ntest == 2)
        ode = new ScalarNonStiffNonLinearTest();
    else if (ntest == 3)
        ode = new NeuronCable();
    else if (ntest == 4)
        ode = new Brusselator();
    else if (ntest == 5)
        ode = new MultiratePDEBrusselator();
    else if (ntest == 6)
        ode = new Krogh10();
    else if (ntest == 7)
        ode = new PopulationDynamics();
    else if (ntest == 8)
        ode = new VanDerPol();
    else if (ntest == 9)
        ode = new ODEIonicModel();
    else if (ntest == 10)
        ode = new MultirateDiffusionRefinedMesh();
    else if (ntest == 11)
        ode = new MultirateInfectiousDiseaseTransmission();
    else if (ntest == 12)
        ode = new MultirateRobertsonChemicalSystem();
    else if (ntest == 13)
        ode = new Oregonator();
    else if (ntest == 14)
        ode = new MultirateCUSP();
    else if (ntest == 15)
        ode = new MultiratemReactionDiffusion2DEquations();
    else if (ntest == 16)
        ode = new MultirateRadiationDiffusion();
    else if (ntest == 17)
        ode = new MultirateIntegroDifferentialEquation();
    else if (ntest == 18)
        ode = new MultirateMonoDomain();
    else
    {
        cout << "Problem not known" << endl;
        return false;
    }

    if (!std::filesystem::is_directory("../results/"))
        std::filesystem::create_directory("../results/");

    string output_dir = "../results/" + ode->get_problem_name();
    if (!std::filesystem::is_directory(output_dir))
        std::filesystem::create_directory(output_dir);

    output_path = "../results/" + ode->get_problem_name() + "/" + output_file;
    if (refsol_file.compare(string("")) != 0)
        refsol_path = "../results/" + ode->get_problem_name() + "/" + refsol_file;

    problem_size = ode->get_system_size();

    return true;
}

bool Parameters::initOdeTimeIntegrator(OdeRungeKuttaIntegrator *&rk, Ode *ode)
{
    if (rk_name == "EE")
        rk = new ExplicitEuler(this, ode);
    else if (rk_name == "IE")
        rk = new ImplicitEuler(this, ode);
    else if (rk_name == "EM")
        rk = new ExplicitMidpoint(this, ode);
    else if (rk_name == "IM")
        rk = new ImplicitMidpoint(this, ode);
    else if (rk_name == "RK4")
        rk = new RungeKutta4(this, ode);
    else if (rk_name == "RKC1")
        rk = new RKC1(this, ode);
    else if (rk_name == "RKC2")
        rk = new RKC2(this, ode);
    else if (rk_name == "ROCK2")
        rk = new ROCK2(this, ode);
    else if (rk_name == "DROCK2")
        rk = new DROCK2(this, ode);
    else if (rk_name == "RKL1")
        rk = new RKL1(this, ode);
    else if (rk_name == "RKL2")
        rk = new RKL2(this, ode);
    else if (rk_name == "RKU1")
        rk = new RKU1(this, ode);
    else if (rk_name == "RKU2")
        rk = new RKU2(this, ode);
    else if (rk_name == "mRKC")
    {
        if (ode->is_multirate())
            rk = new mRKC(this, dynamic_cast<MultirateOde *>(ode));
        else
        {
            cout << "ERROR: you cant use the multirate integrator " << rk_name << " with the non multirate problem " << ode->get_problem_name() << endl;
            return false;
        }
    }
    else if (rk_name == "IERKC")
    {
        if (ode->is_multirate())
            rk = new IERKC(this, dynamic_cast<MultirateOde *>(ode));
        else
        {
            cout << "ERROR: you cant use the multirate integrator " << rk_name << " with the non multirate problem " << ode->get_problem_name() << endl;
            return false;
        }
    }
    else
    {
        cout << "Integrator " << rk_name << " not known." << endl;
        return false;
    }

    return true;
}

bool Parameters::initOdeIntegration(OdeRungeKuttaIntegrator *&rk, Ode *&ode)
{
    if (!initOde(ode))
        return false;

    return initOdeTimeIntegrator(rk, ode);
}

bool Parameters::read_command_line(int argc, char **argv)
{
    GetPot cl(argc, argv); // command line parser

    ntest = cl.follow(ntest, 2, "-ntest", "-test");
    output_file = cl.follow(output_file.c_str(), 2, "-outputfile", "-ofile");
    refsol_file = cl.follow(refsol_file.c_str(), 3, "-refsol", "-ref", "-refsolfile");
    output_freq = cl.follow(output_freq, 2, "-outputfreq", "-ofreq");
    verbose = cl.follow(verbose, 2, "-verbose", "-verb");

    matlab_output = cl.follow(matlab_output, 3, "-matlab_output", "-matlab_out", "-matlab");
    bin_output = cl.follow(bin_output, 3, "-bin_output", "-bin_out", "-bin");
    specific_output = cl.follow(specific_output, 3, "-spec_output", "-spec_out", "-spec");

    rk_name = cl.follow(rk_name.c_str(), 2, "-solver", "-rk");
    dt = cl.follow(dt, 3, "-dt", "-tau", "-h");
    rho_freq = cl.follow(rho_freq, 3, "-rhofreq", "-rho_freq", "-rfreq");

    if (cl.search("-convtest"))
        conv_test = true;
    else
        conv_test = false;
    max_pow = cl.follow(max_pow, 2, "-maxpow", "-max_pow");
    min_pow = cl.follow(min_pow, 2, "-minpow", "-min_pow");

    dtadap = cl.follow(dtadap, 3, "-dtadaptivity", "-dtadap", "-dta");
    rtol = cl.follow(rtol, 2, "-rtol", "-rt");
    if (!cl.search(2, "-atol", "-at"))
        atol = rtol;
    atol = cl.follow(atol, 2, "-atol", "-at");
    if (!cl.search(2, "-rtol", "-rt"))
        rtol = atol;
    ode_contr = static_cast<Controller>(cl.follow(static_cast<int>(ode_contr), "-oec"));
    err_write_data = cl.follow(err_write_data, "-ewd");

    if (conv_test)
    {
        verbose = false;
        dtadap = false;
        output_freq = -1;
        rho_freq = 1;
    }

    if (cl.search("--help"))
    {
        cout << "This code implements several explicit stabilized methods, among which mRKC, RKC, ROCK2.\n"
             << "Standard methods as EE, IE, RK4, midpoint methods are implemented for comparison.\n"
             << "Different problems are hardcoded in the executable. Look into src/OdeProblems.cpp.\n"
             << "Run the code from the ./build folder as ./MultirateIntegrators OPTIONS, where OPTIONS is a combination of the below.\n"
                "Results are stored in the ./results folder.\n"
             << endl;
        cout << "The following options are available:\n"
             << "--- General options:\n"
             << "    -test       : a number in 1-18 specifying the problem that we want to solve. Default: 1\n"
             << "                  The list of problems is given below.\n"
             << "    -ofile      : name of output file. Default: sol\n"
             << "    -refsol     : name of the reference solution (if available). Default: \"\"\n"
             << "                  If available, errors are computed at the end of the simulation.\n"
             << "    -bin        : Writes solution in binary format. Default: false\n"
             << "    -spec       : Writes solution in specific format for the current ode problem. Default: false\n"
             << "    -matlab     : Writes solution in matlab format. Default: false\n"
             << "    -ofreq      : Output frequency. Default: -1\n"
             << "                  - If 0 writes solution only at the end of simulation.\n"
             << "                    In general, used to generate a reference solution (combine it with '-bin true').\n"
             << "                  - If >0 writes solution every ofreq time steps,\n"
             << "                  - If -1 never writes the solution.\n"
             << "                    In general, used to generate a solution and just compare it against a reference solution.\n"
             << "--- Time integration options:\n"
             << "    -rk         : The name of the numerical integrator to use. Default: RKC1\n"
             << "                  The list of integrators is given below.\n"
             << "    -dt         : Time step size. Default: 1e-2\n"
             << "                  When running with an adaptive step size, this is the initial step size.\n"
             << "    -rfreq      : Frequency at which the spectral radii are re-estimated. Default: 10\n"
             << "    -verb       : Enables or disables verbosity. Default: true\n"
             << "    -dtadap     : Enables/disables time step adaptivity. Default: false\n"
             << "    -atol       : sets the absolute tolerance for error control. Default: 1e-2\n"
             << "    -rtol       : sets the relative tolerance. Default: rtol=atol\n"
             << "    -oec        : sets the error controller. Values are 1,2,3 for\n"
             << "                  I=Integral, PI=Proportional I, PPI= Predictive PI controllers\n"
             << "                  respectively. Default: 3\n"
             << "    -convtest   : If provided, performs a time convergence test, i.e. runs several simulations and checks errors. Default: false\n"
             << "                  If a reference solution is not provided, it is computed on the fly..\n"
             << "    -maxpow     : Minimal step size used for the convergence test is dt=tend/2^maxpow. Default: 6\n"
             << "    -minpow     : Maximal step size used for the convergence test is dt=tend/2^minpow. Default: 3\n"
             << "--- List of problems:\n"
             << " This is the list of problems hardcoded in src/OdeProblems.cpp. You choose them via the -test option."
             << " To change a parameter go to src/OdeProblems.cpp, change it and recompile.\n"
             << " The ones marked with (m) are solvable with both a standard and a multirate solver.\n"
             << "     1 (m) : The Dahlquist test equation.\n"
             << "     2     : A scalar non stiff non linear problem.\n"
             << "     3     : The Neuron cable problem.\n"
             << "     4     : The famous ODE Brusselator benchmark.\n"
             << "     5 (m) : The famous PDE Brusselator benchmark.\n"
             << "     6     : The Krogh10 benchmark.\n"
             << "     7     : A population dynamics problem.\n"
             << "     8     : The famous Van der Pol benchmark.\n"
             << "     9     : Solves a ionic model. By default its Courtemanche1998, but\n"
             << "             HodgkinHuxley1952, LuoRudy1991, and Fox2002 are available as well.\n"
             << "    10 (m) : A heat equation on a non uniform mesh.\n"
             << "    11 (m) : An infectious disease transimission problem.\n"
             << "    12 (m) : The Robertson chemical system benchmark.\n"
             << "    13     : The Oregonator benchmark.\n"
             << "    14 (m) : The CUSP benchmark.\n"
             << "    15 (m) : Four coupled reaction-diffusion problems.\n"
             << "    16 (m) : A radiation-diffusion problem.\n"
             << "    17 (m) : An integro-differential equation.\n"
             << "    18 (m) : The monodomain model in 2D.\n"
             << "--- List of numerical integrators:\n"
             << " Those marked with (m) are multirate integrators and can be applied only to (m) problems.\n"
             << " For implicit solvers, the Jacobians are either provided analytically or \n"
             << " computed numerically (expensive) depending on the problem.\n"
             << "     EE      : Explicit Euler.\n"
             << "     IE      : Implicit Euler.\n"
             << "     EM      : Explicit Midpoint.\n"
             << "     IM      : Implicit Midpoint.\n"
             << "     RK4     : Runge-Kutta 4 (Classical RK).\n"
             << "     IERKC   : Splitting method: stiff term with IE, mildly stiff term with RKC.\n"
             << "     IE      : Implicit Euler.\n"
             << "     mRKC (m): Multirate RKC1 method.\n"
             << "     RKC1    : First order Runge-Kutta-Chebyshev method.\n"
             << "     RKC2    : Second order Runge-Kutta-Chebyshev method.\n"
             << "     RKL1    : First order Runge-Kutta-Legendre method.\n"
             << "     RKL2    : Second order Runge-Kutta-Legendre method.\n"
             << "     RKU1    : First order Runge-Kutta-Chebyshev method with second kind polynomials.\n"
             << "     RKU2    : Second order Runge-Kutta-Chebyshev method with second kind polynomials.\n"
             << "     ROCK2   : Second order stabilized method with optimal polynomials.\n"
             << "     DROKC2  : Damped ROCK2.\n"
             << endl;

        return false;
    }

    return true;
}

void Parameters::print_info()
{
    cout << scientific;

    cout << "--------------- ODE Integration ---------------" << endl;
    cout << "Solver: " << rk_name << "." << endl;
    cout << "Step size: " << dt << "." << endl;
    cout << "Step size adaptivity: " << (dtadap ? "yes." : "no.") << endl;
    if (dtadap)
    {
        cout << "Controller type: ";
        if (ode_contr == I)
            cout << "I." << endl;
        else if (ode_contr == PI)
            cout << "PI." << endl;
        else if (ode_contr == PPI)
            cout << "PPI." << endl;
        cout << "Relative tolerance: " << rtol << "." << endl;
        cout << "Absolute tolerance: " << atol << "." << endl;
        cout << "Write controller data: " << (err_write_data ? "yes." : "no.") << endl;
    }
    if (conv_test)
    {
        cout << "Convergence test parameters: min_pow = " << min_pow << ", max_pow = " << max_pow << "." << endl;
    }
    else
    {
        cout << "Output file name: " << output_path << endl;
        cout << "Output frequency: " << output_freq << endl;
        cout << "Verbose: " << (verbose ? "yes." : "no.") << endl;
    }
    cout << "-----------------------------------------------" << endl;
}

void Parameters::print_info(Ode *ode)
{
    cout << scientific;

    cout << "----------------- ODE Problem -----------------" << endl;
    cout << "Problem name: " << ode->get_problem_name() << "." << endl;
    cout << "Problem size: " << ode->get_system_size() << "." << endl;
    cout << "Constant rho: " << (ode->is_rho_constant() ? "yes." : "no.") << endl;
    cout << "Known rho estimation: " << (ode->estimation_rho_known() ? "yes." : "no.") << endl;
    cout << "-----------------------------------------------" << endl;
}
