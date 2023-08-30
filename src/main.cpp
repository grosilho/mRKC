#include "MainHeader.h"
#include "Parameters.h"
#include "OdeRungeKuttaIntegrator.h"


int main(int argc, char** argv)
{
  
//  Default parameters
    Parameters params; // Class for storing the parameters/options and initializing the integrators
    params.ntest = 20;           //problem number
    params.refsol_path = "";    //name of file containing reference solution. If empty, in a convergence test a reference solution is computed on the fly.
    params.verbose = true;      //show step info, otherwise just initial and final info.
    params.eq = ODE;            //ODE or D_SDE (SDE driven by diffusion) or JD_SDE (SDE driven by jump diffusion).
    params.rk_name = "RKC1";     //name of the solver
    params.dt = 1e-2;            //step size, or starting step size in case dtadap=true.
    params.rho_freq = 10;        //estimation of the spectral radius every rho_freq steps
    
    params.output_file = "sol"; //output file name
    params.output_freq = -1;    //-1 means no outputs, 0 just at the end, otherwise every output_freq steps
    params.matlab_output = false; // generates the .m file or not
    params.bin_output = false;    //generates the .bin file or not
    params.specific_output = false;//calls a problem specific output function
    
//  Parameters for convergence experiments
    params.conv_test = false;   //if false, we run once the experiment with step size params.dt
    params.max_pow = 7;         //if true we run a convergence test with step sizes tend/2^k,
    params.min_pow = 6;         //with k=min_pow,...,max_pow.
    
//  Specific for step size adaptivity (only for ODEs)
    params.dtadap = false;      //toggles step size adaptivity.
    params.rtol = 1e-2;         //relative tolerance
    params.atol = 1e-2;         //absolute tolerance
    params.ode_contr = PPI;     //type of error controller: I (integral) or PPI (predictive proportional integral)
    params.err_write_data = false; //write error control data, as estimated errors, step sizes, predicted errors
    
    
    //read input parameters and eventually resolve incompatible choices. Exit if major error.
    if(!params.read_command_line(argc, argv))
        return 0;
    
    params.print_info();
        
    Ode* ode=0;
    OdeRungeKuttaIntegrator* integrator=0;

    if(!params.initOdeIntegration(integrator,ode))
        return 0;

    params.print_info(ode);

    if(params.conv_test)
        integrator->convergence_test();
    else
    {
        integrator->integrate();
        integrator->print_integration_info();
    }

    delete ode;
    delete integrator;
    
    return 0;
}
