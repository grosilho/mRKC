#ifndef INIT_H
#define	INIT_H

#include "MainHeader.h"

class TimeIntegrator;
class OdeRungeKuttaIntegrator;
class Ode;

class GetPot;

class Parameters
{
public:
    Parameters();
    ~Parameters();
    
    bool read_command_line(int argc, char** argv);
    
    bool initOde(Ode*& ode);
    bool initOdeTimeIntegrator(OdeRungeKuttaIntegrator*& rk, Ode* ode);
    bool initOdeIntegration(OdeRungeKuttaIntegrator*& rk, Ode*& ode);
    
    void print_info();
    void print_info(Ode* ode);
    
public:
    int ntest;
    string output_path;
    string refsol_path;
    string output_file;
    string refsol_file;
    int output_freq;
    bool verbose;
    
    bool matlab_output;
    bool bin_output;
    bool specific_output;
    
    string rk_name;
    Real dt;
    unsigned int rho_freq;
    
    bool conv_test;
    unsigned int max_pow;
    unsigned int min_pow;
    
    bool dtadap;
    Real rtol;
    Real atol;
    Controller ode_contr;
    bool err_write_data;
    
    int problem_size;
};


#endif	/* INIT_H */

