/*
 <<<<<<<<<<<<<<< Multi-body Dynamics Simualtion, Estimation and Control Lab - McGill University >>>>>>>>>>>>>>
 
 
 Hierarchy     : Cable Driven Parallel Robots
 Model         : 1 DOF Flexible Cable
 Function Name : Example12a
 Function Type : Dynamic Simulation
 
 Description   : The following function contains dynamic sumulation of a
                 single DOF constrained CDPR.
 
 Example-L-0-1 is named using the following convention:
 L - Lumped-mass model type
 0 - Dynamics simulation
 1 - Variable mass/stiffness lumped-mass method.
 
 Revisions: 
 
 2017-09-06: Template creation
 
 Function Details:
 
 Information on addtional functionality can be added here.
 
 References:
 
 [1] Dynamic Modelling and Control of Cable-actuated systems, Harsh Godbole, Master's Thesis, McGill University, 2017.
 
 Armadillo documentation is available at:
 http://arma.sourceforge.net/docs.html
 
 Tempelate by Harsh Godbole.
 Reference credits to Dr. James Richard Forbes and Ryan Caverly.
 
 */

#include <iostream>
#include <armadillo>

//#include <conio.h>

#include <vector>
#include <fstream>
#include <stdio.h>

using namespace std;
using namespace arma;

#ifndef LO1CONSTANTS_H
#define LO1CONSTANTS_H

//1. Structure of Constants
class constants
{
public:
    
    // Exe Globals:
    std::string GenDir;
    // Simulation Model:
	double pi;
	double m_p;
	double m_w;
	double R;
	double J_w;
	double rho;
	double A;
	long double E;
	double L;
	double w;
    double m_rr;
    
	//Damping Constant
	double c_i;
	
	//number of fragments
	int n;
	
	// External forces
	double F_ext;
	
	// Cable Pre-tension:
	double tau_pt_1;
    double tau_pt_2;

	// Adaptive Controller:
    double mu;
	double Gamma_1;
	double Gamma_2;
	double K_d;
	double lambda_k;

    //Load sharing Parameters:
	double c_1;
	double c_2;

	//Test parameters: 

	int rho_d_switch;
	double omega_d;
	double settle_time;
	int post_processing_switch_off;
	int post_processing_switch_on;
    int adaptive_switch;
    int test_switch;
    int feedforward_switch;
    int regressor_switch;
    int disturbance_switch;
    
    // Constructor:
    constants(std::string EnvGenDir);

};

#endif
