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


#include "L01Constants.hpp"

//1. Structure of Constants
constants::constants(std::string EnvGenDir)
{
    // Global Home Path for generated files:
    GenDir = EnvGenDir;
    
    // Numerical Constants
    pi=3.14159265;
    
    // Physical variables
    m_p = 0.5 ;//* 1.5 ;//* 1.6;
    m_w = 0.5;
    R = 0.04 ;//* 1.2;
    J_w = 1.39*0.00001;
    rho = 2200 ;//* 0.2;
    //rho = 7.8;
    
    
    A = 17.95*0.000001;
    //A = 17.95*0.001;
    
    E = 500 * 1000000 ;//* (0.25);
    //E = 180 * 1000000;// Steel found online
    
    L = 0.5;
    //L=0.5;
    
    w =2;
    
    c_i=0.00;
    
    // Corrosponds to the number of rounds of cable already wound
    /*
    m_p = 0.5;
    R = 2;
    rho = 5;
    A = 1;
    E = 1000;
    L = 5;
    w = 2;
    //Corrosponds to the number of rounds of cable already wound
    */
    
    
    J_w = 1.39 * 0.00001 ;//* 1.5; // + rho*A*(R*R*R)*6.28 * w;
    
    // Options:
    // J_w = 0.5*m_w*R*R + rho*A*(R*R*R)*6.28 * w;
    // J_w =  rho*A*(R*R*R)*6.28 * w;
    
    // J_w = 1.39*0.00001 + rho*A*(R*R*R)*6.28 * w;
    
    
    
    //Number of fragments
    n = 4;
    
    // Constant External forces, if any
    F_ext = 0;
    
    
    // Cable Pre-tension:
    tau_pt_1 = 5.1136364;
    tau_pt_2 = 5.1136364;
    
    // Feedforward term:
    m_rr = (((J_w/(R*R))+((rho*A*L)) + (m_p))*2);
    
    
    //Control Parameters:
    
    // mu-tip constant:
    mu=0.9;
    
    // PD coefficients:
    //K_p=500;
    //K_d=60;
    
    
    //K_p=120;
    //K_d=30;
    
    // Options:
    // K_p=1.6;
    // K_d=5;
    
    
    
    
    // Load sharing parameters:
    
    c_1=0.5;
    c_2=0.5;                // Note: c_1+c_2=1
    
    // Test function switch:
    
    // #1 Options: 1. Step,  2. Steady state, 3. parabolic infinite, 4. Sine oscillations, 5. Steady state 0.2L
    // #2 Value of Omega for  option 3,
    // #5 0. Fixed mass Feedforward 1. Adaptive controller
    
    // #6 0. Unperturbed mass and stiffness 1. Perturbed mass and stiffness
    
    // #7 0. fixed mass FF off 1. Fixed mass FF off
    
    // #8 1. Regressor(2) constant 2. Regreser(2) with \theta
    
    rho_d_switch                = 1;  // #1
    omega_d                     = 2;  // #2
    settle_time                 = 0;  // #3
    post_processing_switch_off  = 0;
    post_processing_switch_on   = 1;  //
    adaptive_switch             = 1;  // #5
    test_switch                 = 0;  // #6
    feedforward_switch          = 1;  // #7
    regressor_switch            = 1;  // #8
    disturbance_switch          = 0;  //


    
    
    // Model Pertubations:
    if(test_switch == 1)
    {
        m_p = 0.5 * 1.3 ;
        
        rho = 2200 * 0.5;
        
        E = 500 * 1000000 * (0.5);
        
    }
    
    //Adaptive controller tuning:
    
    //Generic Tuning Values:
    Gamma_1   = 500;
    Gamma_2   = 1000;
    
    //lambda_k  = 0.25;
    //K_d       = 30;
    
    lambda_k  = 4.4;
    K_d       = 150;
    
    
    if(rho_d_switch == 1)
    {
        Gamma_1=500;
        
        if(regressor_switch == 2)
            Gamma_2 = 650;
        else
            Gamma_2=40500;
        // Gamma_2 = 5000;
        
    }
    
    if(rho_d_switch == 4)
    {
        Gamma_1=500;
        
        if(regressor_switch == 2)
            Gamma_2=550;
        else
            Gamma_2=1000;
        //Gamma_2=550;
    }
    
    
    //Gain collection references:
    
    //For Steady signals
    
    // 1. Older Data
    //K_d=200;
    //Gamma_1=3;
    //Gamma_2=2;
    //lambda_k=4;
    
    // 2. Latest data
    
    //K_d=350;
    //Gamma_1=100;
    //Gamma_2=100;
    //lambda_k=4.5;
    
    //3. Best performance
    
    //K_d=300;
    //Gamma_1=100;
    //Gamma_2=100;
    //lambda_k=0.25;
    
    //For oscillatory disturbances:
    
    //K_d=150;
    //Gamma_1=1000;
    //Gamma_2=1000;
    //lambda_k=4.4;
}

