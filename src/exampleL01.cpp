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
 
 Tempelate by Harsh Godbole. 
 Reference credits to Dr. James Richard Forbes and Ryan Caverly.
 */
//Tempelate by Harsh Godbole. Reference credits Dr. James Richard Forbes

// Generic includes:
#include <iostream>
#include <armadillo>
#include <vector>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>

// Project Specific includes:
#include "L01Constants.hpp"
#include "L01ODE.hpp"
#include "RK4.hpp"
#include "L01PostProcessing.hpp"


using namespace std;
using namespace arma;

//17. Main Function
int main(int argc, char** argv) {
    
    /*
     ------------------------------
     Step 01 : Initialize ODE
     ------------------------------
     */
    
    // Get generated directory path and dump files in correct location
    std::string RelGenPath("/Generated/");
    std::string ProjectHomeDir("");
    
    // Check that getenv has not returned a NULL pointer:
    const char *ProjectHomeDirPtr = getenv("PROJECT_HOME_DIR");
    if (!ProjectHomeDirPtr)
    {
         ProjectHomeDir.assign("//");
         cout<<"Env Variable PROJECT_HOME_DIR is not set, Use L01Runfile to execute code"<<endl;
    }
    else
    {
        ProjectHomeDir.assign(getenv("PROJECT_HOME_DIR"));
    
        //std::string ProjectHomeDir = getenv("PROJECT_HOME_DIR");
        std::string GenDir(ProjectHomeDir + RelGenPath);
    
        // Initializing constants parameters class
        constants cst(GenDir);
    
        // ODE: Initial conditions
        int x_length_1;
        x_length_1 = 4 * (cst.n) + 9;
        
        vec x_IC=zeros<vec>(x_length_1);

        //(a) Generalized co-ordinates
    
        x_IC(1) = 0.18; //theta_1
    
        x_IC(2) = 0.00144; //n+1 states corrosponding to x_n1 and x_p1
        x_IC(3) = 0.00144;
        x_IC(4) = 0.00144;
        x_IC(5) = 0.00144;
        x_IC(6) = 0.00144;
    
        x_IC(7) = 0.00144; //n+1 states corrosponding to x_n2 and x_p2
        x_IC(8) = 0.00144;
        x_IC(9) = 0.00144;
        x_IC(10) = 0.00144;
        x_IC(11) = 0.00144;
    
        x_IC(0) = (-cst.R*x_IC(1)+ sum(x_IC.subvec(2,11)))/cst.R; //theta_2 by constraint equation to keep rho = 0 initially.
    
        //(b) Generalized velocities
    
    
        x_IC(12) = 0.00; //theta_1_dot
        
        x_IC(13) = 0.00; //x_n1_dot and x_p1_dot
        x_IC(14) = 0.00;
        x_IC(15) = 0.00;
        x_IC(16) = 0.00;
        x_IC(17) = 0.00;
    
        x_IC(18) = 0.00; //x_n2_dot and x_p2_dot
        x_IC(19) = 0.00;
        x_IC(20) = 0.00;
        x_IC(21) = 0.00;
        x_IC(22) = 0.00;

        //Adaptive controller states:
    
        //x_IC(23) = 0.05;
        x_IC(23) = (((cst.J_w/(cst.R*cst.R))+((cst.rho*cst.A*cst.L)) + (cst.m_p))*2)*0.1;  //M_rho_cap: Intial estimate is 10 percent of actual
	    x_IC(24) = 0.0;  //d_cap = G(rho)
    
        /*
         -------------------------------------
         Step 02 : Set Simulation Parameters:
         -------------------------------------
         */
        
	    double
	    t_start = 0,          // Start time in seconds
	    t_end   = 4,          // End time in seconds
	    h       = 0.0001;     // Step Size in seconds
    
	    long long int a = (t_end - t_start) / h; //Need better counting Subroutine
	    cout << "Number of Iterations Estimated "<<a<<endl;
	    mat x_out=zeros<mat>(x_length_1, a);
	                                                    //vec x_out(x_length_1);//test
    
        /*
         -------------------------------------
         Step 03 : Solve ODE
         -------------------------------------
         */
    
	    cout<< "RK4 process started"<<endl;
    
	    x_out = RK4(h, t_start, t_end, x_IC, cst, ODE);
	                                                    //x_out=ODE(x_IC,1,cst);//test
    
	    cout << "RK4 process complete"<<endl;
	    
    
         /*
         -------------------------------------
         Step 04 : Postprocessing Results
         -------------------------------------
         */
    
	    // Lagrange multipluer calculation option:
    
         int calc_lagrange_multiplier = 1; //  1 -> on
                                           //  0 -> off
      
    
	    // Post processing and file export
        
        cout<< "Post processing started"<<endl;
    
        post_processing(h,t_start,t_end,x_out,calc_lagrange_multiplier,cst);
    
        cout<< "Post processing completed"<<endl;
     
        
        // **** Plot-scripts in matlab:
        // **** 1. system() Directly execute test_script.m and plot script commands
        
        //system("/Applications/MATLAB_R2016a.app/bin/matlab -nodesktop -r \"run test_script.m\"");
    } // end if (!ProjectHomeDirPtr)

	return 0;
}
// End of Main
