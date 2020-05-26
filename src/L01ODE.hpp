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

#include <iostream>
#include <armadillo>

//#include <conio.h>

#include <vector>
#include <fstream>
#include <stdio.h>
#include "L01Constants.hpp"
#include "L01MassMat.hpp"
#include "L01StiffnessMat.hpp"
#include "L01ControlLaws.hpp"

using namespace std;
using namespace arma;


//Tempelate by Harsh Godbole. Reference credits Dr. James Richard Forbes


// Armadillo documentation is available at:
// http://arma.sourceforge.net/docs.html

//13. ODE function: f'(x)=ODE(x)
vec ODE(vec x, double t, constants cst)
{

	//Extracting states and Preprocessing

	// (a) Declarations
	int n;
	n = cst.n;
	
	vec x_dot=zeros<vec>(4*n+9);
	
	double theta_2,theta_2_dot,theta_1, theta_1_dot, m_i_1, m_i_2;
	vec q_e_1(n+1);
    vec q_e_2(n+1);
    vec q_e_1_dot(n+1);
    vec q_e_2_dot(n+1);

	// (b) Initializations
	
	theta_2 = x(0);//flag
	theta_1 = x(1);//flag
    q_e_1=x.subvec(2,n+2);
    q_e_2=x.subvec(n+3,(2*n+3));

	theta_1_dot = x(2*n + 4);//flag
	q_e_1_dot=x.subvec((2*n+5),(3*n+5));
    q_e_2_dot=x.subvec((3*n+6),(4*n+6));

    // Adaptive Controller state

   	vec a_cap(2);

   	a_cap = x.subvec(4*n+7,4*n+8);


	// Defining Transformation matrix phi
	mat phi=zeros<mat>(2*n+4,2*n+4);
	phi(0,0)=1;
	phi(n+2,1)=1;
	for(int i=0;i<n+1;i++)
	    { 
		    for(int j=0;j<n+1;j++)
			  { if(i==j)
			     phi(i+1,i+2)=1;
			     phi(i+n+3,i+n+3)=1;
			  }
	    }

    //Defining Reduction matrix R (upsilon)
	    mat R=zeros<mat>(2*n+4,2*n+3);

	    double J_theta_1,J_theta_2;
	    J_theta_1=-cst.R;
	    J_theta_2=cst.R;

	    vec J_e_1_trans=ones<vec>(n+1);
	    vec J_e_2_trans=ones<vec>(n+1);
	    J_e_2_trans=(-1)*J_e_2_trans;

	    R(0,0)=1;
        R(1,0)=J_theta_1/J_theta_2;

	    for(int i=0;i<n+1;i++)
	    {
	    	R(1,i+1)=(J_e_1_trans(i))/J_theta_2;
	    	R(1,i+n+2)=-(J_e_2_trans(i))/J_theta_2;
	    }
	    
	    for(int i=0;i<(2*n+2);i++)
	    	{
	    		for(int j=0;j<(2*n+2);j++)
	    	    {
	    		   if(i==j)
	    			   R(i+2,i+1)=1;
	    	    }
	    	}
    

	//First stage of integraton:
	vec theta_2_dot_vec(1);

	theta_2_dot_vec=(((J_theta_1/J_theta_2)*theta_1_dot)+ (((trans(J_e_1_trans)*q_e_1_dot)/J_theta_2))-(((trans(J_e_2_trans)*q_e_2_dot)/J_theta_2)));
    theta_2_dot=theta_2_dot_vec(0);

    vec z(2*n+3);
    vec z_dot(2*n+3);
    
    z(0)=theta_1;
    z.subvec(1,n+1)=q_e_1;
    z.subvec(n+2,2*n+2)=q_e_2;

    z_dot(0)=theta_1_dot;
    z_dot.subvec(1,n+1)=q_e_1_dot;
    z_dot.subvec(n+2,2*n+2)=q_e_2_dot;

    //Filling values of First stage into x_dot
    x_dot(0)=theta_2_dot;
    x_dot.subvec(1,2*n+3)=z_dot;


	//Second Stage of Integration:

	//Defining Kinematic variables:
	vec q_1=zeros<vec>(n+2);
	vec q_2=zeros<vec>(n+2);
	vec nu_1=zeros<vec>(n+2);
	vec nu_2=zeros<vec>(n+2);

	q_1(0)=theta_1;
	q_1.subvec(1,n+1)=q_e_1;
	q_2(0)=theta_2;
	q_2.subvec(1,n+1)=q_e_2;

	nu_1(0)=theta_1_dot;
	nu_1.subvec(1,n+1)=q_e_1_dot;
	nu_2(0)=theta_2_dot;
	nu_2.subvec(1,n+1)=q_e_2_dot;
    
                                                          //q_1.print("q_1=");
	                                                      //q_2.print("q_2=");
	                                                      //nu_1.print("nu_1=");
	                                                      //nu_2.print("nu_2=");


	//Mass variables:
	m_i_1 = cst.rho*cst.A*(cst.L - cst.R*theta_1) / n;       
    m_i_2 = cst.rho*cst.A*(cst.L - cst.R*theta_2) / n;      
                                                          //cout << "m_i_1 = " << m_i_1 << endl;//test
                                                          //cout << "m_i_2 = " << m_i_2 << endl;//test


	//(a) Non Linear Components:
   
    vec f_non=zeros<vec>(2*n+4);
	vec f_non_zz=zeros<vec>(2*n+3);
	
	f_non.subvec(0,n+1) = -(((M_dot(n, theta_1_dot, cst))*nu_1) - (Partial_M(n, nu_1, cst)) + (Partial_k(n, theta_1, q_1, cst)));
	f_non.subvec(n+2,2*n+3)= -(((M_dot(n, theta_2_dot, cst))*nu_2) - (Partial_M(n, nu_2, cst)) + (Partial_k(n, theta_2, q_2, cst)));
  
                                                         //f_non.print("f_non");//test
                                                         //R.print("R:");
                                                         //phi.print("phi:");


    f_non_zz=(((trans(R))*(trans(phi)))*f_non);
                                                         //f_non_zz.print("f_non_zz = ");//test


	//(b) Calculating Mass Matrix and Stiffness Matrices
    
    //Mass Matrix:
	mat M = zeros<mat>(2*n + 4, 2*n + 4);
	mat M_zz = zeros<mat>(2*n + 3, 2*n + 3);
	

	M.submat(0,0,n+1,n+1) = mass_matrix(n, m_i_1, theta_1, cst);
	M.submat(n+2,n+2,2*n+3,2*n+3) = mass_matrix(n, m_i_2, theta_2, cst);
	
	M_zz=(trans(R))*(trans(phi))*M*phi*R;
                                                          //M_zz.print("M_zz = ");//test
    //Stiffness Matrix:
	
	mat K = zeros<mat>(2*n + 4, 2*n + 4);
	mat K_zz = zeros<mat>(2*n + 3, 2*n + 3);
	
	K.submat(0,0,n+1,n+1) = K_matrix(n, theta_1, cst);
	K.submat(n+2,n+2,2*n+3,2*n+3) = K_matrix(n, theta_2, cst);
    
	K_zz=(trans(R))*(trans(phi))*K*phi*R;
                                                          //K_zz.print("K_zz = ");//test
    	                             

	//(c) Damping Matrix:

	mat D = zeros<mat>(2*n + 4, 2*n + 4);
	mat D_zz = zeros<mat>(2*n + 3, 2*n + 3);
	
	D=damp(2*n+2,cst);
	D_zz=(trans(R))*(trans(phi))*D*phi*R;
	                                                      //D_zz.print("D_zz = ");//test
    


	//(d) External forces and Control Laws
	
	vec b_1_cap=zeros<vec>(n+2);
	vec b_2_cap=zeros<vec>(n+2);
	b_1_cap(0) = 1;
	b_2_cap(0) = 1;
    
    if(cst.disturbance_switch == 1)
    {
        b_1_cap(n+1)  = 1.5 * t + 1.2 * sin(t);
        b_2_cap(n+1)  = 1.5 * t + 1.2 * cos(t);
    }
	
	mat B_cap=zeros<mat>(2*n+4,2);
	mat B_zz=zeros<mat>(2*n+3,2);
	

	B_cap.submat(0,0,n+1,0)=b_1_cap;
	B_cap.submat(n+2,1,2*n+3,1)=b_2_cap;  //flag
	B_zz=(trans(R))*(trans(phi))*B_cap;

	//Calculating control torque:

	//Position of payload;

	double rho,rho_dot;

	rho=as_scalar((J_theta_1*theta_1)+((trans(J_e_1_trans))*q_e_1));
	
	//Options: 
	//rho=as_scalar((cst.c_1*(J_theta_1*theta_1) + cst.c_1*((trans(J_e_1_trans))*q_e_1))+((cst.c_2*(J_theta_2*theta_2)+(cst.c_2*(trans(J_e_2_trans))*q_e_2))));
    
    rho_dot=as_scalar((J_theta_1*theta_1_dot)+((trans(J_e_1_trans))*q_e_1_dot));

	
	vec tau_pt(2);   // constant torque corrosponding to pre-tension in the cables to avoid cable slack
    tau_pt(0)=cst.tau_pt_1;     // 1 N-m corrosponds to 25 N tension in the cables for radius of winch 0.04m.  
    tau_pt(1)=cst.tau_pt_2;     // Directions can be verified from the figure. 


    vec tau(2);
    vec tau_out(4);


	tau_out = tau_c_theta(t,rho,rho_dot,theta_1,theta_2,theta_1_dot,theta_2_dot,a_cap,cst.post_processing_switch_off,cst);

	tau = tau_pt + tau_out.subvec(0,1);


	//(d) Equations of Motion
	vec z_ddot(2*n+3);
	
	z_ddot = solve(M_zz,(f_non_zz + (B_zz*tau) - (((K_zz)*z)+((D_zz)*z_dot))));
	
	                                                      
	//Options:
	//nu_dot = pinv(M,tolerance)*(f_non + f_ext - K*(x.subvec(0, n + 1)));
	//nu_dot = inv(M)*(f_non + f_ext - K*(x.subvec(0,n+1)));
	//nu_dot = solve(M,(f_non + f_ext - K*(x.subvec(0,n+1))));
	
	//(e) Reallocationg variable 
	

	x_dot.subvec(2*n+4,4*n+6)=z_ddot;

	//Propogating Adaptive controller feed forward estimator:

	x_dot.subvec(4*n+7,4*n+8)=tau_out.subvec(2,3);       //a_cap_dot = tau_out.subvec(2,3); 
                                                       
                                                         //cout<<"return size"<<x_dot.size()<<endl;//test
    	                                                 //x_dot.print("x_dot =");//test

	return x_dot;

}

//14. Test ODEs function
vec ODES01(vec x, double t, constants cst)
{
	x = 0.5*x;
	return x;
}

