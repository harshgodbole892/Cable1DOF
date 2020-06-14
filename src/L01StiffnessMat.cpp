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

#include "L01StiffnessMat.hpp"

//6. Stiffness matrix dependencies on state variables.
vec Partial_k(int n, double theta_1, vec q_1, constants cst)
{
	mat K_1 = zeros<mat>(n + 2, n + 2);
	for (int i = 0; i < n + 2; i++)
		for (int j = 0; j < n + 2;j++)
			if (i == j &&i != 0)
			{
				K_1(i, j) = 1;
			}
	K_1 = ((n+1)*cst.E*cst.R*cst.A / ((cst.L - cst.R*theta_1) *(cst.L - cst.R*theta_1)))*K_1;
	
	vec K = zeros<vec>(n + 2);
	
	// Partial_M_partial_q turns out to be a row matrix. Its transpose is used in Lagranges Equation.
	// Predefining the matrix M as a column matrix such that the function trans() does not need to be applied twice. 

	K.subvec(0, 0) = 0.5*trans(q_1)*K_1*q_1;
	return K;
}

//7. Stiffness matrix
mat K_matrix(int n, double theta_1, constants cst)
{
	mat K_1 = zeros<mat>(n + 2, n + 2);
	for (int i = 0; i < n + 2; i++)
		for (int j = 0; j < n + 2; j++)
			if (i == j &&i != 0)
			{
				K_1(i, j) = 1;
			}
	K_1 = ((n+1)*cst.E*cst.A / ((cst.L - cst.R*theta_1)))*K_1;

	return K_1;
}

//8. Damping Matrix
mat damp(int n, constants cst)
{
	mat D_1 = zeros<mat>(n + 2, n + 2);
	for (int i = 0; i < n + 2; i++)
		for (int j = 0; j < n + 2; j++)
			if (i == j &&i != 0)
			{
				D_1(i, j) = 1;
			}
	D_1=cst.c_i*D_1;		
	return D_1;
}
