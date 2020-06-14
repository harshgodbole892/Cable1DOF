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

#include "L01MassMat.hpp"

//2. Common mass matrix function for equaly divided variable mass cable segments.
mat M_i_common(int i, int n, constants cst)
{    
	// Creates the mass matrix tempelate used M and M_dot. We can obtain the mass matrix by multiplying by m_i or m_i_dot 
	//This method will also assume that m_i is same for all of the fragments
	mat M_i=zeros<mat>(n + 2, n + 2);
	
	M_i(0, 0) = cst.R*cst.R;
	for (int j = 0; j < n + 2; j++)
		for (int k = 0; k < n + 2; k++)
		{
			if (j <= i && k <= i)
			{
				M_i(j, k) = 1;
				if (j == 0 || k == 0)
					M_i(j, k) = -M_i(j, k)*cst.R;
			}
			else
			{
				M_i(j, k) = 0;
			}
		}
	M_i(0, 0) = cst.R*cst.R;
	return M_i;
}

//3. Mass Matrix function
mat mass_matrix(int n, double m_i, double theta_1, constants cst)
{
	mat sigma_M_i = zeros<mat>(n + 2, n + 2);
	mat sigma_M_p = zeros<mat>(n + 2, n + 2);
	mat sigma_M_w = zeros<mat>(n + 2, n + 2);

	for (int i = 1; i < n + 1; i++)
	{
		sigma_M_i = sigma_M_i + M_i_common(i, n, cst);
	}
	sigma_M_i = m_i*sigma_M_i;

	sigma_M_p = cst.m_p*M_i_common(n + 1, n, cst);

	sigma_M_w(0, 0) = (cst.J_w + (cst.rho*cst.A*((cst.R) *(cst.R) *(cst.R))*theta_1));

	mat M = zeros<mat>(n + 2, n + 2);
	M = sigma_M_i + sigma_M_p + sigma_M_w;
	return M;
}

//4. Rate of change of Mass Matrix
mat M_dot(int n, double theta_1_dot, constants cst)
{
	mat sigma_M_i_dot = zeros<mat>(n + 2, n + 2);

	for (int i = 1; i < n + 1; i++)
	{
		sigma_M_i_dot = sigma_M_i_dot + M_i_common(i, n, cst);
	}
	
	sigma_M_i_dot = ((-cst.rho)*cst.A*cst.R*theta_1_dot / n)*sigma_M_i_dot;
	
	
	/*
	M_w_dot.zeros();
	M_w_dot(0, 0) = cst.rho*cst.A*(cst.R*cst.R*cst.R)*theta_1_dot;
	*/
	
	//Analytically M_dot has a zero at position (1,1)
	sigma_M_i_dot(0,0)=0;

	return sigma_M_i_dot;
}

//5. Mass matrix dependency on state variables
vec Partial_M(int n, vec nu, constants cst)
{
	mat sigma_M_i_partial= zeros<mat>(n + 2, n + 2);
	
	for (int i = 1; i < n + 1; i++)
	{
		sigma_M_i_partial = sigma_M_i_partial + M_i_common(i, n, cst);
	}

	sigma_M_i_partial = ((-cst.rho)*cst.A*cst.R / n)*sigma_M_i_partial;

	//mat M_w_partial = zeros<mat>(n + 2, n + 2);
	//M_w_partial(0, 0) = cst.rho*cst.A*(cst.R*cst.R*cst.R);

	sigma_M_i_partial(0, 0) = 0;

	vec  M = zeros<vec>(n + 2);
	
	// partial_M_partial_q turns out to be a row matrix. Its transpose is used in Lagranges Equation.
	// predefining the matrix M as a column matrix such that the function trans() does not need to be applied twice. 

	M.subvec(0,0) = 0.5*trans(nu)*(sigma_M_i_partial)*nu;
	
	//Analytically M_dot has a zero at position (1,1) // check 
	//sigma_M_i_dot(0, 0)+ M_w_partial(0,0) = 0;
	
	return M;
}
