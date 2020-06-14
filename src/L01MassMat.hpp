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
#include "L01Constants.hpp"

using namespace std;
using namespace arma;

//Mass Matrix Functions Prototype:
#ifndef L01MASSMAT_H
#define L01MASSMAT_H

//2. Common mass matrix function for equaly divided variable mass cable segments.
mat M_i_common(int i, int n, constants cst);

//3. Mass Matrix function
mat mass_matrix(int n, double m_i, double theta_1, constants cst);

//4. Rate of change of Mass Matrix
mat M_dot(int n, double theta_1_dot, constants cst);

//5. Mass matrix dependency on state variables
vec Partial_M(int n, vec nu, constants cst);

#endif
