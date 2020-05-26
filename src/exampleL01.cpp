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
#include "L01ODE.hpp"
#include "RK4.hpp"
#include "L01PostProcessing.hpp"


using namespace std;
using namespace arma;


//Tempelate by Harsh Godbole. Reference credits Dr. James Richard Forbes

//17. Main Function
int main(int argc, char** argv) {
    
    /*
     ------------------------------
     Step 01 : Initialize ODE
     ------------------------------
     */
    
    // Initializing constants parameters class
    constants cst;
    
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

	
	//x_out.print("output:");
	//getch();
    
    
    // **** Plot-scripts in matlab:
    // **** 1. system() Directly execute test_script.m and plot script commands
    
    //system("/Applications/MATLAB_R2016a.app/bin/matlab -nodesktop -r \"run test_script.m\"");
	

	return 0;
}

//18. Test function extension to main
/*

   inside ODE for F_non

    //options :
	//vec f_non_1=zeros<vec>(n + 2);
	//vec f_non_2=zeros<vec>(n + 2);
	
	//vec f_non_1_1=zeros<vec>(n + 2);
	//vec f_non_2_1=zeros<vec>(n + 2);
	
	//vec f_non_1_2=zeros<vec>(n + 2);
	//vec f_non_2_2=zeros<vec>(n + 2);
	
	//vec f_non_1_3=zeros<vec>(n + 2);
	//vec f_non_2_3=zeros<vec>(n + 2);
	
	
	
	//f_non_1_1=-((M_dot(n, theta_1_dot, cst))*nu_1);
	//f_non_2_1=-((M_dot(n, theta_2_dot, cst))*nu_2);
	
	//f_non_1_2=- (Partial_M(n, nu_1, cst));
	//f_non_2_2=- (Partial_M(n, nu_2, cst));
	
	//f_non_1_3=(Partial_k(n, theta_1, q_1, cst));
	//f_non_2_3=(Partial_k(n, theta_2, q_2, cst));

	//f_non_1_1.print("f_non_1_1");
	//f_non_2_1.print("f_non_2_1");
	//f_non_1_2.print("f_non_1_2");
	//f_non_2_2.print("f_non_2_2");
	//f_non_1_3.print("f_non_1_3");
	//f_non_2_3.print("f_non_2_3");
	
	//over



	vec x_IC(12);
	x_IC(0) = 1;
	x_IC(1) = 2;
	x_IC(2) = 3;
	x_IC(3) = 4;
	x_IC(4) = 5;
	x_IC(5) = 6;

	x_IC(6) = 7;
	x_IC(7) = 8;
	x_IC(8) = 9;
	x_IC(9) = 10;
	x_IC(10) = 11;
	x_IC(11) = 12;

	int al, be;
	al = 3;
	be = 4;
	double ga = 1;
	mat x_out(6,6);
	x_out = ODE(x_IC,10,cst);
	x_out.print("Output:");
	
	getch();

	return 0;
	
	}
*/

// Armadillo Code examples for Referance
/* int
main(int argc, char** argv)
  {
  cout << "Armadillo version: " << arma_version::as_string() << endl;
  
  mat A(2,3);  // directly specify the matrix size (elements are uninitialised)
  
  cout << "A.n_rows: " << A.n_rows << endl;  // .n_rows and .n_cols are read only
  cout << "A.n_cols: " << A.n_cols << endl;
  
  A(1,2) = 456.0;  // directly access an element (indexing starts at 0)
  A.print("A:");
  
  A = 5.0;         // scalars are treated as a 1x1 matrix
  A.print("A:");
  
  A.set_size(4,5); // change the size (data is not preserved)
  
  A.fill(5.0);     // set all elements to a particular value
  A.print("A:");
  
  // endr indicates "end of row"
  A << 0.165300 << 0.454037 << 0.995795 << 0.124098 << 0.047084 << endr
    << 0.688782 << 0.036549 << 0.552848 << 0.937664 << 0.866401 << endr
    << 0.348740 << 0.479388 << 0.506228 << 0.145673 << 0.491547 << endr
    << 0.148678 << 0.682258 << 0.571154 << 0.874724 << 0.444632 << endr
    << 0.245726 << 0.595218 << 0.409327 << 0.367827 << 0.385736 << endr;
  
  A.print("A:");
  
  // determinant
  cout << "det(A): " << det(A) << endl;
  
  // inverse
  cout << "inv(A): " << endl << inv(A) << endl;
  
  // save matrix as a text file
  A.save("A.txt", raw_ascii);
  
  // load from file
  mat B;
  B.load("A.txt");
  
  // submatrices
  cout << "B( span(0,2), span(3,4) ):" << endl << B( span(0,2), span(3,4) ) << endl;
  
  cout << "B( 0,3, size(3,2) ):" << endl << B( 0,3, size(3,2) ) << endl;
  
  cout << "B.row(0): " << endl << B.row(0) << endl;
  
  cout << "B.col(1): " << endl << B.col(1) << endl;
  
  // transpose
  cout << "B.t(): " << endl << B.t() << endl;
  
  // maximum from each column (traverse along rows)
  cout << "max(B): " << endl << max(B) << endl;
  
  // maximum from each row (traverse along columns)
  cout << "max(B,1): " << endl << max(B,1) << endl;
  
  // maximum value in B
  cout << "max(max(B)) = " << max(max(B)) << endl;
  
  // sum of each column (traverse along rows)
  cout << "sum(B): " << endl << sum(B) << endl;
  
  // sum of each row (traverse along columns)
  cout << "sum(B,1) =" << endl << sum(B,1) << endl;
  
  // sum of all elements
  cout << "accu(B): " << accu(B) << endl;
  
  // trace = sum along diagonal
  cout << "trace(B): " << trace(B) << endl;
  
  // generate the identity matrix
  mat C = eye<mat>(4,4);
  
  // random matrix with values uniformly distributed in the [0,1] interval
  mat D = randu<mat>(4,4);
  D.print("D:");
  
  // row vectors are treated like a matrix with one row
  rowvec r;
  r << 0.59119 << 0.77321 << 0.60275 << 0.35887 << 0.51683;
  r.print("r:");
  
  // column vectors are treated like a matrix with one column
  vec q;
  q << 0.14333 << 0.59478 << 0.14481 << 0.58558 << 0.60809;
  q.print("q:");
  
  // convert matrix to vector; data in matrices is stored column-by-column
  vec v = vectorise(A);
  v.print("v:");
  
  // dot or inner product
  cout << "as_scalar(r*q): " << as_scalar(r*q) << endl;
  
  // outer product
  cout << "q*r: " << endl << q*r << endl;
  
  // multiply-and-accumulate operation (no temporary matrices are created)
  cout << "accu(A % B) = " << accu(A % B) << endl;
  
  // example of a compound operation
  B += 2.0 * A.t();
  B.print("B:");
  
  // imat specifies an integer matrix
  imat AA;
  imat BB;
  
  AA << 1 << 2 << 3 << endr << 4 << 5 << 6 << endr << 7 << 8 << 9;
  BB << 3 << 2 << 1 << endr << 6 << 5 << 4 << endr << 9 << 8 << 7;
  
  // comparison of matrices (element-wise); output of a relational operator is a umat
  umat ZZ = (AA >= BB);
  ZZ.print("ZZ:");
  
  // cubes ("3D matrices")
  cube Q( B.n_rows, B.n_cols, 2 );
  
  Q.slice(0) = B;
  Q.slice(1) = 2.0 * B;
  
  Q.print("Q:");
  
  // 2D field of matrices; 3D fields are also supported
  field<mat> F(4,3); 
  
  for(uword col=0; col < F.n_cols; ++col)
  for(uword row=0; row < F.n_rows; ++row)
    {
    F(row,col) = randu<mat>(2,3);  // each element in field<mat> is a matrix
    }
  
  F.print("F:");
  getch();
  
  return 0;
  }
  */
