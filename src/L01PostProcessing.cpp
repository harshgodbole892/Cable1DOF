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

#include "L01PostProcessing.hpp"

//16. Post Processing

void post_processing(double h, double t_start, double t_end, mat x_out, int calc_lagrange_multiplier , constants cst)
{
	//Initializing variables

	long long int a = (t_end - t_start) / h; //Need better counting Subroutine 
	
	int n=cst.n;

	vec t_out(a);
	
	double theta_1,theta_1_dot, m_i_1,m_i_2, theta_2, theta_2_dot; 

	mat M = zeros<mat>(2*n + 4, 2*n + 4);
	mat M_zz = zeros<mat>(2*n + 3, 2*n + 3);

	mat K = zeros<mat>(2*n + 4, 2*n + 4);
	mat K_zz = zeros<mat>(2*n + 3, 2*n + 3);

    
	vec z(2*n+3);
    vec z_dot(2*n+3);
    vec H(a);
    
    vec rho(a);
    vec rho_2(a);
    vec rho_dot(a);
    vec rho_dot_2(a);

    vec rho_des(a);
    vec rho_des_dot(a);
    
    vec rho_mu(a);
    vec rho_mu_dot(a);

    vec M_rho_e(a);
    vec M_rho_e_error(a);
    vec M_rho_act(a);
    
    vec G_rho_e(a);
    vec G_rho_e_error(a);
    vec G_rho_act(a);

    vec t(a);

    vec tau_1(a);
    vec tau_2(a);
    vec tau_c(a);

    vec theta_1_out(a);
    vec theta_2_out(a);

    vec T_k1(a);
    vec T_k2(a);

    vec W_a_cap=zeros<vec>(a);
    vec W_a_tilde=zeros<vec>(a);
    double rho_r_pp=0;      
    



    // required if (calc_lagrange_multiplier==1)	
    

    vec f_non=zeros<vec>(2*n+4);
	vec f_non_eta=zeros<vec>(2*n+4);
	
	vec b_1_cap=zeros<vec>(n+2);
	vec b_2_cap=zeros<vec>(n+2);
	b_1_cap(0)=1;
	b_2_cap(0)=1;
	
	mat B_cap=zeros<mat>(2*n+4,2);
	mat B_eta=zeros<mat>(2*n+4,2);
	

	B_cap.submat(0,0,n+1,0)=b_1_cap;
	B_cap.submat(n+2,1,2*n+3,1)=b_2_cap;  //flag
	

	double rho_pp,rho_dot_pp, s_mu_pp;  //  Temporary variables for post processing

    vec lambda_1(a);
    vec lambda_1_temp(2*n+4);

    vec eta(2*n+4);
    vec eta_dot(2*n+4);
    vec eta_ddot(2*n+4);
    vec eta_dot_minus = zeros<vec>(2*n+4);


    vec q_1=zeros<vec>(n+2);
	vec q_2=zeros<vec>(n+2);
	vec nu_1=zeros<vec>(n+2);
	vec nu_2=zeros<vec>(n+2);
    

    vec a_cap_pp=zeros<vec>(2);

    

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


	// External torque: 
	B_eta=((trans(phi))*B_cap);
    


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

	// Above Matrices being ocnstant through post processing are defined before the loop


	//Post processing:
	    	
	for(int i=0; i<a; i++)
	
	{
		//Generating time step
		if(i==0)
			t(i)= t_start;
		else
			t(i)=t(i-1)+h;

		//Extracting required variables for post processing from x_out

		theta_1=x_out(1,i);
		theta_2=x_out(0,i);
		a_cap_pp = x_out.submat(23,i,24,i);

		m_i_1=cst.rho*cst.A*(cst.L-cst.R*theta_1)/n;      
        m_i_2=cst.rho*cst.A*(cst.L-cst.R*theta_2)/n;

        

        //(a) Mass Matrix
       
	    M.submat(0,0,n+1,n+1) = mass_matrix(n, m_i_1, theta_1, cst);
	    M.submat(n+2,n+2,2*n+3,2*n+3) = mass_matrix(n, m_i_2, theta_2, cst);
	
	    M_zz=(trans(R))*(trans(phi))*M*phi*R;
    
        //(b) Stiffness Matrix

	    K.submat(0,0,n+1,n+1) = K_matrix(n, theta_1, cst);
	    K.submat(n+2,n+2,2*n+3,2*n+3) = K_matrix(n, theta_2, cst);
    
	    K_zz=(trans(R))*(trans(phi))*K*phi*R;

	    //(c) Hameltonian

	    z=x_out.submat(1,i,2*n+3,i);
	    z_dot=x_out.submat(2*n+4,i,4*n+6,i);
	   
	    H(i)=as_scalar((0.5*((trans(z_dot))*M_zz*z_dot)+(0.5*((trans(z))*K_zz*z))));

	    


	    //Lagrange multiplier calculation starts here:

	    if(calc_lagrange_multiplier==1)
	    {

	    theta_1_dot = x_out(2*n + 4,i);//flag

	    q_1(0)=theta_1;
	    q_1.subvec(1,n+1)=x_out.submat(2,i,n+2,i);
	    q_2(0)=theta_2;
	    q_2.subvec(1,n+1)=x_out.submat(n+3,i,(2*n+3),i);

	    nu_1(0)=theta_1_dot;
	    nu_1.subvec(1,n+1)=x_out.submat((2*n+5),i,(3*n+5),i);
	    nu_2(0)=theta_2_dot;
	    nu_2.subvec(1,n+1)=x_out.submat((3*n+6),i,(4*n+6),i);


	    theta_2_dot=as_scalar((((J_theta_1/J_theta_2)*theta_1_dot)+ (((trans(J_e_1_trans)*(nu_1.subvec(1,n+1))/J_theta_2))-(((trans(J_e_2_trans)*(nu_2.subvec(1,n+1))/J_theta_2))))));	    	    

	    eta(0)=theta_1;
	    eta(1)=theta_2;
	    eta.subvec(2,2*n+3)=x_out.submat(2,i,2*n+3,i);
	    
	    eta_dot(0)=theta_1_dot;
	    eta_dot(1)=theta_2_dot;
	    eta_dot.subvec(2,2*n+3)=x_out.submat(2*n+5,i,4*n+6,i);

	    eta_ddot=(1/h)*(eta_dot-eta_dot_minus);

	    

	    eta_dot_minus=eta_dot;

	    

	    // Non Linear Components:
   
    
    	f_non.subvec(0,n+1) = -(((M_dot(n, theta_1_dot, cst))*nu_1) - (Partial_M(n, nu_1, cst)) + (Partial_k(n, theta_1, q_1, cst)));
	    f_non.subvec(n+2,2*n+3)= -(((M_dot(n, theta_2_dot, cst))*nu_2) - (Partial_M(n, nu_2, cst)) + (Partial_k(n, theta_2, q_2, cst)));

	    
  
                                                         //f_non.print("f_non");//test
                                                         //R.print("R:");
                                                         //phi.print("phi:");


        f_non_eta=(((trans(phi)))*f_non);
                                                         //f_non_zz.print("f_non_zz = ");//test

        


	    // External forces and Control Laws
	

	     rho_pp=as_scalar((J_theta_1*theta_1)+((trans(J_e_1_trans))*q_1.subvec(1,n+1)));
         rho_dot_pp=as_scalar((J_theta_1*theta_1_dot)+((trans(J_e_1_trans))*nu_1.subvec(1,n+1)));

	     
         vec tau_pt(2);    // constant torque corrosponding to pre-tension in the cables to avoid cable slack
         tau_pt(0)=-1;     //1 N-m corrosponds to 25 N tension in hte cables for radius of winch 0.04m  
         tau_pt(1)=1;

        
	     vec tau_out(4);

         vec tau(2);
	     tau_out = tau_c_theta(t(i),rho_pp,rho_dot_pp,theta_1,theta_2,theta_1_dot,theta_2_dot,a_cap_pp, cst.post_processing_switch_on, cst);

	     tau = tau_pt + tau_out.subvec(0,1);

	     tau_1(i)= tau(0);
         tau_2(i)= tau(1);
         tau_c(i)= ((tau_out(0)/(cst.c_1*J_theta_1)));

         theta_1_out(i)=theta_1;
         theta_2_out(i)=theta_2;

         W_a_cap (i)   = as_scalar(tau_out(2));
         rho_r_pp      = as_scalar(tau_out(3));
            
         //Tensions in the Cable
         T_k1(i)=(cst.E*cst.A/(cst.L- cst.R*theta_1))*sum(q_1.subvec(1,n+1));
         T_k2(i)=(cst.E*cst.A/(cst.L- cst.R*theta_2))*sum(q_2.subvec(1,n+1));


	    lambda_1_temp=(((trans(phi))*M*phi*eta_ddot)+((trans(phi))*K*phi*eta)+((trans(phi))*(damp(2*n+2,cst))*phi*eta_dot)-(f_non_eta)-(B_eta*tau));
	    lambda_1(i)=(1/J_theta_1)*lambda_1_temp(0);

	    }


	  //Calculating payload position rho;

     	rho(i)=as_scalar((J_theta_1*theta_1)+((trans(J_e_1_trans))*x_out.submat(2,i,n+2,i)));
     	rho_2(i)=as_scalar((J_theta_2*theta_2)+((trans(J_e_2_trans))*x_out.submat((n+3),i,(2*n+3),i)));

     	rho_dot(i)=as_scalar((J_theta_1*theta_1_dot)+((trans(J_e_1_trans))*nu_1.subvec(1,n+1)));
     	rho_dot_2(i)=as_scalar((J_theta_2*theta_2_dot)+((trans(J_e_2_trans))*nu_2.subvec(1,n+1)));

     	
     	rho_des(i)=as_scalar(rho_d(t(i),cst));
     	rho_des_dot(i)=as_scalar(rho_d_dot(t(i),cst));

     	rho_mu(i) = ((cst.mu)*rho(i) + ((1-cst.mu) * ( (cst.c_1*J_theta_1*theta_1) + (cst.c_2*J_theta_2*theta_2))));
        rho_mu_dot(i) = ((cst.mu)*rho_dot(i) + ((1-cst.mu) * ( (cst.c_1*J_theta_1*theta_1_dot) + (cst.c_2*J_theta_2*theta_2_dot))));
     	
     	s_mu_pp=as_scalar((rho_mu_dot(i)-rho_des_dot(i)) + (cst.lambda_k)*(rho_mu(i)-rho_des(i)));
     
     // Calculating Estimated Quantities:
     
     M_rho_e(i) = a_cap_pp(0);
     M_rho_act(i) = (((cst.J_w/(cst.R*cst.R))+((cst.rho*cst.A*cst.L)) + (cst.m_p))*2);

     M_rho_e_error(i) = M_rho_act(i)-M_rho_e(i);
     
        if (cst.regressor_switch == 1)
        {
            G_rho_e(i) = a_cap_pp(1);
        }
        if (cst.regressor_switch == 2)
        {
            G_rho_e(i) = a_cap_pp(1) * (1 / (((cst.L-cst.R*theta_1)) * ((cst.L-cst.R*theta_1))));
        }
     
        
     G_rho_act(i) = as_scalar(0.5*(cst.n+1)*cst.E*cst.A*cst.R*((((trans(q_1.subvec(1,n+1)))*(q_1.subvec(1,n+1))) / (J_theta_1*(cst.L-cst.R*theta_1)*(cst.L-cst.R*theta_1))) + (((trans(q_2.subvec(1,n+1)))*(q_2.subvec(1,n+1))) / (J_theta_2*(cst.L-cst.R*theta_2)*(cst.L-cst.R*theta_2)))));

     G_rho_e_error(i) = G_rho_act(i) - G_rho_e(i);

     W_a_tilde(i)  = M_rho_e_error(i) * rho_r_pp + G_rho_e_error(i)-cst.K_d*s_mu_pp;

	}



	//Saving relevant data:

	x_out.save(cst.GenDir + "x_out.txt", raw_ascii);
	H.save(cst.GenDir + "H.txt", raw_ascii);
    
    
    if(cst.test_switch == 1)
    {
        if(cst.adaptive_switch == 1)
        {
            rho.save(cst.GenDir + "rho_3.txt", raw_ascii);
            //Tensions in the Cable
            T_k1.save(cst.GenDir + "T_k1_3.txt" , raw_ascii);
            T_k2.save(cst.GenDir + "T_k2_3.txt" , raw_ascii);
            
            tau_1.save(cst.GenDir + "tau_1_3.txt" , raw_ascii);
            tau_2.save(cst.GenDir + "tau_2_3.txt" , raw_ascii);
            
            
            M_rho_e.save(cst.GenDir + "M_rho_e_3.txt", raw_ascii);
            M_rho_act.save(cst.GenDir + "M_rho_act_3.txt", raw_ascii);
            
            G_rho_e.save(cst.GenDir + "G_rho_e_3.txt", raw_ascii);
            G_rho_act.save(cst.GenDir + "G_rho_act_3.txt", raw_ascii);
        }
        else
        {
            rho.save(cst.GenDir + "rho_4.txt", raw_ascii);
            //Tensions in the Cable
            T_k1.save(cst.GenDir + "T_k1_4.txt" , raw_ascii);
            T_k2.save(cst.GenDir + "T_k2_4.txt" , raw_ascii);
            
            tau_1.save(cst.GenDir + "tau_1_4.txt" , raw_ascii);
            tau_2.save(cst.GenDir + "tau_2_4.txt" , raw_ascii);
        }
    }
    else
    {
        if(cst.adaptive_switch == 1)
        {
            rho.save(cst.GenDir + "rho.txt", raw_ascii);
            //Tensions in the Cable
            T_k1.save(cst.GenDir + "T_k1_1.txt" , raw_ascii);
            T_k2.save(cst.GenDir + "T_k2_1.txt" , raw_ascii);
            
            tau_1.save(cst.GenDir + "tau_1_1.txt" , raw_ascii);
            tau_2.save(cst.GenDir + "tau_2_1.txt" , raw_ascii);
            
            M_rho_e.save(cst.GenDir + "M_rho_e_1.txt", raw_ascii);
            M_rho_act.save(cst.GenDir + "M_rho_act_1.txt", raw_ascii);
            
            G_rho_e.save(cst.GenDir + "G_rho_e_1.txt", raw_ascii);
            G_rho_act.save(cst.GenDir + "G_rho_act_1.txt", raw_ascii);
        }
        else
        {
            rho.save(cst.GenDir + "rho_2.txt", raw_ascii);
            //Tensions in the Cable
            T_k1.save(cst.GenDir + "T_k1_2.txt" , raw_ascii);
            T_k2.save(cst.GenDir + "T_k2_2.txt" , raw_ascii);
            
            tau_1.save(cst.GenDir + "tau_1_2.txt" , raw_ascii);
            tau_2.save(cst.GenDir + "tau_2_2.txt" , raw_ascii);
        }
    }
    
    if (cst.regressor_switch == 1)
    {
        rho.save(cst.GenDir + "rho_5.txt", raw_ascii);
    }
    if (cst.regressor_switch == 2)
    {
        rho.save(cst.GenDir + "rho_6.txt", raw_ascii);
    }
    if (cst.regressor_switch == 3)
    {
        rho.save(cst.GenDir + "rho_7.txt", raw_ascii);
    }
    
    //rho_2.save(cst.GenDir + "rho_2.txt", raw_ascii);
    
    rho_dot.save(cst.GenDir + "rho_dot.txt", raw_ascii);
    rho_dot_2.save(cst.GenDir + "rho_dot_2.txt", raw_ascii);

    rho_des.save(cst.GenDir + "rho_des.txt", raw_ascii);
    rho_des_dot.save(cst.GenDir + "rho_des_dot.txt", raw_ascii);
    
    rho_mu.save(cst.GenDir + "rho_mu.txt", raw_ascii);
    rho_mu_dot.save(cst.GenDir + "rho_mu_dot.txt", raw_ascii);

    M_rho_e.save(cst.GenDir + "M_rho_e.txt", raw_ascii);
    M_rho_act.save(cst.GenDir + "M_rho_act.txt", raw_ascii);
    M_rho_e_error.save(cst.GenDir + "M_rho_e_error.txt", raw_ascii);

    G_rho_e.save(cst.GenDir + "G_rho_e.txt", raw_ascii);
    G_rho_act.save(cst.GenDir + "G_rho_act.txt", raw_ascii);
    G_rho_e_error.save(cst.GenDir + "G_rho_e_error.txt", raw_ascii);

    
    t.save(cst.GenDir + "t.txt" , raw_ascii);

    tau_1.save(cst.GenDir + "tau_1.txt" , raw_ascii);
    tau_2.save(cst.GenDir + "tau_2.txt" , raw_ascii);
    tau_c.save(cst.GenDir + "tau_c.txt" , raw_ascii);

    theta_1_out.save(cst.GenDir + "theta_1_out.txt" , raw_ascii);
    theta_2_out.save(cst.GenDir + "theta_2_out.txt" , raw_ascii);
   
    //Tensions in the Cable
    T_k1.save(cst.GenDir + "T_k1.txt" , raw_ascii);
    T_k2.save(cst.GenDir + "T_k2.txt" , raw_ascii);

    //Errors in nonlinear forces:
    W_a_tilde.save(cst.GenDir + "W_a_tilde.txt" , raw_ascii);
    W_a_cap.save(cst.GenDir + "W_a_cap.txt" , raw_ascii);

    if(calc_lagrange_multiplier==1)
	{    
	  lambda_1.save(cst.GenDir + "lambda_1.txt", raw_ascii);
	}    


    cout << "C++ : Text saving process complete" << endl;
}
