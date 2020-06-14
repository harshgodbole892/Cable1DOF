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

#include "L01ControlLaws.hpp"

//Function definations for Control Laws start here: 

//9. Desired trajectory
double rho_d(double t_act, constants cst)
{
	double rho_d_t;
	double rho_i;
	double rho_f;
	double t_f;
	double t_i;
	int t_temp;
	double t;

	t = t_act;
 if(t_act>cst.settle_time)
 {
  if(cst.rho_d_switch==1)
  {  
  	//For propogating in time:

  	 t_temp = t_act/2;

  	 t = t_act-(t_temp*2);

    //Traverse to rho_f
	if( t >= 0 && t < 0.1)
	{
		rho_i=0;
		rho_f=0.2*cst.L;
		
		t_i=t;   //start time of step transition
		t_f=0.1; //end time of step transition
		
	}

	//Stationary at rho_f 

	if( t >= 0.1 && t < 0.5)
	{
		rho_i=0;
		rho_f=0.2*cst.L;
		t_f=t;
		t_i=t;
	}

    //Traverse to rho_f
	if( t >= 0.5 && t < 0.6)
	{
		rho_i=0.2*cst.L;
		rho_f=-0.2*cst.L;
		
		t_i=t-0.5;
		t_f=0.1;
	}

	//Stationary at rho_f 
	if( t >= 0.6 && t < 1)
	{
		rho_i=0.2*cst.L;
		rho_f=-0.2*cst.L;
		t_i=t;
		t_f=t;
	}

	//Traverse to rho_f
	if( t >= 1 && t < 1.1)
	{
		rho_i=-0.2*cst.L;
		rho_f=0.2*cst.L;
		
		t_i=t-1.0;
		t_f=0.1;
	}

	//Stationary at rho_f 
	if( t >= 1.1 && t < 1.5)
	{
		rho_i=0;
		rho_f=0.2*cst.L;
		t_i=t;
		t_f=t;
	}

	//Traverse to rho_f
	if( t >= 1.5 && t < 1.6)
	{
		rho_i=0.2*cst.L;
		rho_f=0;
		t_i=t-1.5;
		t_f=0.1;
	}

	//Stationary at rho_f 
	if( t >= 1.6 && t <= 2)
	{
		rho_i=0.2*cst.L;
		rho_f=0;
		t_i=t;
		t_f=t;
	}


    rho_d_t= ( ( (10*(t_i/t_f)*(t_i/t_f)*(t_i/t_f)) - (15*(t_i/t_f)*(t_i/t_f)*(t_i/t_f)*(t_i/t_f)) + (6*(t_i/t_f)*(t_i/t_f)*(t_i/t_f)*(t_i/t_f)*(t_i/t_f)) )*(rho_f-rho_i) ) + rho_i;
  }     

  if(cst.rho_d_switch==2)
  { 
  	rho_d_t=cst.L*0.2;
  }

  if(cst.rho_d_switch==3)
  {
  	rho_d_t=0.5*t*t;
  }

  if(cst.rho_d_switch==4)
  {
  	rho_d_t=(((cst.L)/6)*sin(2*cst.pi*cst.omega_d*t));
  }

  
  if(cst.rho_d_switch==5)
  {
     t_temp = t_act;
     
     t = t - t_temp;

     if(t<0.5)
   
   	 rho_d_t=0.2*cst.L;
   
     else
     
     rho_d_t=-0.2*cst.L;
  }
 }
 else
 {
  rho_d_t = 0;
 }

  	//options: 
    //rho_d_t=cst.L*0.25;
    //rho_d_t=0.5*t*t;
    //rho_d_t=(((cst.L)/4)*sin(2*cst.pi*0.1*t));
    
    return rho_d_t;
}

//10. Velocity of desired trajectory
double rho_d_dot(double t_act, constants cst)
{
	double rho_d_t;
	double rho_i;
	double rho_f;
	double t_f;
    double t_i;
    int t_temp;
    double t;

    t = t_act;
 if(t_act>cst.settle_time)
 {	
  if(cst.rho_d_switch==1)
  {
  	// For propogating in time

  	 t_temp = t_act/2;

  	 t = t_act-(t_temp*2);


    //Traverse to rho_f
	if( t >= 0 && t < 0.1)
	{
		rho_i=0;
		rho_f=0.2*cst.L;
		
		t_i=t;   //start time of step transition
		t_f=0.1; //end time of step transition
		
	}

	//Stationary at rho_f 

	if( t >= 0.1 && t < 0.5)
	{
		rho_i=0;
		rho_f=0.2*cst.L;
		t_f=t;
		t_i=t;
	}

    //Traverse to rho_f
	if( t >= 0.5 && t < 0.6)
	{
		rho_i=0.2*cst.L;
		rho_f=-0.2*cst.L;
		
		t_i=t-0.5;
		t_f=0.1;
	}

	//Stationary at rho_f 
	if( t >= 0.6 && t < 1)
	{
		rho_i=0.2*cst.L;
		rho_f=-0.2*cst.L;
		t_i=t;
		t_f=t;
	}

	//Traverse to rho_f
	if( t >= 1 && t < 1.1)
	{
		rho_i=-0.2*cst.L;
		rho_f=0.2*cst.L;
		
		t_i=t-1.0;
		t_f=0.1;
	}

	//Stationary at rho_f 
	if( t >= 1.1 && t < 1.5)
	{
		rho_i=0;
		rho_f=0.2*cst.L;
		t_i=t;
		t_f=t;
	}

	//Traverse to rho_f
	if( t >= 1.5 && t < 1.6)
	{
		rho_i=0.2*cst.L;
		rho_f=0;
		t_i=t-1.5;
		t_f=0.1;
	}

	//Stationary at rho_f 
	if( t >= 1.6 && t <= 2)
	{
		rho_i=0.2*cst.L;
		rho_f=0;
		t_i=t;
		t_f=t;
	}

    rho_d_t= ( ( (10*(3/t_f)*(t_i/t_f)*(t_i/t_f)) - (15*(4/t_f)*(t_i/t_f)*(t_i/t_f)*(t_i/t_f)) + (6*(5/t_f)*(t_i/t_f)*(t_i/t_f)*(t_i/t_f)*(t_i/t_f)) )*(rho_f-rho_i) );
  } 


  if(cst.rho_d_switch==2)
  {
  	rho_d_t=0;
  }

  if(cst.rho_d_switch==3)
  {
  	rho_d_t=t;
  }

  if(cst.rho_d_switch==4)
  {
  	rho_d_t=(((cst.L)*2*cst.pi*cst.omega_d/6)*cos(2*cst.pi*cst.omega_d*t));
  }
  	

  if(cst.rho_d_switch==5)
  {
     t_temp = t_act;
     
     t = t - t_temp;

     if(t<0.5)
   
   	 rho_d_t=0;
   
     else
     
     rho_d_t=0;
  }
 }
 else
 {
	rho_d_t=0;
 } 

	//Options:
    //rho_d_t=0;
    //rho_d_t=t;
    //rho_d_t=(((cst.L)*2*cst.pi*0.1/4)*cos(2*cst.pi*0.1*t));

    return rho_d_t;   
}

//11. Acceleration of desired trajectory
double rho_d_ddot(double t_act, constants cst)
{
    double rho_d_t;
	double rho_i;
	double rho_f;
	double t_f;
	double t_i;
	int t_temp;
	double t;

	t = t_act;
 if(t_act>cst.settle_time)
 {  	
  if(cst.rho_d_switch==1)
  {
  	//For propogating in time:

  	 t_temp = t_act/2;

  	 t = t_act-(t_temp*2);


    //Traverse to rho_f
	if( t >= 0 && t < 0.1)
	{
		rho_i=0;
		rho_f=0.2*cst.L;
		
		t_i=t;   //start time of step transition
		t_f=0.1; //end time of step transition
		
	}

	//Stationary at rho_f 

	if( t >= 0.1 && t < 0.5)
	{
		rho_i=0;
		rho_f=0.2*cst.L;
		t_f=t;
		t_i=t;
	}

    //Traverse to rho_f
	if( t >= 0.5 && t < 0.6)
	{
		rho_i=0.2*cst.L;
		rho_f=-0.2*cst.L;
		
		t_i=t-0.5;
		t_f=0.1;
	}

	//Stationary at rho_f 
	if( t >= 0.6 && t < 1)
	{
		rho_i=0.2*cst.L;
		rho_f=-0.2*cst.L;
		t_i=t;
		t_f=t;
	}

	//Traverse to rho_f
	if( t >= 1 && t < 1.1)
	{
		rho_i=-0.2*cst.L;
		rho_f=0.2*cst.L;
		
		t_i=t-1.0;
		t_f=0.1;
	}

	//Stationary at rho_f 
	if( t >= 1.1 && t < 1.5)
	{
		rho_i=0;
		rho_f=0.2*cst.L;
		t_i=t;
		t_f=t;
	}

	//Traverse to rho_f
	if( t >= 1.5 && t < 1.6)
	{
		rho_i=0.2*cst.L;
		rho_f=0;
		t_i=t-1.5;
		t_f=0.1;
	}

	//Stationary at rho_f 
	if( t >= 1.6 && t <= 2)
	{
		rho_i=0.2*cst.L;
		rho_f=0;
		t_i=t;
		t_f=t;
	}


    rho_d_t= ( ( (10*(3/t_f)*(2/t_f)*(t_i/t_f)) - (15*(4/t_f)*(3/t_f)*(t_i/t_f)*(t_i/t_f)) + (6*(5/t_f)*(4/t_f)*(t_i/t_f)*(t_i/t_f)*(t_i/t_f)) )*(rho_f-rho_i) );
  } 
  
    if(cst.rho_d_switch==2)
  {
  	rho_d_t=0;
  }

  if(cst.rho_d_switch==3)
  {
  	rho_d_t=1;
  }

  if(cst.rho_d_switch==4)
  {
  	rho_d_t=-(((cst.L)*(2*cst.pi*cst.omega_d)*(2*cst.pi*cst.omega_d)/6)*sin(2*cst.pi*cst.omega_d*t));
  }

  if(cst.rho_d_switch==5)
  {
     t_temp = t_act;
     
     t = t - t_temp;

     if(t<0.5)
   
   	 rho_d_t=0;
   
     else
     
     rho_d_t=0;
  } 
 }
 else
 {
	rho_d_t=0;
 }

    //options
    //rho_d_t=0;
    //rho_d_t=1;
    //rho_d_t=-(((cst.L)*(2*cst.pi*0.1)*(2*cst.pi*0.1)/4)*sin(2*cst.pi*0.1*t));

    return rho_d_t;
}

//12. Main PD Feedback and Feed Forward control law
vec tau_c_theta(double t,double rho, double rho_dot, double theta_1, double theta_2, double theta_1_dot, double theta_2_dot,vec a_cap, int post_processing_switch, constants cst)
{
    // Adaptive control law starts here:
    
	vec tau_c_theta_out(4);

	int n;
	n=cst.n;

	
	double tau_FF_rho,tau_FB_rho;
	
   
    //Feedback Control Law:

      //Extracted variables will be: 

	  //theta_1=q(0);
	  //theta_1_dot=q_dot(0);
	  //q_e_1=q.subvec(1,n+1);
	  //q_e_1_dot=q_dot.subvec(1,n+1);

	double J_theta_1;
	double J_theta_2;

    J_theta_1=-cst.R;
    J_theta_2=cst.R;
	
	//vec J_e_1_trans=ones<vec>(n+1);
	//vec J_e_2_trans=ones<vec>(n+1);
    
    //J_e_2_trans=(-1)*J_e_2_trans;
    

    double rho_mu,rho_mu_dot;
    
    rho_mu = ((cst.mu)*rho + ((1-cst.mu) * ( (cst.c_1*J_theta_1*theta_1) + (cst.c_2*J_theta_2*theta_2))));
    rho_mu_dot = ((cst.mu)*rho_dot + ((1-cst.mu) * ( (cst.c_1*J_theta_1*theta_1_dot) + (cst.c_2*J_theta_2*theta_2_dot))));

    double rho_mu_tilda,rho_mu_tilda_dot;
    
    rho_mu_tilda= (rho_mu) - (((cst.mu)*rho_d(t,cst) + ((1-cst.mu) * ( (cst.c_1*J_theta_1*theta_1) + (cst.c_2*J_theta_2*theta_2)))));
    rho_mu_tilda_dot= ((rho_mu_dot) - (rho_d_dot(t,cst)));


    //Feedback control torque:

    double s_mu;

    s_mu = rho_mu_tilda_dot + (cst.lambda_k * rho_mu_tilda );

    tau_FB_rho = - (cst.K_d * s_mu);
    

    //Adaptive Feed Forward Control law:
		 
	 //Compensator:

    double rho_r_ddot;

    rho_r_ddot = rho_d_ddot(t,cst) - (cst.lambda_k * rho_mu_tilda_dot);

    vec W_trans=zeros<vec>(2);
    
    W_trans(0)=rho_r_ddot;
    //W_trans(1)=1;
    

    //Options:
    if (cst.regressor_switch == 1)
    {
        W_trans(1) = 1;
    }
    if (cst.regressor_switch == 2)
    {
        W_trans(1) = (1 / (((cst.L-cst.R*theta_1)) * ((cst.L-cst.R*theta_1))));
    }
    if (cst.regressor_switch == 3)
    {
        W_trans(0) = 0;
        W_trans(1) = 1;
    }

    
	tau_FF_rho = as_scalar((trans(W_trans))*a_cap);


    //Total Control torque: 
    double tau_c_rho;

    if(cst.adaptive_switch == 1)
    {
        tau_c_rho = tau_FF_rho + tau_FB_rho; //tau_FB is negative due to the nature of positive sign of theta_1 which causes displacement in the negative x direction.
    }
    else
    {
        if(cst.feedforward_switch == 1)
            tau_c_rho = tau_FB_rho + (cst.m_rr * rho_r_ddot);
        else
            tau_c_rho = tau_FB_rho;
        
    }
    tau_c_theta_out(0)=(cst.c_1)*J_theta_1*tau_c_rho;
    tau_c_theta_out(1)=(cst.c_2)*J_theta_2*tau_c_rho;

    // Propogating the Controller states:
    mat Gamma=zeros<mat>(2,2);
    
    Gamma(0,0)=cst.Gamma_1;
    Gamma(1,1)=cst.Gamma_2;
    

    if(post_processing_switch==0)
    {
    	tau_c_theta_out.subvec(2,3)= - Gamma * W_trans * s_mu;

    }

    if(post_processing_switch==1)
    {
       tau_c_theta_out(2)=  tau_FF_rho;
       tau_c_theta_out(3)=  rho_r_ddot;
    }
    

    return tau_c_theta_out;
}
