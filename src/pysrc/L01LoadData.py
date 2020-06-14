#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun  7 19:56:14 2020

@author: Harsh Godbole
"""
import os
import numpy             as np
import matplotlib.pyplot as plt 
import config.L01Env     as Env

# Plot script here: 
class L01DataArray: 
    
    def namestr(obj, namespace):
        return [name for name in namespace if namespace[name] is obj]

    def Plot2D(self,var1,var2,name1,name2):
        plt.figure()
        plt.plot(var1,var2)
        plt.xlabel(name1)
        plt.ylabel(name2)
        plt.grid(True)
        plt.show()
        
    # Constructor loads all data from the generated folder
    def __init__(self):
    
        self.G_rho_act      =  np.loadtxt(Env.GenDir + "G_rho_act.txt")
        self.G_rho_act_1    =  np.loadtxt(Env.GenDir + "G_rho_act_1.txt")
        self.G_rho_e        =  np.loadtxt(Env.GenDir + "G_rho_e.txt")
        self.G_rho_e_1      =  np.loadtxt(Env.GenDir + "G_rho_e_1.txt")
        self.G_rho_e_error  =  np.loadtxt(Env.GenDir + "G_rho_e_error.txt")
        self.H              =  np.loadtxt(Env.GenDir + "H.txt")
        self.M_rho_act      =  np.loadtxt(Env.GenDir + "M_rho_act.txt")
        self.M_rho_act_1    =  np.loadtxt(Env.GenDir + "M_rho_act_1.txt")
        self.M_rho_e        =  np.loadtxt(Env.GenDir + "M_rho_e.txt")
        self.M_rho_e_1      =  np.loadtxt(Env.GenDir + "M_rho_e_1.txt")
        self.M_rho_e_error  =  np.loadtxt(Env.GenDir + "M_rho_e_error.txt")
        self.T_k1           =  np.loadtxt(Env.GenDir + "T_k1.txt")
        self.T_k1_1         =  np.loadtxt(Env.GenDir + "T_k1_1.txt")
        self.T_k2           =  np.loadtxt(Env.GenDir + "T_k2.txt")
        self.T_k2_1         =  np.loadtxt(Env.GenDir + "T_k2_1.txt")
        self.W_a_cap        =  np.loadtxt(Env.GenDir + "W_a_cap.txt")
        self.W_a_tilde      =  np.loadtxt(Env.GenDir + "W_a_tilde.txt")
        self.lambda_1       =  np.loadtxt(Env.GenDir + "lambda_1.txt")
        self.rho            =  np.loadtxt(Env.GenDir + "rho.txt")
        self.rho_5          =  np.loadtxt(Env.GenDir + "rho_5.txt")
        self.rho_des        =  np.loadtxt(Env.GenDir + "rho_des.txt")
        self.rho_des_dot    =  np.loadtxt(Env.GenDir + "rho_des_dot.txt")
        self.rho_dot        =  np.loadtxt(Env.GenDir + "rho_dot.txt")
        self.rho_dot_2      =  np.loadtxt(Env.GenDir + "rho_dot_2.txt")
        self.rho_mu         =  np.loadtxt(Env.GenDir + "rho_mu.txt")
        self.rho_mu_dot     =  np.loadtxt(Env.GenDir + "rho_mu_dot.txt")
        self.t              =  np.loadtxt(Env.GenDir + "t.txt")
        self.tau_1          =  np.loadtxt(Env.GenDir + "tau_1.txt")
        self.tau_1_1        =  np.loadtxt(Env.GenDir + "tau_1_1.txt")
        self.tau_2          =  np.loadtxt(Env.GenDir + "tau_2.txt")
        self.tau_2_1        =  np.loadtxt(Env.GenDir + "tau_2_1.txt")
        self.tau_c          =  np.loadtxt(Env.GenDir + "tau_c.txt")
        self.theta_1_out    =  np.loadtxt(Env.GenDir + "theta_1_out.txt")
        self.theta_2_out    =  np.loadtxt(Env.GenDir + "theta_2_out.txt")
        self.x_out          =  np.loadtxt(Env.GenDir + "x_out.txt")
    
    # Generic function to plot all time histories of variables: 
    def PlotAll(self):
    
        self.Plot2D(self.t, self.G_rho_act     , "time",   "G_rho_act")
        self.Plot2D(self.t, self.G_rho_act_1   , "time",   "G_rho_act_1")
        self.Plot2D(self.t, self.G_rho_e       , "time",   "G_rho_e")
        self.Plot2D(self.t, self.G_rho_e_1     , "time",   "G_rho_e_1")
        self.Plot2D(self.t, self.G_rho_e_error , "time",   "G_rho_e_error")
        self.Plot2D(self.t, self.H             , "time",   "H")
        self.Plot2D(self.t, self.M_rho_act     , "time",   "M_rho_act")
        self.Plot2D(self.t, self.M_rho_act_1   , "time",   "M_rho_act_1")
        self.Plot2D(self.t, self.M_rho_e       , "time",   "M_rho_e")
        self.Plot2D(self.t, self.M_rho_e_1     , "time",   "M_rho_e_1")
        self.Plot2D(self.t, self.M_rho_e_error , "time",   "M_rho_e_error")
        self.Plot2D(self.t, self.T_k1          , "time",   "T_k1")
        self.Plot2D(self.t, self.T_k1_1        , "time",   "T_k1_1")
        self.Plot2D(self.t, self.T_k2          , "time",   "T_k2")
        self.Plot2D(self.t, self.T_k2_1        , "time",   "T_k2_1")
        self.Plot2D(self.t, self.W_a_cap       , "time",   "W_a_cap")
        self.Plot2D(self.t, self.W_a_tilde     , "time",   "W_a_tilde")
        self.Plot2D(self.t, self.lambda_1      , "time",   "lambda_1")
        self.Plot2D(self.t, self.rho           , "time",   "rho")
        self.Plot2D(self.t, self.rho_5         , "time",   "rho_5")
        self.Plot2D(self.t, self.rho_des       , "time",   "rho_des")
        self.Plot2D(self.t, self.rho_des_dot   , "time",   "rho_des_dot")
        self.Plot2D(self.t, self.rho_dot       , "time",   "rho_dot")
        self.Plot2D(self.t, self.rho_dot_2     , "time",   "rho_dot_2")
        self.Plot2D(self.t, self.rho_mu        , "time",   "rho_mu")
        self.Plot2D(self.t, self.rho_mu_dot    , "time",   "rho_mu_dot")
        self.Plot2D(self.t, self.t             , "time",   "t")
        self.Plot2D(self.t, self.tau_1         , "time",   "tau_1")
        self.Plot2D(self.t, self.tau_1_1       , "time",   "tau_1_1")
        self.Plot2D(self.t, self.tau_2         , "time",   "tau_2")
        self.Plot2D(self.t, self.tau_2_1       , "time",   "tau_2_1")
        self.Plot2D(self.t, self.tau_c         , "time",   "tau_c")
        self.Plot2D(self.t, self.theta_1_out   , "time",   "theta_1_out")
        self.Plot2D(self.t, self.theta_2_out   , "time",   "theta_2_out")
        #self.Plot2D(self.t, self.x_out         , "time",   "x_out")
        
    #Specific plots: 
    def PlotCustom(self):
        
        # rho/rho_des
        plt.figure(1)
        plt.grid(True)
        plt.plot(self.t, self.rho_des, color='k')
        plt.plot(self.t, self.rho, color='g')
        plt.ylabel("rho, rho_des")
        plt.xlabel("time")
        
         # M_rho_act/M_rho_e
        plt.figure(2)
        plt.grid(True)
        plt.plot(self.t, self.M_rho_act, color='k')
        plt.plot(self.t, self.M_rho_e,   color='g')
        plt.ylabel("M, M_est")
        plt.xlabel("time")
        
        # Show all plots 
        plt.show()
        
        


#plt.figure()
#plt.plot(rho,H)
