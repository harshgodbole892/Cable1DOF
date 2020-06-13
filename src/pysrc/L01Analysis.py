#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun  7 19:56:14 2020

@author: Harsh Godbole
"""
import os
os.chdir("../../")


from src.pysrc.L01CppRunner import L01CppRunner as CppR
from src.pysrc.L01LoadData  import L01DataArray as DatArr

# Created an instance of CPP runner to build CPP solution: 
runner = CppR()
runner.Build()
runner.RunTest()

# Plot script here: 

op = DatArr()
#op.PlotAll()
op.PlotCustom()
