#!/bin/bash

#This is an example script for calling DecomposeSolver.py
#The user should change the full paths
#The user will likely wish to change the name of the folder for results (2nd argument) and the
#tumor name (3rd argument)
#go the LLSolver that contains the main codes
cd '/pghbio/cure/hylei/TumorDecompose/LLSolver/'

#1st argv: parent directory that will contain the result folder to save the result
#2nd argv: specify a folder name to store the results
#3rd argv: tumor name, GBM07 or GBM33
#4th argv: tumor number, in this serial, we do 3, 6, 9
#5th argv: regularization parameter for modified NMF: will do 0.0, 0.2, 0.4, 0.6, 0.8, 1.0
#6th argv: regularization parameter for ILP tree structure using gurobi: will do 0.0, 0.2, 0.4, 0.6, 0.8, 1.0
#7th argv: regularization parameter for SCIP, not available yet, put 0.0
#8th argv: solver name: nmf or gurobi

#call to deconvolve 3 samples of GBM33 using NMF
python DecomposeSolver.py '/pghbio/cure/hylei/TumorDecompose/' 9_28 GBM33 3 0.2 0.2 0.0 nmf
