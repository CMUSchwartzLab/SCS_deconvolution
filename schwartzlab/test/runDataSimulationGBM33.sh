#!/bin/bash

#go to the directory that has the code

cd '/pghbio/cure/hylei/TumorDecompose/LLSolver'

:<<'COMMENT'
set up the argument for variable
1s argv: string, the directory for the parent directory that contains subfolders such as 
            code/
			simulation/ 
			results/ 
			data/ 
			etc..
2nd argv: string, specify a folder name to store the simualtion, 
          I prefer to use the date to name the folder
3nd argv: string, this is to choose tumor name, GBM07 or GBM33
4th argv: integer, tumor number to simulate
5th argv: list, alpha for Dirichlet Distribution, we test [10, 0.1, 0.01] is a good valid
		   parameter  for Dirichlet alpha
6th argv: integer, number of simulated data, which is also the number of experiments for later
          in the current paper, we simulate 40 replicate simulated data for each tumor sample case
7th argv: Boolean, if (True) or not (False) cap the copy numbers that are larger than 10 to 10
COMMENT

#The call is to simulate 3, 6, 9 tumor samples from the SCS data from GBM33
python DataSimulation.py '/pghbio/cure/hylei/TumorDecompose/' 9_28 GBM33 3 [10,0.1,0.01] 40 True
python DataSimulation.py '/pghbio/cure/hylei/TumorDecompose/' 9_28 GBM33 6 [10,0.1,0.01] 40 True
python DataSimulation.py '/pghbio/cure/hylei/TumorDecompose/' 9_28 GBM33 9 [10,0.1,0.01] 40 True

