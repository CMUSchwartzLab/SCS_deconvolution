#!/bin/bash
#This is an example script to call DataSimulation.py for tumor GBM33
#go to the directory that has the code, this directory has to be changed by the user

cd '/pghbio/cure/hylei/TumorDecompose/LLSolver'

:<<'COMMENT'
this is a description of the seven arguments
1s argv: string, the full path to the parent directory that contains subfolders such as 
            LLSolver/
			simulation/ 
			results/ 
			data/ 
			etc..
2nd argv: string, specify a folder name to store the simualtion, 
          we suggest to use the date to name this folder
3nd argv: string, this is to choose tumor name, GBM07 or GBM33
4th argv: integer, tumor number to simulate
5th argv: list, alpha for Dirichlet Distribution, we suugest [10, 0.1, 0.01] as a good valid
		   parameter  set for Dirichlet alpha
6th argv: integer, number of simulated replicates, which is also the number of experiments for later
          in the current paper, we simulate 40 replicates for each tumor sample case
7th argv: Boolean, if (True) or not (False) cap the copy numbers that are larger than 10 to 10
COMMENT

#The call is to simulate 3, 6, 9 tumor samples from the SCS data from GBM33
#At least the directories will need to be changed by the user depending on the user's file naming conventions
python DataSimulation.py '/pghbio/cure/hylei/TumorDecompose/' 9_28 GBM33 3 [10,0.1,0.01] 40 True
python DataSimulation.py '/pghbio/cure/hylei/TumorDecompose/' 9_28 GBM33 6 [10,0.1,0.01] 40 True
python DataSimulation.py '/pghbio/cure/hylei/TumorDecompose/' 9_28 GBM33 9 [10,0.1,0.01] 40 True

