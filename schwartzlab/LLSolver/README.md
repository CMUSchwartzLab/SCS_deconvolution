# LLSolver

This directory is part of the repository
https://github.com/CMUSchwartzLab/SCS_deconvolution.git
that contains software to solve several formulations of the problem of deconvolving bulk tumor data into subpopulations,
as described in (Lei et al., in preparation).

The manuscript describes two methods called "phylogeny-free' and "phylogeny-based". Much of the documentation below assumes
that the user has at least skimmed the manuscript and is familiar with the terminology therein.

The present structure of the repository is that, except for README.md, all files are in the subdirectories
* schwartzlab/LLSolver
* schwartzlab/data
* schwartzlab/test

The setup assumes that the user will create another subdirectory
`schwartzlab/simulation`.

The exact naming of the subdirectories can vary, but it is inherent and some of the code and
documentation that the four subdirectories {LLSolver, data, test, simulation} are parallel, at the same level.

This LLsolver subdirectory contains the main prorgams to 
1) simulate data
2) solve the deconvolution problem in two different ways
3) generate the figures for the manuscript

Programs were written initially by Haoyun Lei and Bchuan Lyu. Programs were modified by Haoyun Lei and E. Michael Gertz. Programs were tested by
E. Michael Gertz, Haoyun Lei, Bchuan Lyu, and Alejandro Schaffer.

All programs are written in python3 (not python2). Some programs assume the avilability of the Gurobi package
(http://www.gurobi.com/downloads/download-center)
to solve optimization problems. We are also testing SCIP (https://scip.zib.de/index.php#download) as an alternative to Gurobi, but all analyses in the
manuscript that used an optimization package were done with Gurobi.
Several programs assume the availability of the numpy and scipy packages.
Example shell scripts call
  python
with the assumption that on the user's system python defaults to python 3. If that is not the case, then
it is necessary to change calls to python to calls to python3.

The main programs that a user may want to try are:
* **SimulateSCS.py**      [simulate single-cell dataresembling the observed data]
* **DataSimulation.py**   [simulate many replicates of the hypothetical bulk data from the output of SimulateSCS.py or from the observed single-cell data]
* **DecomposeSolver.py**  [solve the deconvolution problem]

DecomposeSolver.py will in turn use one of 
  * NMF_solver.py
  * GurobiILP_solver.py
  * SCIP_solver.py

to solve the optimization problems in whichever of the two formulations is chosen. NMF_solver.py solves the optimization
problem in the phylogeny-free method, while GurobiILP_solver and SCIP_solver solve the optimization problem in the
phylogeny-based method.

Several utility functions are in the code file testFunction.py


-------------------
## SimulateSCS.py

The purpose of SimulateSCS.py is to simulate realistic single-cell data that is similar to the observed single-cell data.
The program SimulateSCS.py is needed because the observed data are human subjects data, which cannot be redistributed.
SimulateSCS.py uses summary statistics from he observed data, which can be found in subdirectory
schwartzlab/data

SimulateSCS.py takes four arguments:
1st: The full path to the data directory
2nd: The tumor name, which is GBM07 or GBM33
3rd: The integer mean of the Poissson distribution used in the simulation (denoted by lam, short for lambda, in the code)
4th: The depth of the binary tree of single cells that defines the clone structure

Example outputs are shown in
* schwartzlab/data/simulated_GBM07_integer_CNV.csv
* schwartzlab/data/simulated_GBM33_integer_CNV.csv

and these were obtained via the calls
* `python SimulateSCS.py 'PATH/TO/Cell/Information' GBM07 115 6`
* `python SimulateSCS.py 'PATH/TO/Cell/Information' GBM33 45 6`

The lam (3rd argument values) of 115 and 45 are recommended for GBM07 and GBM33 respectively.

If the user wishes to change the seed for the random number generator, then it is necessary to modify SimulateSCS.py
and make an assignment to variable
seed
at the top of the program.

---------------------
## DataSimulation.py 
  
This program simulates bulk tumor data with desired number of samples from the single cell sequencing data
Simulation is based on the Geometric Model described in the manuscript and uses a Dirichlet distribution

  The arguments to DataSimulation.py are as follows:

1.  ParentDirectory: specify a directory that contains the LLSolver folder
2.  DateFolder: allows for different simulations run on different dates to be stored n different subdirectories (a.k.a. folders)
3.  TumorName: pick a tumor from which you choose the single cell data; the name of the input files should be ../data/<TumorName>_Integer_CNV.csv
4.  tumor_number: choose how many tumor samples you want to simulate, in the main experiments in the manuscript, we chose 3, 6, or 9 samples
5.  alpha: the alpha parameter for the Dirichlet distribution, describe in the manuscript
6.  N: the number of simulated replicates (the manuscript uses N=40)
7.  Cap: True of False, whether or not to cap the copy numbers larger than 10 at 10 (this was set to True for all runs in the manuscript and must be set to True for the example data provided)

Example calls to DataSimulation.py can be found in the scripts
schwartzlab/test/runDataSimulationGBM07.sh
schwartzlab/test/runDataSimulationGBM33.sh
but those are for the observed data.
If one wants to use the simulated data, simulated_GBM07_integer_CNV.csv or simulated_GBM07_integer_CNV.csv, then
the third argument, TumorName, should be specified as
simulated_GBM07
or
simulated_GBM33
respectively.

To reuse the scripts runDataSimulationGBM07.sh and runDataSimulationGBM33.sh, the user must also replace the example full paths give, with other full paths that are suitable
for the user's installation of this repository.

After running DataSimulation.py the folder/subdirectory structure should be abstractly as follows:
/some/path/to/the/ParentDirectory:
                                  /data
                                  /LLSolver/DataSimulation.py
                                           /DecomposeSolver.py
                                           /....
                                           
                                  /simulation/DateFolder/GBM07/3/simulateData1
                                                                /simulateData2
                                                                ...
                                                                /simulateDataN
                                                                
                                  /test/GTest/<TestCase1>.sh
                                             /<TestCase2>.sh
                                             ...
                                             /<TestCaseN>.sh
                                       /NTest/<TestCase1>.sh
                                             /<TestCase2>.sh
                                             ...
                                             /<TestCaseN>.sh

Here, GBM07 is the tumor name selected and could be GBM33 instead). Here, 3 is number of samples.
GTest refers to the phylogeny-based method, which uses the Gurobi (hence the G) python library, while Ntest refers to the phylogeny-free method, which is also called non-negative
matrix factorization (NMF, hence the N).

`<TestCase1>.sh` through `<TestCaseN>.sh` are abstract names for the shell scripts that call DecomposeSolver.py.

DataSimulation.py contains two auxiliary functions

  * ImportSCData(ImportSCData(ParentDirectory, TumorName, IntCNV=True, Cap=False)
    #this is to import the data from single cell sequencing data
    IntCNV: if or not choose the integer Copy number
  
  * SimulateData(N=40, Cap=Cap)
  #This is to simulate the desired number of simulated data and save them in an independent folder as described above

----------------------
## DecomposeSolver.py
  This is to decompose the bulk tumor data to resolve the Copy number in loci in each fundamental cell type. 
  The code will retrieve the simulation data if you set up the input arguments correctly (see below)
  The code will solve all simulateData when you specify the number of tumor samples, e.g. if the number of tumor samples is 3, it will       solve all the simulateData saved in the folder and save the results for each of the simulateData
    
  * ParentDirectory: specify a directory that contains the LLSolver folder
  * DateFolder: you may do different simulation at different time, this is just for you to record different date, consistent with the one               you set up with DataSimulation.py
  * TumorName: pick a tumor from which you choose the single cell data, now the available single cell data are from GBM07 or GBM33
  * TumorNUmber: specify a tumor sample number so it can get data from that subfolder of the simulation folder
  * reg1: regularization parameter of the penalty in nmf, only will be effecitve if the solver is nmf
  * alpha: regularization parameter of the penalty in ILP, only will be effective if the solver is gurobi
  solver: choose nmf or gurobi to solve the problem
  
  The folder structure will be as following:
  /some/path/to/the/ParentDirectory:
                                  /data/single cell sequencing data
                                  /LLSolver/DataSimulation.py
                                           /DecomposeSolver.py
                                           /....
                                           
                                  /simulation/DateFolder/GBM07/3/<simulateData1>
                                                                /<simulateData2>
                                                                ...
                                                                /<simulateDataN>
                                                                
                                  /test/GTest/<TestCase1>.sh
                                             /<TestCase2>.sh
                                             ...
                                             /<TestCaseN>.sh
                                       /NTest/<TestCase1>.sh
                                             /<TestCase2>.sh
                                             ...
                                             /<TestCaseN>.sh
                                             
                                  /results/DataFolder/GBM07/3/nmf/result_for_simulateData1
                                                                 /result_for_simulateData2
                                                                 ...
                                                                 /result_for_simulateDataN
                                                              
                                                             /gurobi/result_for_simulateData1
                                                                    /result_for_simulateData2
                                                                    ...
                                                                    /result_for_simulateDataN

DecomposeSolver.py includes the following auxiliary functions:

extractValue(directory)
  
> This is to retrieve the data in the replicate in the simulated data generated by DataSimulation.py
> directory: the path to each simulateData
    
SolveDecomposition(AllDataPaths, solver):
  
> This is to solve the problem with desired solver
> AllDataPaths: a list contains the path to each simulateData : [/path/to/simulateData1, /path/to/simulateData2, ...,                                                                          /path/to/simulateDataN ] solver: choose nmf or gurobi to solve the problem

and imports additional functions from NMF_solver.py

In the git repository, the directories
* schwartzlab/test/GTest
* schwartzlab/test/NTest

each contain one example of what <TestCase1>.sh could look like.


----------------------  
NMF_solver.py is called via DecomposeSolver.py, not directly by the user.

NMF_solver.py includes one auxiliary function:

`decompose(tumorMat, major_F, major_cellMat, refer_cellMat, initial_cellMat,
           threshold=1, reg1=0.002, k=6, diploidRatio=0, IncludeInitial=True, seedNum=None)`

----------------------  
GurobiILP_solver.py is called via DecomposeSolver.py, not directly by the user

GurobiILP_solver.py includes the following auxiliary functions:

  * updateProportion(B, Ctotal, cells=None, root=0, vType='C', dirichlet=False, dirA=None, beta=0.0):
  * getDistanceMatrix(Ctotal):
  * updateTree(B, Ctotal, cells=None, alpha=0.1, root=0):
  * calcObjVal(B, F, C, metric='L1'):
  * updateCopyNum(B, F, S, CRefer, cells, alpha=0.1, root=0, vType='C', stopVal=None, Cap=False):
  * makeSurePath(directory):

----------------------
  
SCIP_solver.py is called via DecomposeSolver.py, not directly by the user.

SCIP_solver.py includes the same auxilliary functions as GurobiILP_solver.py.




Reference:
Haoyun Lei, Bochuan Lyu, E. Michael Gertz, Alejandro A. Schaffer, Xulian Shi, Kui Wu, Guibo Li, Liqin Xu, Yong Hou, Michael Dean, Russell Schwartz,
_Tumor Copy Number Deconvolution Integrating Bulk and Single-Cell Sequencing Data_, in preparation.
