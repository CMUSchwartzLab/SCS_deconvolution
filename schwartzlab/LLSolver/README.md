# LLSolver

This directory is part of the repository

https://github.com/CMUSchwartzLab/SCS_deconvolution.git

that contains software to solve several formulations of the problem of
deconvolving bulk tumor data into subpopulations, as described in (Lei
et al., in preparation).  Much of the documentation below assumes that
the user has at least skimmed the manuscript and is familiar with the
terminology therein.

In short, this repository contains implementations of the
"phylogeny-free' and "phylogeny-based" described in the manuscript,
and simulated data that may be used to test these problems.  It also
contains scripts that _could_ be used to produce the semi-simulated
data and figures in the paper.  These scripts are provided as
documentation of what was done, but the semi-simulated data is based
on human subjects data, which cannot be redistributed.  Thus, most
users will probably want to try the deconvolution algorithms on
fully-simulated data.

The main programs that a user may want to try are:

* **SimulateSCS.py** [simulate single-cell data resembling the observed
    data]

* **DataSimulation.py** [simulate many replicates of the hypothetical
    bulk data from the output of SimulateSCS.py or from the observed
    single-cell data]

* **DecomposeSolver.py**  [solve the deconvolution problem]

The DecomposeSolver.py script will typically be the best place to
start.

Programs were written initially by Haoyun Lei and Bochuan Lyu. Programs
were modified by Haoyun Lei and E. Michael Gertz. Programs were tested
by E. Michael Gertz, Haoyun Lei, Bochuan Lyu, and Alejandro Schaffer.

----------
# Installation

Users should clone the git repository, possibly by typing,
```
git clone https://github.com/CMUSchwartzLab/SCS_deconvolution.git
```
Instructions in this README assume a GNU Linux command line or a
Macintosh terminal.  The `git` command above will create a
subdirectory named `SCS_deconvolution`.  The instructions assume that 
the user's current directory is
`SCS_deconvolution/schwartzlab/LLSolver`, in other words the directory that
contains the `README.md` that you are currently reading.

The setup assumes that the user will create another subdirectory
`SCS_deconvolution/schwartzlab/simulation`.  It is inherent to the
code and documentation that the four subdirectories {LLSolver, data,
test, simulation} are parallel, at the same level.  For scripts that
require a path, the user should always specify an absolute path (a
path that starts with `/`).  We cannot include such paths in this 
document, because they are specific to the user's file system.

The scripts in that end in `.py` must be run using python3, not
python2.  The scripts that end in `.sh` are intended to be run using
bash.  Some programs assume the availability of the Gurobi package
(http://www.gurobi.com/downloads/download-center) to solve
optimization problems. We are also testing SCIP
(https://scip.zib.de/index.php#download) as an alternative to Gurobi,
but all analyses in the manuscript that used an optimization package
were done with Gurobi.  Several programs assume the availability of
the python3 numpy and scipy packages.

----------------------
## DecomposeSolver.py

This is the main program that decomposes the bulk tumor data to
resolve the copy number in 9934 loci in each fundamental cell type.
The code will retrieve the simulation data if you set up the input
arguments correctly (see below)

DecomposeSolver.py will use one of three methods, chosen by the user,
to solve the tumor decomposition problem.  The code for these methods
are in

  * NMF_solver.py
  * GurobiILP_solver.py
  * SCIP_solver.py

though the user should not invoke these files directly.
`NMF_solver.py` solves the deconvolution optimization problem in the
phylogeny-free method, while `GurobiILP_solver` and `SCIP_solver`
solve the optimization problem in the phylogeny-based method.

The code will solve all the problmes in the `simulateData1` to
`simulateDataN` subfolders when you specify the number of tumor
samples, e.g. if the number of tumor samples is 3, it will solve all
the simulateData saved in the folder and save the results for each of
the simulateData

  * ParentDirectory: specify a directory that contains the LLSolver folder

  * DateFolder: you may do different simulation at different time,
    this is just for you to record different date, consistent with the
    one you set up with DataSimulation.py

  * TumorName: pick a tumor from which you choose the single cell
    data, now the available single cell data are from GBM07 or GBM33

  * TumorNumber: specify a tumor sample number so it can get data from
    that subfolder of the simulation folder

  * reg1: regularization parameter of the penalty in nmf, only will be
    effective if the solver is nmf

  * alpha: regularization parameter of the penalty in ILP, only will
    be effective if the solver is gurobi solver: choose nmf or gurobi
    to solve the problem

The folder structure will be as following:
```
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
```


--------------
In the git repository, the directories
* schwartzlab/test/GTest
* schwartzlab/test/NTest

each contain one example of what `<TestCase1>.sh` could look like.
Breifly, a call to `DecomposeSolver.py` should look something like
```
python DecomposeSolver.py '/home/some_user/SCS_decomposition/schwartzlab/' 9_28 GBM07 3 0.2 0.0 gurobi
```

where `/home/some_user/SCS_decomposition/schwartzlab/` would be
replaced with an appropriate path on the user's file system.

-------------------
## SimulateSCS.py

The purpose of SimulateSCS.py is to simulate realistic single-cell
data that is similar to the observed single-cell data.  The program
SimulateSCS.py is needed because the observed data are human subjects
data, which cannot be redistributed.  SimulateSCS.py uses summary
statistics from the observed data, which can be found in subdirectory
schwartzlab/data.

SimulateSCS.py takes four arguments:
1. The full path to the data directory
2. The tumor name, which is GBM07 or GBM33

3. The integer mean of the Poisson distribution used in the
simulation (denoted by lam, short for lambda, in the code)

4. The depth of the binary tree of single cells that defines the clone
structure

Example outputs are shown in
* schwartzlab/data/simulated_GBM07_integer_CNV.csv
* schwartzlab/data/simulated_GBM33_integer_CNV.csv

and these were obtained via the calls
* `python SimulateSCS.py 'PATH/TO/Cell/Information' GBM07 115 6`
* `python SimulateSCS.py 'PATH/TO/Cell/Information' GBM33 45 6`

The lam (3rd argument values) of 115 and 45 are recommended for GBM07
and GBM33 respectively.

If the user wishes to change the seed for the random number generator,
then it is necessary to modify SimulateSCS.py and make an assignment
to variable seed at the top of the program.

---------------------
## DataSimulation.py 
  
This program simulates bulk tumor data with desired number of samples
from the single cell sequencing data Simulation is based on the
Geometric Model described in the manuscript and uses a Dirichlet
distribution.  As written, this method requires access to the single
cell data described in the manuscript, which is human subjects data
and cannot be redistributed.  The method is included in this distribution
as supplemental documentation to the manuscript.

The arguments to DataSimulation.py are as follows:

1.  ParentDirectory: specify a directory that contains the LLSolver folder

2.  DateFolder: allows for different simulations run on different
dates to be stored n different subdirectories (a.k.a. folders)

3.  TumorName: pick a tumor from which you choose the single cell
data; the name of the input files should be
../data/<TumorName>_Integer_CNV.csv

4.  tumor_number: choose how many tumor samples you want to simulate,
in the main experiments in the manuscript, we chose 3, 6, or 9 samples

5.  alpha: the alpha parameter for the Dirichlet distribution,
describe in the manuscript

6.  N: the number of simulated replicates (the manuscript uses N=40)

7.  Cap: True of False, whether or not to cap the copy numbers larger
than 10 at 10 (this was set to True for all runs in the manuscript and
must be set to True for the example data provided)

Example calls to DataSimulation.py can be found in the scripts

* schwartzlab/test/runDataSimulationGBM07.sh
* schwartzlab/test/runDataSimulationGBM33.sh

but those are for the observed data.

If one wants to use the simulated data,
`simulated_GBM07_integer_CNV.csv` or `simulated_GBM07_integer_CNV.csv`,
then the third argument, TumorName, should be specified as

> simulated_GBM07

or

> simulated_GBM33

respectively.

To reuse the scripts runDataSimulationGBM07.sh and
runDataSimulationGBM33.sh, the user must also replace the example full
paths give, with other full paths that are suitable for the user's
installation of this repository.

After running DataSimulation.py the folder/subdirectory structure
should be abstractly as follows:
```
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
```

Here, GBM07 is the tumor name selected and could be GBM33
instead). Here, 3 is number of samples.  GTest refers to the
phylogeny-based method, which uses the Gurobi (hence the G) python
library, while Ntest refers to the phylogeny-free method, which is
also called non-negative matrix factorization (NMF, hence the N).

`<TestCase1>.sh` through `<TestCaseN>.sh` are abstract names for the
shell scripts that call DecomposeSolver.py.

Reference:

Haoyun Lei, Bochuan Lyu, E. Michael Gertz, Alejandro A. Schaffer,
Xulian Shi, Kui Wu, Guibo Li, Liqin Xu, Yong Hou, Michael Dean,
Russell Schwartz, _Tumor Copy Number Deconvolution Integrating Bulk
and Single-Cell Sequencing Data_, in preparation.