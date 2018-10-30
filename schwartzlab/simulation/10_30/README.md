In this specific DateFolder, we included small amount of replicates of simulated Bulk data. Each tumor directory contain 3, 6, 9, 33, 99 samples case, each sample case contain replicates for that sample case.

If the user works on bulk data simulated from Simulated SCS data and specify the TumorName with the prefix `simulated_`, the directory would be like:
```
simulation/DateFolder/simulated_GBM07/3/simulate_data_0
                                       /simulate_data_1
                                     /6/simulate_data_0
                                       /simulate_data_1
                                     /9/simulate_data_0
                                       /simulate_data_1
                                    ...
                     /simulated_GBM33/3/simulate_data_0
                                       /simulate_data_1
                                     /6/simulate_data_0
                                       /simulate_data_1
                                     /9/simulate_data_0
                                       /simulate_data_1
                                    ...                               
```
If the user has the real observed SCS data, the TumorName then should be specifed as `GBM07` or `GBM33`, and the directory structure would be like:
```
simulation/DateFolder/GBM07/3/...
                           /6/...
                           /9/...
                           ....
                     /GBM33/3/...
                           /6/...
                           /9/...
                           ....
```

While we cannot redistribute human subject data, we created `GBM07` and `GBM33` with empty subdectories for user's reference. Also, it would be ok to have `simulated_GBM07` and `GBM07` in the same directory, we expect the user to specify the corret argument when calling DecomposeSolver.py 
Much more information can be found in
[../../LLSolver/README.md](../../LLSolver/README.md)