import numpy as np
import glob
import sys
import scipy.io as sio
import time
import testFunction

ParentDirectory = sys.argv[1]   #get to directory that contain other subfolder such as code/ test/ data/ etc.
DateFolder = str(sys.argv[2])   #specify a folder to retrive the data, it also save the result in a folder with corresponding name
TumorName = str(sys.argv[3])    #pick a tumor, GBM07 or GBM33
TumorNumber = int(sys.argv[4])  #specify a tumor sample number so it can get data from that subfolder
reg1 = float(sys.argv[5])       #regularization parameter for nmf
alpha = float(sys.argv[6])      #regularization parameter for gurobi
beta = float(sys.argv[7])       #regularization parameter for SCIP
solver = str(sys.argv[8])       #define what solver to be used, we used NMF and groubi for now, SCIP will be available later
#get the directory of each simulated data
#AllDataPaths = glob.glob(ParentDirectory +'simulation/*.mat')
AllDataPaths = glob.glob('%ssimulation/%s/%s/%s/*.mat'%(ParentDirectory, DateFolder, TumorName, str(TumorNumber)))
def extractValue(directory):
    data = sio.loadmat(directory)
    return data['CIndex'], data['CRefer'], data['CReferIndex'], data['CInitial'],\
            data['CTrue'], data['FTrue'], data['FTrueAll'], data['TumorSample'], data['dirA'], data['COrigin']



def SolveDecomposition(AllDataPaths, solver):
    N = len(AllDataPaths)
    #start_time = time.time()
    #accurate_in_cells = np.zeros((N, k))
    #accurate_rows = np.zeros(N)
    for z in range(N):
        CIndex, CRefer, CReferIndex, CInitial, CTrue, FTrue, FTrueAll, TumorSample, dirA, COrigin=extractValue(AllDataPaths[z])
        TumorNumber = TumorSample.shape[1]
        #choose which solver you want to use:
        if solver == 'nmf':
            import NMF_solver as NS
            #make a directory to save the results
            result_path = '%sresults/%s/%s/%s/%s/' % (
                ParentDirectory, DateFolder, TumorName, TumorNumber, solver)
            testFunction.CheckDirectory(result_path)
            cells = CTrue.shape[1]
            print("No.%s experiment(s) using %s"%(z+1, solver))
            print(" From %s tumor samples, infer %s cells" %
                  (str(TumorNumber), str(cells)))
            iter_nn, dist, accuracy, right_row, rmsd_c, rmsd_f, rms_c, rms_f, InferC, InferF=\
                NS.decompose(TumorSample, FTrue, CTrue, CRefer, CInitial, reg1=reg1, k=cells)
            meanAcc = np.sum(accuracy)/cells
            #accurate_in_cells[z, :] = accuracy[:, 0]
            #accurate_rows[z] = right_row
            sio.savemat(result_path + 'result' + str(z) + 'alpha' + str(reg1) + '.mat',
                             {'CTrue': CTrue, 'CInferred': InferC, 'CRefer': CRefer, 'CIndex': CIndex,
                              'CReferIndex': CReferIndex, 'FTrue':FTrue, 'FInferred': InferF, 'FTrueAll': FTrueAll,
                              'Accuracy': accuracy, 'rmsdC': rmsd_c, 'rmsdF': rmsd_f, 'rmsInC': rms_c, 'rmsInF': rms_f,
                              'Step': iter_nn, 'totalAcc': right_row, 'meanAcc': meanAcc})

            #RunTime = time.time() - start_time
            
            #print("Task run %0.2f hours." % (RunTime / 3600.0))
            #print('\n')

        elif solver=='gurobi':
            import GurobiILP_solver as GS
            #make a directory to save the results
            result_path = '%sresults/%s/%s/%s/%s/' % (
                ParentDirectory, DateFolder, TumorName, TumorNumber, solver)
            testFunction.CheckDirectory(result_path)

            print("No.%s experiment(s) using %s" % (z+1, solver))
            #start_time = time.time()
            CRefer = CRefer.T
            CInitial = CInitial.T
            CTrue = CTrue.T
            COrigin = COrigin.T
            FTrue = FTrue.T
            TumorSample = TumorSample.T

            cellsList = [2, 2, 2]
            cellsNoiseList = [23, 23, 23]
            CSel = np.array(cellsList) + np.array(cellsNoiseList)
            majorIndex = testFunction.findMajorIndex(cellsList, cellsNoiseList)
            dirA = dirA[:, majorIndex]

            cells = CTrue.shape[0]
            cellsObserv = CRefer.shape[0]

            CRefer = np.concatenate(
                (CRefer, 2 * np.ones([1, CRefer.shape[1]])), axis=0)     #add one diploid row to the end as the root
            Ctotal = np.zeros(
                [cells+cellsObserv+1, COrigin.shape[1]], dtype=np.float)
            Ctotal[range(cells),:] = CInitial
            Ctotal[cells:(cells+cellsObserv+1), :] = CRefer

            # #################################################
            # Calculate inferred single-cell components
            oldObj = [0, 0, 0]
            thresholdI = 10 ** (-4)
            step = 1

            Cprev = np.matrix(Ctotal, dtype=np.float)
            print(" regularization parameter=%s, from %s tumor samples, infer %s cells." % (str(alpha), str(TumorNumber),str(cells)))
            while(1):
                #print("Step:", step)
                [F, objVal1] = GS.updateProportion(
                    TumorSample, Ctotal, cells, root=cells+cellsObserv, dirA=dirA)
                [S, objVal2] = GS.updateTree(
                    TumorSample, Ctotal, cells, alpha=alpha, root=cells+cellsObserv)
                step += 1
                [CUnknown, objVal] = GS.updateCopyNum(TumorSample, F, S, CRefer, cells, alpha=alpha, root=cells+cellsObserv,
                                                vType='I', Cap=True)
                Ctotal[0:cells, :] = CUnknown
                change = abs(oldObj[2] - objVal)
                change1 = abs(oldObj[0] - objVal1)
                change2 = abs(oldObj[1] - objVal2)
                #print('objVal:', objVal)
                oldObj[2] = objVal
                oldObj[0] = objVal1
                oldObj[1] = objVal2
                if (change < thresholdI or step > 100):
                    break
            acc = testFunction.calcAccuracy(CUnknown, CTrue, CellsInCol=False)
            [CUnknown, order] = testFunction.arrangeC(
                CUnknown, CTrue, CellsInCol=False)

            totalAcc = testFunction.calcAccuracyByRow(
                CUnknown, CTrue, CellsInCol=False)
            F = F[:, order]
            F = np.matrix(F)
            
            rmsdC = testFunction.calcRMSD(CUnknown, CTrue)
            rmsdF = testFunction.calcRMSD(F, FTrue)
            meanAcc = np.sum(acc)/cells
            rms_c = testFunction.calcRMSInCell(CUnknown, CTrue, CellsInCol=False)
            rms_f = testFunction.calcRMSInCell(F, FTrue, Cell=False, CellsInCol=False)

            sio.savemat(result_path + 'result' + str(z) + 'alpha' + str(alpha) + '.mat',
                        {'meanAcc': meanAcc, 'Accuracy': acc, 'CTrue': CTrue.T, 'CRefer': CRefer[0:6,:].T,
                        'totalAcc': totalAcc, 'CInferred': CUnknown.T, 'FTrueAll': FTrueAll, 'FTrue': FTrue.T,
                         'FInferred': F.T, 'Step': step, 'rmsdC': rmsdC, 'rmsdF': rmsdF, 'CIndex': CIndex,
                         'CReferIndex': CReferIndex, 'rmsInC': rms_c, 'rmsInF': rms_f, 'TreeStr': S})

        elif solver=='scip':
            import SCIP_solver as SP
            #make a directory to save the results
            result_path = '%sresults/%s/%s/%s/%s/' % (
                ParentDirectory, DateFolder, TumorName, TumorNumber, solver)
            testFunction.CheckDirectory(result_path)
            print("No.%s experiment(s) using %s" % (z+1, solver))

            CRefer = CRefer.T
            CInitial = CInitial.T
            CTrue = CTrue.T
            COrigin = COrigin.T
            FTrue = FTrue.T
            TumorSample = TumorSample.T

            cellsList = [2, 2, 2]
            cellsNoiseList = [23, 23, 23]
            CSel = np.array(cellsList) + np.array(cellsNoiseList)
            majorIndex = testFunction.findMajorIndex(cellsList, cellsNoiseList)
            dirA = dirA[:, majorIndex]

            cells = CTrue.shape[0]
            cellsObserv = CRefer.shape[0]

            CRefer = np.concatenate(
                (CRefer, 2 * np.ones([1, CRefer.shape[1]])), axis=0)  # add one diploid row to the end as the root
            Ctotal = np.zeros(
                [cells+cellsObserv+1, COrigin.shape[1]], dtype=np.float)
            Ctotal[range(cells), :] = CInitial
            Ctotal[cells:(cells+cellsObserv+1), :] = CRefer

            # #################################################
            # Calculate inferred single-cell components
            oldObj = [0, 0, 0]
            thresholdI = 10 ** (-4)
            step = 1

            Cprev = np.matrix(Ctotal, dtype=np.float)
            print(" regularization parameter=%s, from %s tumor samples, infer %s cells." % (
                str(alpha), str(TumorNumber), str(cells)))
            while(1):
                #print("Step:", step)
                [F, objVal1] = SP.updateProportion(
                    TumorSample, Ctotal, cells, root=cells+cellsObserv, dirA=dirA, beta=beta)
                [S, objVal2] = SP.updateTree(
                    TumorSample, Ctotal, cells, alpha, root=cells+cellsObserv)
                step += 1
                [CUnknown, objVal] = SP.updateCopyNum(TumorSample, F, S, CRefer, cells, beta, root=cells+cellsObserv,
                                                vType='I', Cap=True)
                Ctotal[0:cells, :] = CUnknown
                change = abs(oldObj[2] - objVal)
                change1 = abs(oldObj[0] - objVal1)
                change2 = abs(oldObj[1] - objVal2)
                #print('objVal:', objVal)
                oldObj[2] = objVal
                oldObj[0] = objVal1
                oldObj[1] = objVal2
                if (change < thresholdI or step > 100):
                    break

            acc = testFunction.calcAccuracy(CUnknown, CTrue, CellsInCol=False)
            [CUnknown, order] = testFunction.arrangeC(
                CUnknown, CTrue, CellsInCol=False)

            totalAcc = testFunction.calcAccuracyByRow(
                CUnknown, CTrue, CellsInCol=False)
            F = F[:, order]
            F = np.matrix(F)

            rmsdC = testFunction.calcRMSD(CUnknown, CTrue)
            rmsdF = testFunction.calcRMSD(F, FTrue)
            meanAcc = np.sum(acc)/cells
            rms_c = testFunction.calcRMSInCell(
                CUnknown, CTrue, CellsInCol=False)
            rms_f = testFunction.calcRMSInCell(
                F, FTrue, Cell=False, CellsInCol=False)

            sio.savemat(result_path + 'result' + str(z) + 'alpha' + str(beta) + '.mat',
                        {'meanAcc': meanAcc, 'Accuracy': acc, 'CTrue': CTrue.T, 'CRefer': CRefer[0:6, :].T,
                         'totalAcc': totalAcc, 'CInferred': CUnknown.T, 'FTrueAll': FTrueAll, 'FTrue': FTrue.T,
                         'FInferred': F.T, 'Step': step, 'rmsdC': rmsdC, 'rmsdF': rmsdF, 'CIndex': CIndex,
                         'CReferIndex': CReferIndex, 'rmsInC': rms_c, 'rmsInF': rms_f, 'TreeStr': S})

        else:
            print('Solver Not Available, please choose nmf, gurobi or scip')
            pass

                
if __name__ == '__main__':
    start_time = time.time()
    SolveDecomposition(AllDataPaths=AllDataPaths, solver=solver)
    RunTime = time.time() - start_time
    print("\nTotal Runtime is %0.2f hours." % (RunTime / 3600.0))
    print("Simulate tumor from %s"%TumorName) 
    print("Simulate %s tumor sample(s) in each region, %s tumors in total." % (TumorNumber/3, TumorNumber))
    print("Using %s as solver" % solver)
    if solver == 'gurobi':
        print("Regularization paramter of penalty is %s" % str(alpha))
    elif solver == 'nmf':
        print("Regularization paramter of penalty is %s" % str(reg1))
    elif solver == 'scip':
        print("Regularization paramter of penalty is %s" % str(beta))
