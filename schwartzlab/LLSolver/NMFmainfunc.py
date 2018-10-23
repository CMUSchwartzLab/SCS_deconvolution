import decompose
import numpy as np
import time
import sys
import os
import errno
import ast
import scipy.io


date_folder = str(sys.argv[1]) 
tumor_name = str(sys.argv[2])
tumor_number = int(sys.argv[3])
alpha = ast.literal_eval(sys.argv[4])
N = int(sys.argv[5])
k = int(sys.argv[6])
threshod = float(sys.argv[7])
reg1 = float(sys.argv[8])
diploidRatio=float(sys.argv[9])
IncludeInitial = bool(sys.argv[10])

code_path = '/home/haoyunl/Schwartz Lab/Projects/Tumor Matrix Decomposition/Data/7.11_test/'
path = '/home/haoyunl/Schwartz Lab/Projects/Tumor Matrix Decomposition/Data/7.11_test/all_SCS/'
result_path = '/home/haoyunl/Schwartz Lab/Projects/Tumor Matrix Decomposition/Results/' + \
    date_folder +'_'+ str(tumor_number)+'tumor' + '/' + tumor_name + '/'
#sys.path.insert(0, code_path + 'Main_Code/')
import testFunction

if not os.path.exists(result_path):
    try:
        os.makedirs(result_path)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

'''
This is the function to save all the parameters information
into file after each experiments.
'''
def printVar(cellList, cellsNoiseList, CReferSel, DiAalpha, tumor_number, RunTime,
             N=100, k=6, reg1=0.002, threshod=1, diploidRatio=0.0, IncludeInitial=True):
    with open(result_path + 'thre_' + str(threshod) + '_reg1_' + str(reg1) +
              '_variables.txt', 'w') as f:

        f.write('Use %s cell(s) as funtamental components in each region.\n' % cellList)
        f.write('Use %s cell(s) as trivial components in each region.\n' %
                cellsNoiseList)
        f.write('Use %s cell(s) as reference in each region.\n' % CReferSel)
        f.write(
            'The parameter for Dirthlet distrition for mixtrue fraction is %s.\n' % str(alpha))
        f.write("Simulate %s tumor sample(s) in each region, %s tumors in total.\n" %
                (tumor_number / 3, tumor_number))
        f.write('Run %s experiments in total.\n' % N)
        f.write('Use %s as the regularization of the penalty.\n' % reg1)
        f.write('The threshod for convergence is %s.\n' % threshod)
        f.write('If or not use the real cell data as initialization: %s \n' %
                str(IncludeInitial))
        f.write('Containing %s diploid in the initilazation' % diploidRatio)
        f.write("Task run %0.2f hours.\n" % (RunTime / 3600.0))
    f.close()


'''
    Algorithm for decomposing B to C and F
'''


def main_function(N=100):
    start_time = time.time()
    accurate_in_cells = np.zeros((N, k))
    accurate_rows = np.zeros(N)

    
    for z in range(N):
        print("No.", (z + 1), "experiment(s)...")
        '''prepare original cell matrix, referred cell matrix, initial cell matrix
           preare F for each region
        '''
        [F_orign, dirA] = testFunction.generateF(cellList, cellsNoiseList, tumor_samples, alpha=alpha)
        [C_orign, C_refer, C_orign_id, C_refer_id] = \
        testFunction.generateCandCRefer(allSC, CSel, CReferSel, tumorType=tumor_name)
        init_C = testFunction.initialC(allSC, cells=6, usedList = C_orign_id+C_refer_id)
        
        B = np.dot(C_orign, F_orign)
        
        major_index = testFunction.findMajorIndex(cellList, cellsNoiseList)
        major_F = decompose.extraMajorComponent(F_orign.T, major_index)
        major_F = major_F / np.sum(major_F, axis=0)
        major_C = decompose.extraMajorComponent(C_orign, major_index)
        major_C = major_C.T
        
        iter_nn, dist, accuracy, right_row, rmsd_c, rmsd_f, rms_c, rms_f, InferC, InferF = \
        decompose.decompose(B,major_F, major_C, C_refer, init_C, threshod=threshod,\
                            reg1=reg1, k=k, diploidRatio=diploidRatio, IncludeInitial=IncludeInitial)
        
        accurate_in_cells[z, :] = accuracy[:,0]
        accurate_rows[z] = right_row
        '''
        save the result of every experiment in a matlab file,
        so there will be N mat file for in the output.

        alpha: regularization parameter for penalty
        CTrue: the fundamental cells used to simulate tumor, a subset of All True Cells
        CInferred: inferred Cells
        CRefer: the cell used as the observed ones in the penalty 
        CIndex, CReferIndex: The column number of the All True Cells and refer cells
            for all true cells, we selected the first two from each region to simulate tumor
        FTrue, FInferred, FTrueAll: fractions of fundamental cells, inferred cells, and All True Cells
        Accuracy: The accuracy in each cell
        rmsdC: RMSD between true and inferred cell matrix
        rmsdF: RMSD between true and inferred fraction matrix
        rmsInC: RMSD beween true and inferred cell in each cell component
        rmsInF: RMSD beween true and inferred fraction in each cell component
        Step: number of iteration before convergence or break
        totalAcc: accurate row (too strict, just as a reference, might not be useful at all)
        '''
        scipy.io.savemat(result_path + 'result' + str(z) + 'alpha'+ str(reg1) + '.mat', 
                         {'CTrue': major_C, 'CInferred': InferC, 'CRefer': C_refer, 'CIndex':C_orign_id,
                         'CReferIndex': C_refer_id, 'FTrue': major_F, 'FInferred': InferF, 'FTrueAll': F_orign,
                         'Accuracy': accuracy, 'rmsdC': rmsd_c, 'rmsdF': rmsd_f, 'rmsInC':rms_c, 'rmsInF':rms_f,
                         'Step': iter_nn, 'totalAcc': right_row })

    
    RunTime = time.time() - start_time
    print('\n\n')
    print('Use %s cell(s) as funtamental components in each region.' % cellList[0])
    print('Use %s cell(s) as trivial components in each region.' % cellsNoiseList[0])
    print('Use %s cell(s) as reference in each region.' % CReferSel[0]) 
    print('The parameter for Dirthlet distrition for mixtrue fraction is %s.' % str(alpha))
    print("Simulate %s tumor sample(s) in each region, %s tumors in total." % (tumor_number / 3, tumor_number))
    print('Run %s experiments in total.' % N)
    print('Use %s as the regularization of the penalty.' % reg1)
    print('The threshod for convergence is %s.' % threshod)
    print('If or not use the real cell data as initialization: %s ' % str(IncludeInitial))
    print('Containing %s diploid in the initilazation' % diploidRatio)
    print("Task run %0.2f hours." % (RunTime / 3600.0))


    
    printVar(cellList=cellList[0], cellsNoiseList=cellsNoiseList[0], CReferSel=CReferSel[0],
             DiAalpha=alpha, tumor_number=tumor_number, RunTime=RunTime,
             N=N, k=k, reg1=reg1, threshod=threshod, diploidRatio=diploidRatio,
             IncludeInitial=IncludeInitial)
    
if __name__ == '__main__':
    main_function(N=N)
    
    
    
    
