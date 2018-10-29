# -*- coding: utf-8 -*-
"""
Created on Thur Sep 27 10:18:58 2018

@author: Haoyun Lei
"""

import ast
import sys

import numpy as np
import scipy.io

import testFunction

ParentDirectory = sys.argv[1]    # Directory to the folder that contains subfolders for
                                 # code, data, result, simulation etc.
DateFolder = str(sys.argv[2])    # specify a folder to save the results
TumorName = str(sys.argv[3])     # pick a tumor, GBM07 or GBM33
tumor_number = int(sys.argv[4])  # the total number of bulk tumor samples
alpha = ast.literal_eval(sys.argv[5]) # Dirichlet distribution parameters
N = int(sys.argv[6])             #how many replicates you want to simulate
Cap = bool(sys.argv[7])          # the largest permitted copy-number; larger
                                 # numbers will be set equal to Cap

#check and/or make the directory to save the simulated data
#save the simulated data with different tumor samples in different data
output_dir = '%ssimulation/%s/%s/%s'%(ParentDirectory, DateFolder, TumorName, str(tumor_number))
testFunction.CheckDirectory(output_dir)

'''
import single cell data
'''

# single cell data should be stored in a folder that 
# is a subfolder of the parent folder where the code stored
# For exmaple: codes are stored in:
#   ~/ParentDirectory/code/
#   Then the data should be saved in:
#   ~/ParentDirectory/data/
        
def ImportSCData(ParentDirectory, TumorName, IntCNV=True, Cap=False):
    '''Read the single-cell copy number data for 'TumorName' from
    'ParentDir'/data.  Return a numpy array with markers as rows
    and single-cell samples as columns.

    If 'IntCNV' is true, read the integer data; otherwise read the
    fractional data.  If Cap is True, make the largest copy number be
    10.

    '''

    if IntCNV:  #get the integer copy number
        DataPath = ParentDirectory + 'data/' +  '%s_integer_CNV.csv' % TumorName 
    else:
        DataPath = ParentDirectory + 'data/' + '%s_fractional_CNV.csv' % TumorName
    
    allSC = np.genfromtxt(DataPath, delimiter=',')
    if Cap:    #Cap the copy number larger than 10 to 10
        allSC[allSC > 10] = 10
    
    return allSC


"""
Set up variables. This can be changed according to the User
Current example is to randomly get 25 cells from each region
make the first 2 of the 25 cells as fundanmental cell components
"""

# In each region there are "dominant" cells and "noise" cells. The
# dominant cells are intended to simulate dominant clones, and will be
# used (stocastically) with larger mixture fractions than the noise
# cells when simulating bulk tumor data.

# cellsList: an array of length 3 with cellsList[i] equal to the
# number of dominant cells to choose from region i
cellsList = [2, 2, 2]

# cellsNoiseList: an array of length 3 with cellsNoiseList[i] equal to the
# number of dominant cells to choose from region i
cellsNoiseList = [23, 23, 23]

# The number of bulk components to simulate for each region, usually
# 1, 2, or 3 (for a total of 3, 6 or 9 bulk components.)
tumor_samples = [int(tumor_number/3), int(tumor_number/3), int(tumor_number/3)]
CSel = [sum(x) for x in zip(cellsList, cellsNoiseList)]
CReferSel = [2, 2, 2]

"""
Function for simulate the desired numbers of data
these data then can be used as the same input to test each method
"""
def SimulateData(N=40, Cap=Cap):    # N is the number of how many data to simulate
    '''Simulate bulk tumor data for N replicates.

    This function relies on several global parameters, including the
    tumor name, the number of bulk components simulated, and the
    relative weights, given as Dirichlet parameters, to use when
    computing the frequencies in the simulation.
    '''
    # The single-cell data, with columns as single cells and
    # rows as genes.
    allSC = ImportSCData(ParentDirectory, TumorName, Cap=Cap)

    for i in range(N):
        # Frequency matrix and Dirichlet parameters.
        [F_orign, dirA] = testFunction.generateF(
            cellsList, cellsNoiseList, tumor_samples, alpha=alpha)
        # Cells used in the simulation and some reference cells.
        [C_orign, C_refer, C_orign_id, C_refer_id] = \
            testFunction.generateCandCRefer(
                allSC, CSel, CReferSel, tumorType=TumorName)

        # Find an initial C to use in (some) deconvolution algorithms.
        init_C = testFunction.initialC(
            allSC, cells=6, usedList=C_orign_id+C_refer_id)

        B = np.dot(C_orign, F_orign)

        # Create copy number and mixture fraction matrices using only
        # 'dominant' cells.
        major_index = testFunction.findMajorIndex(cellsList, cellsNoiseList)
        major_F = F_orign[major_index, :]
        # Normalize the mixture fractions, since we are omitting some cells.
        major_F = major_F / np.sum(major_F, axis=0)
        major_C = C_orign[:, major_index]

        output_file = '%s/simulate_data_%d.mat' % (output_dir, i)
        # Save a dictionary representing the simulation, with the following
        # fields
        scipy.io.savemat(
            output_file,
            {
                'TumorSample':B, # the simulated bulk matrix (an np.array with genes as rows)
                'CTrue':major_C, # copy number of the dominant cells only (transposed)
                'CRefer':C_refer, # several reference cells used in some solvers
                'CIndex':C_orign_id, # the index among single cells of the columns of COrigin
                'dirA':dirA,         # the Dirichlet parameters used to generate B
                'CReferIndex':C_refer_id, # the index among single cells of the columns of CRefer
                'FTrue':major_F,          # Mixture fractions of the dominant cells only.
                'FTrueAll':F_orign, # Mixture fractions used to generated B.
                'CInitial':init_C, # a (possible) starting point for deconvolution algorithms
                'COrigin':C_orign  # Copy numbers used to generate B
            })


if __name__ == '__main__':
    SimulateData(N=N, Cap=Cap)
