import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import sys, glob
import testFunction
import ast
import networkx as nx
import pandas as pd
import scipy.io as sio

# get to directory that contain other subfolders such as code test data etc.
ParentDirectory = sys.argv[1]                  # will create a folder named figures under this directory
DateFolder = str(sys.argv[2])                  # specify a folder to retrieve the results
TumorName = str(sys.argv[3])                   # pick a tumor, GBM07 or GBM33
TumorNumbers = ast.literal_eval(sys.argv[4])   # a list of different number of tumor samples
solvers = ast.literal_eval(sys.argv[5])        # a list of different solver used 
SavedFolder = str(sys.argv[6])                 # # specify a folder to save the figures

AllDataPaths = glob.glob('%sresults/%s/%s/*.mat' %
                         (ParentDirectory, DateFolder, TumorName))

#check and/or make the directory to save the figures
testFunction.CheckDirectory('%sfigures/%s/%s/' %
                            (ParentDirectory, SavedFolder, TumorName))
#set the style of the figures
sns.set(style="ticks", palette="pastel")

#extract the result of the test, such average accuarcy, RMSD etc
def extractValue(directory, key):
    data = sio.loadmat(directory)
    return data[key]

"""
get the result for different solver
then calculathe the average from all the test cases
saved in an array, row is different tumor samples, 
                column is average result in different regularization parameter 
"""
def GetAllData(solver, Date, TumorNumbers, RegPara, key):
    # get the average data from the result
    DateFolder=Date
    Data = []
    for i in TumorNumbers:
        AverValue = []
        for item in RegPara:
            SingleValue = []
            Paths = glob.glob('%sresults/%s/%s/%s/%s/result*alpha%s.mat' % (ParentDirectory, DateFolder,
                                                                 TumorName, i, solver, item))
            for j in Paths:
                SingleValue.append(extractValue(directory=j, key=key))
            AverValue.append(np.average(SingleValue))
        Data.append(AverValue)
    Data = np.array(Data)
    Data = pd.DataFrame(Data, columns=RegPara)
    Data['Tumor Samples'] = TumorNumbers
    return Data

"""
Function to get the average data if guess a diploid matrix
get the average measurement or the measurement in inividual column
"""
def GetCtrData(solver, Date,TumorNumbers, RegPara, key, average=True):
    DateFolder = Date
    #returned result is single average value: either accuracy or RMSD
    #the value will be used as a dashed line in the bar chart
    exData = []
    for i in TumorNumbers:
        AverValue = []
        for item in RegPara:
            SingleValue = []
            Paths = glob.glob('%sresults/%s/%s/%s/%s/result*alpha%s.mat'
                              % (ParentDirectory, DateFolder, TumorName, i, solver, item))
            for j in Paths:
                TrueCell = extractValue(j, 'CTrue')
                m, n = TrueCell.shape
                DiploidMat = np.zeros((m, n))
                DiploidMat[:] = 2
                if key == 'meanAcc':
                    DiploidAcc = testFunction.calcAccuracy(TrueCell, DiploidMat)
                elif key == 'rmsdC':
                    DiploidAcc = testFunction.calcRMSD(TrueCell, DiploidMat)
                if average:
                    SingleValue.append(np.average(DiploidAcc))
                else:
                    SingleValue.append(DiploidAcc)
            AverValue.append(SingleValue)
        exData.append(AverValue)
    if average:
        return np.average(np.array(exData))
    else:
        exData = np.average(np.array(exData), axis=0)
        exData = np.average(exData, axis=0)
        return exData[:,:,0]

"""
grouped barplot: 
    a grouped barplot to show the average performance 
    in random initialization and w/ penalty
"""
def GroupBarPlot(solver, groupedData, key, kw, DipData=None, NoInitialData=None):
    #merge the dataframe by pivoiting the TumorSamples columns
    df = pd.melt(groupedData, id_vars=['Tumor Samples'],
                 var_name='Penalty', value_name='Average %s' % key)
    #plt.figure(figsize=(8, 6))
    b = sns.barplot(x='Tumor Samples', y='Average %s' % key,
                    hue='Penalty', data=df)
    b.set_xlabel('Tumor Samples', fontsize=16)
    b.set_ylabel('Average %s' % key, fontsize=16)
    if DipData != None:
        # plot a dashed line of average  result of Diploid guess
        plt.axhline(DipData, color='k', linestyle='--')
    if NoInitialData != None:
        # plot a dashed line of average  result of random initialization
        plt.axhline(NoInitialData, color="#e74c3c", linestyle='--')
    plt.title('Average %s of %s' % (key, kw), fontsize=18)
    plt.legend(title='ReguPara', loc='center right',
               bbox_to_anchor=(1.22, 0.65))
    #plt.show()
    plt.savefig('%sfigures/%s/%s/%s_%s_%s_barplot.png' %
                (ParentDirectory, SavedFolder, TumorName, solver, key, kw),
                bbox_inches='tight', dpi=600)



def DataIndexforCNV(solver, TumorName, TumorNumbers, RegPara, key='Accuracy', rule='good'):
    PathsInds = []
    for i in TumorNumbers:
        Temp = []
        Paths = glob.glob('%sresults/%s/%s/%s/%s/result*alpha%s.mat' %
                          (ParentDirectory, DateFolder, TumorName, i, solver, RegPara))
        for j in range(len(Paths)):
            if key == 'Accuracy':
                Temp.append(extractValue(Paths[j], key)[:, 0])
            elif key == 'rmsInC':
                Temp.append(extractValue(Paths[j], key)[0, :])
        Temp = np.array(Temp)
        #get the index of the maxium or minimum value in the 2-D array
        if key == 'Accuracy':
            if rule == 'good':
                Ind = np.unravel_index(np.argmax(Temp, axis=None), Temp.shape)
            else:
                Ind = np.unravel_index(np.argmin(Temp, axis=None), Temp.shape)
        elif key == 'rmsInC':
            if rule == 'good':
                Ind = np.unravel_index(np.argmin(Temp, axis=None), Temp.shape)
            else:
                Ind = np.unravel_index(np.argmax(Temp, axis=None), Temp.shape)
        PathIndex = tuple((Paths[Ind[0]], Ind[1]))
        PathsInds.append(PathIndex)
        
    return PathsInds

def DataforCNplot(solver, RegPara, key='Accuracy', rule='good'):
    PathAndInd = DataIndexforCNV(solver=solver, TumorName='GBM07', 
                    TumorNumbers=TumorNumbers, RegPara=RegPara,
                    key=key, rule=rule)
    CTrueAll = []
    CInferredAll = []
    for item in PathAndInd:
        CTrue = extractValue(item[0], 'CTrue')
        CInferred = extractValue(item[0], 'CInferred')
        CTrueAll.append(CTrue)
        CInferredAll.append(CInferred)
    
    return CTrueAll, CInferredAll

"""
Function to get the tree matrix from desired files
"""
def DataforTree(solver, RegPara, key='Accuracy', rule='good'):
    PathAndInd = DataIndexforCNV(solver=solver, TumorName='GBM07',
                                 TumorNumbers=TumorNumbers, RegPara=RegPara,
                                 key=key, rule=rule) 
    Tree = []
    for item in PathAndInd:
        Tree.append(extractValue(item[0], 'TreeStr'))
    return np.array(Tree)
"""
gene distribution plot to show the real CN
in loci, a good example and a bad example
"""
def PlotCNV(dataLst1, dataLst2, solver):
    fig, axs = plt.subplots(1, len(dataLst1), figsize=(
        22, 4), sharex=True, sharey=True)
    for i in range(len(dataLst1)):
        #dataLst2[i][dataLst2[i] > 10] = 10
        axs[i].plot(dataLst1[i], lw=8, color='C1', label='True Copy Number')
        axs[i].plot(dataLst2[i], ls='--', color='k',
                    label='Inferred Copy Number')
        axs[i].set_title('%s tumor samples' % TumorNumbers[i], fontsize=19)
    handles, labels = axs[0].get_legend_handles_labels()
    fig.legend([handles[0], handles[-1]], [labels[0], labels[-1]],
               loc='center right', bbox_to_anchor=(1, 0.07))

    #fig.suptitle('Copy Number Comparison in Infered and True Cell',  y=1.05, fontsize=14)
    fig.text(0.5, -0.06, 'Genomic Loci', ha='center', fontsize=17)
    fig.text(-0.005, 0.5, 'Copy Number', va='center',
             rotation='vertical', fontsize=17)
    plt.tight_layout()
    plt.savefig('%sfigures/%s/%s/%s_CNV.png' %
                (ParentDirectory, SavedFolder, TumorName, solver),
                bbox_inches='tight', dpi=600)


"""
Function to prepare the data for comparison
each bar will be plot with the best parameter
"""
def BestResult(solver, Date, TumorNumbers, RegPara, key, average=True):
    DateFolder = Date
    Data = []
    for i in range(len(TumorNumbers)):
        SingleData = GetAllData(solver, DateFolder,[TumorNumbers[i]], [RegPara[i]], key).values
        Data.append(SingleData[0])
    Data = pd.DataFrame(np.array(Data), columns=[
                        'average values', 'Tumor Samples'])
    return Data


"""
Function to compare the results between solvers
Now only works for NMG and Gurobi
"""
def BarPlotOfSolver(NMFdata, gurobiData, value_name, kw):
    #combine the two dataframes to one
    m, n = NMFdata.shape
    NMFdata['solver'] = ['phylogeny-free'] * m
    gurobiData['solver'] = ['phylogeny-based'] * m
    data = pd.concat([NMFdata, gurobiData])
    data['Tumor Samples'] = data['Tumor Samples'].astype(int)
    #merge the dataframe by pivoiting the TumorSample columns
    df = pd.melt(data, id_vars=['Tumor Samples',
                                'solver'], value_name=value_name)
    #plt.figure(figsize=(10,10))
    sns.barplot(x='Tumor Samples', y=value_name, hue='solver', data=df)
    plt.ylabel('Average %s' % value_name, fontsize=16)
    plt.xlabel('Tumor Samples', fontsize=16)
    plt.title("%s %s" % (kw, value_name), fontsize=18)
    plt.legend(title='Method', loc='center right', bbox_to_anchor=(1.4, 0.65))
    #plt.show()
    plt.savefig('%sfigures/%s/%s/%s_%s_compare.png' %
                (ParentDirectory, SavedFolder, TumorName, value_name, kw),
                bbox_inches='tight', dpi=600)


"""
Function to plot the box of guessing a diploid cell matrix
"""  
def DiplodBoxPlot(IndividualDiploid, AverageDiploid, kw):
    df = pd.DataFrame(
        IndividualDiploid, columns=["cell 1", "cell 2", "cell 3", "cell 4", "cell 5", "cell 6 "])

    sns.boxplot(data=df)
    plt.axhline(y=AverageDiploid, color='k', ls='--')
    plt.ylabel('%s' % kw)
    plt.title('%s in Each Cell Component' % kw, fontsize=14)
    plt.savefig('%sfigures/%s/%s/diploidBOX_%s.png' %
                (ParentDirectory, SavedFolder, TumorName, kw),
                bbox_inches='tight', dpi=600)

"""
Function of boxplot for the measurement in each cell component
for the tested cases
"""
def GetDataInCell(solver, Date, TumorNumbers, RegPara, key):
    DateFolder = Date
    TotalValue = pd.DataFrame()
    for i in TumorNumbers:
        TumorValue = pd.DataFrame()
        for item in RegPara:
            ReguValue = []
            Paths = glob.glob('%sresults/%s/%s/%s/%s/result*alpha%s.mat' % (ParentDirectory, DateFolder,
                                                                            TumorName, i, solver, item))
            for j in Paths:
                ReguValue.append(extractValue(
                    directory=j, key=key).reshape(1, 6)[0])
            ReguValue = pd.DataFrame(np.array(ReguValue), columns=[
                "cell 1", "cell 2", "cell 3", "cell 4", "cell 5", "cell 6"])
            ReguValue['ReguPara'] = item
            TumorValue = TumorValue.append(
                ReguValue, ignore_index=True, sort=False)
            TumorValue['Tumor Samples'] = i
        TotalValue = TotalValue.append(TumorValue, ignore_index=True, sort=False)
    return TotalValue


def BoxPlotInCell(solver, dataframe, kw1, kw2):
    df = pd.melt(dataframe, id_vars=[
                 'Tumor Samples', 'ReguPara'], var_name='Cell', value_name='Average %s' % kw1)
    unique_samples = np.unique(df['Tumor Samples'])
    fig, axs = plt.subplots(1, unique_samples.shape[0], figsize=(
        20, 4), sharex=True, sharey=True)
    for i in range(unique_samples.shape[0]):
        ax = sns.boxplot(x='Cell', y='Average %s' % kw1, hue='ReguPara',
                         data=df[df['Tumor Samples'] == unique_samples[i]], ax=axs[i])
        ax.set_title('%s tumor samples' % str(unique_samples[i]), fontsize=19)
        ax.legend().set_visible(False)

    for ax in axs:
        ax.set_xlabel('')
        ax.set_ylabel('')

    handles, labels = axs[-1].get_legend_handles_labels()
    fig.suptitle(
        '%s of %s in Each Cell Component' % (kw1, kw2), y=1.05, fontsize=21)
    fig.text(0.5, -0.005, 'Cell Components', ha='center', fontsize=17)
    fig.text(-0.005, 0.5, '%s' % kw1, va='center',
             rotation='vertical', fontsize=17)
    axs[-1].legend(handles, labels, title='ReguPara',
                   loc='center right', bbox_to_anchor=(1.2, 0.5))
    plt.tight_layout()
    plt.savefig('%sfigures/%s/%s/%s_%s_%s_boxpplot.png' %
                (ParentDirectory, SavedFolder, TumorName, solver, kw1, kw2),
                bbox_inches='tight', dpi=600)

"""
Function to boxplot the comparision between Tree-free and tree-based methods
"""
def GetBestDataForBoxComparison(solver, Date, TumorNumbers, RegPara, key):
    DateFolder = Date
    TumorValue = pd.DataFrame()
    for i in range(len(TumorNumbers)):
        ReguValue = []
        Paths = glob.glob('%sresults/%s/%s/%s/%s/result*alpha%s.mat' % 
                            (ParentDirectory, DateFolder, TumorName, TumorNumbers[i], solver, RegPara[i]))
        for j in Paths:
            ReguValue.append(extractValue(
                directory=j, key=key).reshape(1, 6)[0])
        ReguValue = pd.DataFrame(np.array(ReguValue), columns=[
            "cell 1", "cell 2", "cell 3", "cell 4", "cell 5", "cell 6"])
        ReguValue['Tumor Samples'] = TumorNumbers[i]
        TumorValue = TumorValue.append(
            ReguValue, ignore_index=True, sort=False)
        
    return TumorValue


def BoxPlotComparison(NMFdata, gurobiData, kw1, kw2):
    m, n = NMFdata.shape
    NMFdata['solver'] = ['phylogeny-free'] * m
    gurobiData['solver'] = ['phylogeny-based'] * m
    data = pd.concat([NMFdata, gurobiData])
    df = pd.melt(data, id_vars=[
                 'Tumor Samples', 'solver'], var_name='Cell', value_name='Average %s' % kw1)
    unique_samples = np.unique(df['Tumor Samples'])
    fig, axs = plt.subplots(1, unique_samples.shape[0], figsize=(
        20, 4), sharex=True, sharey=True)
    for i in range(unique_samples.shape[0]):
        ax = sns.boxplot(x='Cell', y='Average %s' % kw1, hue='solver',
                         data=df[df['Tumor Samples'] == unique_samples[i]], ax=axs[i])
        ax.set_title('%s tumor samples' % str(unique_samples[i]), fontsize=19)
        ax.legend().set_visible(False)

    for ax in axs:
        ax.set_xlabel('')
        ax.set_ylabel('')

    handles, labels = axs[-1].get_legend_handles_labels()
    fig.suptitle(
        '%s of %s in Each Cell Component' % (kw1, kw2), y=1.05, fontsize=21)
    fig.text(0.5, -0.005, 'Cell Components', ha='center', fontsize=17)
    fig.text(-0.005, 0.5, '%s' % kw1, va='center',
             rotation='vertical', fontsize=17)
    axs[-1].legend(handles, labels, title='ReguPara',
                   loc='center right', bbox_to_anchor=(1.2, 0.5))
    plt.tight_layout()
    plt.savefig('%sfigures/%s/%s/%s_%s_compareboxpplot.png' %
                (ParentDirectory, SavedFolder, TumorName, kw1, kw2),
                bbox_inches='tight', dpi=600)

"""
Functions to plot all the tree structure
"""

#This is to take the Tree Matrix in the result as input 
#and output the tree
def plotTree(S, root, filename='example.png'):
    G = nx.DiGraph()
    fromList = []
    toList = []

    G.add_node(root)
    for i in range(S.shape[0]):
        if (i != root):
            G.add_node(i)

    for i in range(S.shape[0]):
        for j in range(S.shape[1]):
            if (S[i, j] != 0):
                fromList.append(i)
                toList.append(j)

    for i in range(len(fromList)):
        G.add_edge(fromList[i], toList[i], weight=1)
    elarge = [(u, v) for (u, v, d) in G.edges(data=True)]
    p = nx.drawing.nx_pydot.to_pydot(G)

    p.write_png(filename)

#This is to plot all the tree in the result, 
# and store them as indepent .png files
def PlotAllTrees(SavePath, TumorNumber, ReguPara, solver='gurobi'):
    for item in ReguPara:
        Paths = glob.glob('%sresults/%s/%s/%s/%s/result*alpha%s.mat' % (ParentDirectory, DateFolder,
                                                                        TumorName, TumorNumber, solver, item))
        testFunction.CheckDirectory('%s%s/%s/' % (SavePath, TumorNumber, item))
        for i in range(len(Paths)):
            plotTree(S=extractValue(Paths[i], 'TreeStr'), root=12,
                                    filename='%s%s/%s/' % (SavePath, TumorNumber, item)+'%s.png' % i)
