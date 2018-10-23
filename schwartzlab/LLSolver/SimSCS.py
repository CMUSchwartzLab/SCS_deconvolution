#%%
import numpy as np
import sys
path = str(sys.argv[1])
TumorType = str(sys.argv[2])
lam = float(sys.argv[3])
depth = int(sys.argv[4])


def SimSCS(TumorType='GBM07', lam=50, depth=7, seed=None):
    if TumorType == "GBM07":
        regions = [0, 52, 127, 195]
        SCS = np.genfromtxt(path + 'GBM07_integer_CNV.csv', delimiter=',')
    elif TumorType == "GBM33":
        regions = [0, 56, 132, 198]
        SCS = np.genfromtxt(path+'GBM33_integer_CNV.csv', delimiter=',')
    else:
        print("Only GBM07 and GBM33 availalbe")
        return

    if seed != None:
        np.random.seed(seed)

    SCS[SCS > 10] = 10
    regionNum = len(regions) - 1
    CNvariations = list(range(11))

    SimSingleCells = np.empty((SCS.shape[0], 1))
    for i in range(regionNum):
        SCSinregion = []
        subSCS = SCS[:, regions[i]:regions[i+1]+1]

        NodeInLv = []
        for ii in range(depth+1):
            NodeInLv.append(np.power(2, ii))
        AllMat = np.zeros((SCS.shape[0], np.sum(NodeInLv)))
        AllMat[:] = 2

        #get the frequency of each possible copy number in the original data
        CNrate = []
        for item in CNvariations:
            CNrate.append(np.count_nonzero(subSCS == item) /
                          float(subSCS.shape[0]*subSCS.shape[1]))

        #get the frequency of each position
        POrate1 = []
        for z in range(SCS.shape[0]):
            POrate1.append(np.count_nonzero(
                subSCS[z, :] != 2)/float(subSCS.shape[1]))
        POrate = [x/float(np.sum(POrate1)) for x in POrate1]

        #simulate the copy number from the diploid root
        #through a binary tree
        for item in NodeInLv[0:-1]:
            flag = item - 1
            node = [(x+flag) for x in range(item)]
            for a in node:
                for b in range(2):
                    #simulate the copy numbers of the two children of one parent node
                    size = np.random.poisson(lam)
                    loci = np.random.choice(
                        SCS.shape[0], size=size, p=POrate, replace=False)
                    VariedCN = np.random.choice(11, size=size, p=CNrate)

                    #children node should be same as their parent node
                    AllMat[:, (2*a+b+1)] = AllMat[:, 2*a]
                    #mutate the copy number at some loci to get the children node
                    for c in range(len(loci)):
                        AllMat[loci[c], (2*a+b+1)] = VariedCN[c]

        #random pick the number of cells that equals to number of cells in each region
        cellIndex = np.random.choice(
            np.power(2, depth), regions[i+1]-regions[i], replace=False)
        for item in cellIndex:
            SCSinregion.append(AllMat[:, item])
        SCSinregion = np.array(SCSinregion).T
        SimSingleCells = np.append(SimSingleCells, SCSinregion, axis=1)
        
    np.savetxt(path+'simulated%s_integer_CNV.csv'%TumorType, SimSingleCells, delimiter=',')


if __name__ == '__main__':
    SimSCS(TumorType=TumorType, lam=lam, depth=depth)
