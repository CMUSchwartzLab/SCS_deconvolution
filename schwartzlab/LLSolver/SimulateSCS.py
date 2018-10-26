import numpy as np
import sys
path = str(sys.argv[1])
TumorType = str(sys.argv[2])
lam = float(sys.argv[3])
depth = int(sys.argv[4])

def GetData(filename):
    return np.genfromtxt('%s%s'%(path, filename), delimiter=',')


def SimulateCN(TumorType='GBM07', lam=115, depth=6, seed=None):
    if TumorType == "GBM07" or TumorType == 'GBM33':
        regions = GetData(TumorType+'RegionDivide.csv')
        regions = regions.astype(int)
        CNrateTotal = GetData(TumorType+'CopyNumberRate.csv')
        POrateTotal = GetData(TumorType+'PositionRate.csv')
    else:
        print("Only GBM07 and GBM33 availalbe")
        return

    if seed != None:
        np.random.seed(seed)

    regionNum = regions.shape[0]-1
    positonNum = POrateTotal.shape[0]
    CNvariations = [0, 1, 3, 4, 5, 6, 7, 8, 9, 10]

    SimSingleCells = np.empty((positonNum, 1))
    for i in range(regionNum):
        SCSinregion = []
        NodeInLv = []
        for ii in range(depth+1):
            NodeInLv.append(np.power(2, ii))
        AllMat = np.zeros((positonNum, np.sum(NodeInLv)))
        AllMat[:] = 2

        for item in NodeInLv[0:-1]:
            flag = item - 1
            node = [(x+flag) for x in range(item)]
            for a in node:
                for b in range(2):
                    #simulate the copy numbers of the two children of one parent node
                    size = np.random.poisson(lam)
                    loci = np.random.choice(
                        positonNum, size=size, p=POrateTotal[:, i], replace=False)
                    VariedCN = np.random.choice(
                        CNvariations, size=size, p=CNrateTotal[:, i])

                    #children node should be same as their parent node
                    AllMat[:, (2*a+b+1)] = AllMat[:, 2*a]
                    #mutate the copy number at some loci to get the children node
                    for c in range(len(loci)):
                        AllMat[loci[c], (2*a+b+1)] = VariedCN[c]
        #random pick the number of cells that equals to number of cells in each region
        cellIndex = np.random.choice(
            np.power(2, depth+1)-1, regions[i+1]-regions[i], replace=False)

        for item in cellIndex:
            SCSinregion.append(AllMat[:, item])
        SCSinregion = np.array(SCSinregion).T
        SimSingleCells = np.append(SimSingleCells, SCSinregion, axis=1)

    np.savetxt('%ssimulated_%s_integer_CNV.csv' % (path, TumorType),
               SimSingleCells[:, 1:], delimiter=',', fmt='%i')


if __name__ == '__main__':
    SimulateCN(TumorType=TumorType, lam=lam, depth=depth)
