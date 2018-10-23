import numpy as np
import sys
from scipy.optimize import linear_sum_assignment
import random
import testFunction


def decompose(tumorMat, major_F, major_cellMat, refer_cellMat, initial_cellMat,
              threshold=1, reg1=0.002, k=6, diploidRatio=0, IncludeInitial=True, seedNum=None):
    print(" diploidRatio = %s, regularization 1 = %s, real cell as initial? %s " % (
        diploidRatio, reg1, IncludeInitial))
    m, n = tumorMat.shape
    B = tumorMat
    '''
    Initialize cell component
    If we decide to refer to observed cell, we set w0 to Initial Cells
    Else, we random select integer from 0 to 10 to fill the matrix
        or we gradually fill the matrix with diploid.
    '''
    if seedNum != None:
        np.random.seed(int(seedNum))   #seed the random generator.

    if IncludeInitial:
        w0 = initial_cellMat
    else:
        diplodRow = int(m * diploidRatio)
        RowInd = np.random.choice(m, diplodRow, replace=False) #will not pick the row index that already picked
        #RowInd = random.sample(list(range(m)), diplodRow)
        if diploidRatio == 0:
            #just simple random entries from 1-10, including 2
            w0 = np.random.randint(11, size=(m, k))
        else:
            w0 = np.zeros((m,k))
            for i in range(m):
                for j in range(k):
                    #no entries will be 2
                    w0[i,j] = np.random.choice([k for k in range(11) if k!=2],1)[0]

            for item in RowInd:
                #change the whole row to 2 if the row was picked
                w0[item, :] = 2

    h0 = np.random.rand(k, n)
    h0 = h0 / np.sum(h0, axis=0)

    flag = 1000

    dnorm = 1
    dnorm0 = 0

    # What follows here is a straight NMF solver, without a constraint
    # that w0 be integral or that the columns of h sum to one.
    iter_nn = 0
    dist = []
    while flag > threshold:
        #starting condition
        dnorm = np.linalg.norm(B - np.dot(w0, h0))

        # if the updated L2 norm began to diverge, stop
        if (dnorm0 - dnorm > 0 and iter_nn > 10) or (iter_nn > 100000):
            break
        #update h
        #np.spacing is to add small value to avoid being divided by 0
        numer = np.transpose(w0).dot(B)
        h = h0 * \
            (numer / (np.dot(np.dot(w0.T, w0), h0) + np.spacing(numer)))
        h[h < 0] = 0
        #normalize columns of h
        h = h / np.sum(h, axis=0)

        #update w
        numer = B.dot(h.T)
        w = w0 * (numer / (np.dot(w0, np.dot(h, h.T)) + reg1 * (w0 - refer_cellMat)
                           + np.spacing(numer)))
        w[w < 0] = 0
        #make the CN to integer and cap the large CN to 10
        w = np.round(w)
        w[w > 10] = 10

        #update w0 and h0 for new iteration
        w0 = w
        h0 = h

        dnorm0 = np.linalg.norm(B - np.dot(w0, h0))
        flag = np.linalg.norm(B - np.dot(w0, h0))

        #track the update
        iter_nn += 1
        dist.append(dnorm0)

    # Normalize the columns of h0 to sum to one.
    #h0 = h / np.sum(h, axis=0)
    # Do a *new* least-squares estimate of w3
    #w3 = np.transpose(np.linalg.lstsq(h0.T, B.T, rcond=-1)[0])
    # and make the components positive.
    #w3[w3 < 0] = 0

    # rearrange the columns of w3/rows of h(0?) so that they best match
    # the true solution.
    arr_list = testFunction.arrangeC(w0, major_cellMat)
    # *now* round the rearranged columns (why now?)
    w4 = np.round(arr_list[0])
    #cap the copy number larger than 10 to 10
    #w4[w4>10] = 10
    h_updated = []
    for item in arr_list[1]:
        # Should be h0 here? h was not normalized.
        #I think I made a typo here, should be h0, not h --Haoyun
        h_updated.append(h0[item, :])
    h_updated = np.array(h_updated)
    # h0 = np.transpose(np.linalg.lstsq(w4, B, rcond=-1)[0])
    # h0 = h0 / np.sum(h0, axis=0)

    '''compare the inferred result with the original major compoment
       information
    '''
    average_accuracy = testFunction.calcAccuracy(w4, major_cellMat)
    right_row = testFunction.calcAccuracyByRow(w4, major_cellMat)
    rmsd_c = testFunction.calcRMSD(w4, major_cellMat)
    rmsd_f = testFunction.calcRMSD(h_updated, major_F)
    rms_c = testFunction.calcRMSInCell(w4, major_cellMat)
    rms_f = testFunction.calcRMSInCell(h_updated, major_F, Cell=False)

    return iter_nn, dist, average_accuracy, right_row, rmsd_c, rmsd_f, rms_c, rms_f, w4, h_updated
