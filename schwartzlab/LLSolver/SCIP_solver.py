from pyscipopt import Model
from pyscipopt import quicksum
import numpy as np
import csv
import sys
import time
import random
import pickle
from multiprocessing import Process
from multiprocessing import Pool
import scipy.io as sio
import testFunction as test
import os

'''
Update Fraction Matrix.
'''
def updateProportion(B, Ctotal, cells=None, root=0, vType='C', dirA=None, beta=0.0):
    tumors = B.shape[0]
    cellsTotal = Ctotal.shape[0]
    copyNum = Ctotal.shape[1]

    if (cells == None):
        cells = cellsTotal

    m = Model('phylo')
    F = {}
    for i in range(tumors):
        for j in range(cells):
            F[i, j] = m.addVar(vtype=vType, name="F(%s,%s)" % (i, j))

    bDelta = {}
    for i in range(tumors):
        for j in range(copyNum):
            bDelta[i, j] = m.addVar(vtype='C', name='bDelta(%s,%s)' % (i, j))

    for p in range(tumors):
        for s in range(copyNum):
            expr = B[p, s]
            for k in range(cells):
                expr += -(F[p, k]*Ctotal[k, s])
            m.addCons(bDelta[p, s]-(expr) >= 0,
        name='residu(%s,%s)' % (p, s))
            m.addCons(bDelta[p, s]+(expr) >= 0,
        name='residl(%s,%s)' % (p, s))

    for p in range(tumors):
        for k in range(cells):
            m.addCons(F[p, k] <= 1, name='Fupp(%s,%s)' % (p,k))
            m.addCons(F[p, k] >= 10 ** (-4),
        name='Flow(%s,%s)' % (p,k))

    for p in range(tumors):
        m.addCons(quicksum(F[p, k] for k in range(cells)) == 1,
        'Fsum(%s)' % (p, ))

    m.setObjective(quicksum(bDelta[i, j] for i in range(
        tumors) for j in range(copyNum)), "minimize")
    m.optimize()

    #print('Total Obj:', m.getObjVal())

    Fresult = np.zeros([tumors, cells], dtype=np.float)
    for p in range(tumors):
        for k in range(cells):
            Fresult[p, k] = m.getVal(F[p, k])

    objVal = m.getObjVal()
    return [Fresult, objVal]


def getDistanceMatrix(Ctotal):
    cellsTotal = Ctotal.shape[0]
    distance = np.zeros([cellsTotal, cellsTotal])

    for i in range(cellsTotal):
        for j in range(cellsTotal):
            distance[i, j] = np.sum(abs(Ctotal[i, :] - Ctotal[j, :]))

    return distance


'''
Update topological tree structure of phylogenetic tree.
'''
def updateTree(B, Ctotal, cells=None, alpha=0.1, root=0):
    tumors = B.shape[0]
    cellsTotal = Ctotal.shape[0]
    copyNum = Ctotal.shape[1]

    if (cells == None):
        cells = cellsTotal

    m = Model('phylo')

    g, S = {}, {}
    for i in range(cellsTotal):
        for j in range(cellsTotal):
            S[i, j] = m.addVar(vtype='B', name='S(%s,%s)' % (i, j))
            for k in range(cellsTotal):
                g[i, j, k] = m.addVar(
                    vtype='B', name='g(%s,%s,%s)' % (i, j, k))

    for t in range(cellsTotal):
        for u in range(cellsTotal):
            if u != t and u != root:
                expr = 0
                for v in range(cellsTotal):
                    expr += (g[u, v, t] - g[v, u, t])
                m.addCons(expr == 0, name='conserv(%s,%s)' % (u, t))

    for t in range(cellsTotal):
        if t != root:
            m.addCons(
                quicksum(g[i, t, t] for i in range(cellsTotal)) == 1,
                name='sink(%s)' % (t, ))
        m.addCons(
            quicksum(g[t, i, t] for i in range(cellsTotal)) == 0,
            name='leaf(%s)' % (t, ))

    m.addCons(
        quicksum(g[i, root, j] for i in range(cellsTotal)
                 for j in range(cellsTotal)) == 0,
        name='root')
    for t in range(cellsTotal):
        if t != root:
            m.addCons(
                quicksum(g[root, i, t] for i in range(cellsTotal)) == 1,
                name='source(%s)' % (t, ))

    for t in range(cellsTotal):
        for u in range(cellsTotal):
            for v in range(cellsTotal):
                m.addCons(
                    S[u, v] - g[u, v, t] >= 0,
                    name='present(%s,%s,%s)' % (u, v, t))

    wDelta = getDistanceMatrix(Ctotal)

    objExpr = quicksum(S[i, j] for i in range(cellsTotal)
                       for j in range(cellsTotal))
    for u in range(cellsTotal):
        for v in range(cellsTotal):
            objExpr += wDelta[u, v] * S[u, v]

    m.setObjective(alpha * objExpr, "minimize")
    m.optimize()

    #print('Total Obj:', m.getObjVal())

    Sresult = np.zeros([cellsTotal, cellsTotal], dtype=np.int8)
    for u in range(cellsTotal):
        for v in range(cellsTotal):
            Sresult[u, v] = m.getVal(S[u, v])

    objVal = m.getObjVal()
    return [Sresult, objVal]


def calcObjVal(B, F, C, metric='L1'):
    diff = B - np.matmul(F, C)
    if metric == 'L1':
        return(np.sum(abs(diff)))
    elif metric == 'L2':
        return(np.sum(np.multiply(diff, diff)))


'''
Update CNV in single-cell matrix.
'''
def updateCopyNum(B, F, S, CRefer, cells, alpha=0.1, root=0, vType='C', Cap=False):
    tumors = B.shape[0]
    cellsTotal = cells + CRefer.shape[0]
    copyNum = B.shape[1]
    m = Model('phylo')

    C = {}
    for i in range(cells):
        for j in range(copyNum):
            C[i, j] = m.addVar(vtype=vType, name='C(%s,%s)' % (i, j))

    bDelta = {}
    for i in range(tumors):
        for j in range(copyNum):
            bDelta[i, j] = m.addVar(vtype='C', name='bDelta(%s,%s)' % (i, j))

    for p in range(tumors):
        for s in range(copyNum):
            expr = B[p, s]
            for k in range(cells):
                expr += -(F[p, k] * C[k, s])
            m.addCons(bDelta[p, s] - (expr) >= 0, name='resu(%s,%s)' % (p, s))
            m.addCons(bDelta[p, s] + (expr) >= 0, name='resl(%s,%s)' % (p, s))
    for k in range(cells):
        for s in range(copyNum):
            m.addCons(C[k, s] >= 0, name='pos(%s,%s)' % (k, s))

    if (Cap == True):
        for k in range(cells):
            for s in range(copyNum):
                m.addCons(C[k, s] <= 10, name='cap(%s,%s)' % (k, s))

    wDelta = {}
    for i in range(cellsTotal):
        for j in range(cellsTotal):
            for k in range(copyNum):
                wDelta[i, j, k] = m.addVar(
                    vtype='C', name='wDelta(%s,%s,%s)' % (i, j, k))

    for u in range(cellsTotal):
        for v in range(cellsTotal):
            for i in range(copyNum):
                if (u < cells):
                    Cui = C[u, i]
                else:
                    Cui = CRefer[u-cells, i]
                if (v < cells):
                    Cvi = C[v, i]
                else:
                    Cvi = CRefer[v - cells, i]
                m.addCons(
                    wDelta[u, v, i] - Cui + Cvi >= 0,
                    name='delu(%s,%s,%s)' % (u, v, i))
                m.addCons(
                    wDelta[u, v, i] - Cvi + Cui >= 0,
                    name='dell(%s,%s,%s)' % (u, v, i))

    objExpr = quicksum(
        bDelta[i, j] for i in range(tumors) for j in range(copyNum))
    for u in range(cellsTotal):
        for v in range(cellsTotal):
            objExpr += alpha * quicksum(
                wDelta[u, v, i] for i in range(copyNum)) * S[u, v]

    m.setObjective(objExpr, 'minimize')
    m.optimize()

    #print('Total Obj:', m.getObjVal())
    # 2 is OPTIMAL

    Cresult = np.zeros([cells, copyNum], dtype=np.float)
    for k in range(cells):
        for s in range(copyNum):
            Cresult[k, s] = m.getVal(C[k, s])

    objVal = m.getObjVal()
    return [Cresult, objVal]


def makeSurePath(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)
    return
