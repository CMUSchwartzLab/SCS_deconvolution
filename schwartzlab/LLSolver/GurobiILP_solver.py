from gurobipy import *
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


def updateProportion(B, Ctotal, cells=None, root=0, vType='C', dirichlet=False, dirA=None, beta=0.0):
    tumors = B.shape[0]
    cellsTotal = Ctotal.shape[0]
    copyNum = Ctotal.shape[1]

    if (cells == None):
        cells = cellsTotal

    m = Model('phylo')
    m.setParam("OutputFlag", 0)
    F = m.addVars(tumors, cells, name='F', vtype=vType)
    bDelta = m.addVars(tumors, copyNum, name='bDelta', vtype='C')
    for p in range(tumors):
        for s in range(copyNum):
            expr = B[p, s]
            for k in range(cells):
                expr += -(F[p, k]*Ctotal[k, s])
            m.addConstr(bDelta[p, s]-(expr) >= 0,
        name='residu(%s,%s)' % (p, s))
            m.addConstr(bDelta[p, s]+(expr) >= 0,
        name='residl(%s,%s)' % (p, s))

    for p in range(tumors):
        for k in range(cells):
            m.addConstr(F[p, k] <= 1, name='Fupp(%s,%s)' % (p,k))
            m.addConstr(F[p, k] >= 10 ** (-4),
        name='Flow(%s,%s)' % (p,k))
    m.addConstrs((F.sum(p, '*') == 1 for p in range(tumors)),
    name='Fsum')
    if (dirichlet):
        m.setObjective(bDelta.sum('*', '*'), GRB.MINIMIZE)
        # ptsX = [10 ** (-4), 10 **(-3), 0.01, 0.1, 1]
        ptsX = np.linspace(10 ** (-4), 1, 300)
        for p in range(tumors):
            for k in range(cells):
                ptsY = -beta * (dirA[p, k]-1) * np.log(ptsX)
                m.setPWLObj(F[p, k], ptsX, ptsY)
    else:
        m.setObjective(bDelta.sum('*', '*'), GRB.MINIMIZE)
    m.write('update_proportion.lp')
    m.optimize()

    if (m.status != 2):
        print('Update with no F penalty.')
        print(m.status)
        return updateProportion(B, Ctotal, cells, root, vType, dirichlet=False)

    Fresult = np.zeros([tumors, cells], dtype=np.float)
    for p in range(tumors):
        for k in range(cells):
            Fresult[p, k] = F[p, k].getAttr(GRB.Attr.X)

    objVal = m.objVal
    m.reset()
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
    m.setParam("OutputFlag", 0)

    g = m.addVars(cellsTotal, cellsTotal, cellsTotal, vtype='B', name='g')
    S = m.addVars(cellsTotal, cellsTotal, vtype='B', name='S')

    for t in range(cellsTotal):
        for u in range(cellsTotal):
            if u != t and u != root:
                expr = 0
                for v in range(cellsTotal):
                    expr += (g[u, v, t] - g[v, u, t])
                m.addConstr(expr == 0, name='conserv(%s,%s)' % (u, t))

    for t in range(cellsTotal):
        if t != root:
            m.addConstr(g.sum('*', t, t) == 1, name='sink(%s)' % (t, ))
    m.addConstrs(
        (g.sum(t, '*', t) == 0 for t in range(cellsTotal)), name='leaf')
    m.addConstr(g.sum('*', root, '*') == 0, name='root')
    for t in range(cellsTotal):
        if t != root:
            m.addConstr(g.sum(root, '*', t) == 1, name='source(%s)' % (t, ))

    for t in range(cellsTotal):
        for u in range(cellsTotal):
            for v in range(cellsTotal):
                m.addConstr(
                    S[u, v] - g[u, v, t] >= 0,
                    name='present(%s,%s,%s)' % (u, v, t))

    wDelta = getDistanceMatrix(Ctotal)

    objExpr = S.sum('*', '*')
    for u in range(cellsTotal):
        for v in range(cellsTotal):
            objExpr += wDelta[u, v] * S[u, v]
    objExpr = alpha * objExpr

    m.setObjective(objExpr, GRB.MINIMIZE)
    m.optimize()

    if (m.status != 2):
        return None

    Sresult = np.zeros([cellsTotal, cellsTotal], dtype=np.int8)
    for u in range(cellsTotal):
        for v in range(cellsTotal):
            Sresult[u, v] = S[u, v].getAttr(GRB.Attr.X)

    objVal = m.objVal
    m.reset()
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


def updateCopyNum(B, F, S, CRefer, cells, alpha=0.1, root=0, vType='C', stopVal=None, Cap=False):
    tumors = B.shape[0]
    cellsTotal = cells + CRefer.shape[0]
    copyNum = B.shape[1]
    m = Model('phylo')
    m.setParam("OutputFlag", 0)
    if stopVal != None:
        m.setParam("BestObjStop", stopVal)
    C = m.addVars(cells, copyNum, name='C', vtype=vType)

    bDelta = m.addVars(tumors, copyNum, name='bDelta', vtype='C')
    for p in range(tumors):
        for s in range(copyNum):
            expr = B[p, s]
            for k in range(cells):
                expr += -(F[p, k] * C[k, s])
            m.addConstr(
                bDelta[p, s] - (expr) >= 0, name='resu(%s,%s)' % (p, s))
            m.addConstr(
                bDelta[p, s] + (expr) >= 0, name='resl(%s,%s)' % (p, s))
    m.addConstrs(
        (C[k, s] >= 0 for k in range(cells) for s in range(copyNum)),
        name='pos')

    if (Cap == True):
        m.addConstrs(
            (C[k, s] <= 10 for k in range(cells) for s in range(copyNum)),
            name='cap')

    wDelta = m.addVars(
        cellsTotal, cellsTotal, copyNum, vtype='C', name='wDelta')

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
                m.addConstr(
                    wDelta[u, v, i] - Cui + Cvi >= 0,
                    name='delu(%s,%s,%s)' % (u, v, i))
                m.addConstr(
                    wDelta[u, v, i] - Cvi + Cui >= 0,
                    name='dell(%s,%s,%s)' % (u, v, i))

    objExpr = bDelta.sum('*', '*')
    for u in range(cellsTotal):
        for v in range(cellsTotal):
            objExpr += alpha * wDelta.sum(u, v, '*') * S[u, v]

    m.setObjective(objExpr, GRB.MINIMIZE)
    m.optimize()

    # print('Total Obj:', m.objVal)
    # 2 is OPTIMAL
    # print('Status:', m.status)

    if (m.status != 2 and m.status != 15):
        return None

    Cresult = np.zeros([cells, copyNum], dtype=np.float)
    for k in range(cells):
        for s in range(copyNum):
            Cresult[k, s] = C[k, s].getAttr(GRB.Attr.X)

    objVal = m.objVal
    m.reset()
    return [Cresult, objVal]


def makeSurePath(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)
    return
