#!/usr/bin/env python

import pandas as pd
import numpy as np

pctThres = 0.1

def load_lr_pairs(fname):
    lrPairs = pd.read_csv(fname)
    lrPairs.columns = ['L','R']
    return lrPairs

def load_lr_list(fname):
    lrList = []
    f = open(fname,'r')
    for line in f:
        lrList.append(line.strip())
    f.close()
    return lrList

def load_matrix(fname):
    return pd.read_csv(fname,index_col=0)

def load_ident(fname):
    identList = pd.read_csv(fname,index_col=0)
    identList.columns = ['ident']
    return identList

def gen_sig_pairs(scDataFrame,identList,lrPairs,lrList,iterNum):
    """
    format of pctMat:
        clust1    clust2
    gene1   0   0
    gene2   0   0
    format of eachPair:
        [identL,identR,L,R]
    """
    scDataFrame = scDataFrame.join(identList)
    identNum = scDataFrame['ident'].unique()
    pctMat = pd.DataFrame(0,index=lrList,columns=identNum)
    for ident in identNum:
        logiMat = scDataFrame.loc[scDataFrame['ident']==ident,lrList]>0
        pctMat.loc[:,ident] = logiMat.mean()
    expPairs = []
    for eachIndex,eachRow in lrPairs.iterrows():
        for identL in identNum:
            for identR in identNum:
                if (pctMat.loc[eachRow['L'],identL] > pctThres) & (pctMat.loc[eachRow['R'],identR] > pctThres):
                    expPairs.append([identL,identR,eachRow['L'],eachRow['R']])
    print(len(expPairs))
    expRealMat = pd.DataFrame(0,index=lrList,columns=identNum)
    for ident in identNum:
        expRealMat.loc[:,ident] = scDataFrame.loc[scDataFrame['ident']==ident,lrList].mean()
    randomMeanMat = []
    for i in range(iterNum):
        randomMeanList = []
        scDataFrame['ident'] = identList.sample(frac=1,replace=True)
        randomIdentList = identList.sample(frac=1,replace=True)
        randomIdentList.reset_index(inplace=True,drop=True)
        randomIdentList.index = scDataFrame.index
        scDataFrame.update(randomIdentList)
        expRandomMat = pd.DataFrame(0,index=lrList,columns=identNum)
        for ident in identNum:
            expRandomMat.loc[:,ident] = scDataFrame.loc[scDataFrame['ident']==ident,lrList].mean()
        j = 0
        for eachExpPair in expPairs:
            lMean = expRandomMat.loc[eachExpPair[2],eachExpPair[0]]
            rMean = expRandomMat.loc[eachExpPair[3],eachExpPair[1]]
            randomMeanList.append(lMean * rMean)
            j += 1
        randomMeanMat.append(randomMeanList)
        print(i)
    randomMeanDF = pd.DataFrame(randomMeanMat)
    sigPairs = []
    for i in range(len(expPairs)):
        eachExpPair = expPairs[i]
        realExpMult = expRealMat.loc[eachExpPair[2],eachExpPair[0]]*expRealMat.loc[eachExpPair[3],eachExpPair[1]]
        logiMat = randomMeanDF.loc[:,i] > realExpMult
        pval = logiMat.mean()
        FC = realExpMult / randomMeanDF.loc[:,i].mean()
        if pval < 0.01:
            eachExpPair[0] = eachExpPair[0][1:]
            eachExpPair[1] = eachExpPair[1][1:]
            sigPairs.append(eachExpPair+[FC,realExpMult,pval])
    sigPairsDF = pd.DataFrame(sigPairs)
    return sigPairsDF




lrList = load_lr_list('lungmc.data.lrlist.txt')
lrPairs = load_lr_pairs('LR.pairs.core_lungmc.csv')

scDataFrame = load_matrix('lungmcT0.data.csv')
identList = load_ident('lungmcT0.data.ident2.txt')
sigPairsDF = gen_sig_pairs(scDataFrame,identList,lrPairs,lrList,10000)
sigPairsDF.to_csv('sigPairsDF.T0.csv')





