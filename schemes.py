# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 17:59:26 2015

@author: jasperlinthorst
"""

from intervaltree import IntervalTree

global minlength, minscore, pcutoff
minlength=20
minscore=0
pcutoff=1e-3
ts=IntervalTree()
interval2sampleid=dict()

#take multi-mum that is observed in all samples until it drops below maxl
#then, check calculate significance for the multimum if it's not significant
#take the next multi-mum that is not observed in all samples but is significant
def mumpicker(multimums,idx):
    maxscore=0
    maxl=0
    nomum=None
    #nomum=(0,idx,0,[])
    bestmum=nomum
    for multimum in multimums:
        l,n,sp=multimum
        if n==idx.nsamples:
            if l>maxl:
                maxl=l
                bestmum=(l,idx,n,sp)
    if bestmum==None:
        return None
    p=(1-((1-(0.25**bestmum[0]))**(((idx.n/float(idx.nsamples))**2))))**(idx.nsamples-1)
    if p>pcutoff:
        bestmum=(0,idx,0,[])
        multimums.sort(key=lambda l: l[0]*l[1],reverse=True)
        for multimum in multimums:
            l,n,sp=multimum
            p=(1-((1-(0.25**l))**(((idx.n/idx.nsamples)**2))))**(idx.nsamples-1)
            if (p<=pcutoff):
                return (l,idx,n,sp)
        #not a single significant multi-mum
        return nomum
    else:
        #one multimum for all samples index
        return bestmum

#take the largest multimum that occurs in all samples, until it drops below threshold
#then take largest multimum in less samples that is above thresholds
#TODO: calculate score for multi-mums!
def mumpicker2(multimums,idx):
    bestmum_by_n={}
    bestmum=None    
    for multimum in multimums:
        l,n,sp=multimum
        if n in bestmum_by_n:
            if l>bestmum_by_n[n][0]:
                bestmum_by_n[n]=(l,idx,n,0,sp)
        else:
            bestmum_by_n[n]=(l,idx,n,0,sp)
    
    for key in sorted(bestmum_by_n,reverse=True):
        if bestmum_by_n[key][0]>=minlength:
            bestmum=bestmum_by_n[key]
            break
        #else:
        #    cluster(multimums,idx)
    return bestmum

#take the best scoring multimum that occurs in all samples, until it drops below threshold
#then take largest multimum in less samples that is above thresholds
#TODO: calculate score for multi-mums!
def mumpicker3(multimums,idx):
    bestmum_by_n={}
    bestmum=None
    for multimum in multimums:
        l,n,sp=multimum
        ds=[start-ts[start].pop()[0] for start in sp]
        ads=sum(ds)/n
        spenalty=sum([abs(p-ads) for p in ds])
        de=[ts[start].pop()[1]-(start+l) for start in sp]
        ade=sum(de)/n
        epenalty=sum([abs(p-ade) for p in de])
        penalty=min([spenalty,epenalty])
        score=(l*(n**2))-penalty

        if n in bestmum_by_n:
            if score>bestmum_by_n[n][3]:
                bestmum_by_n[n]=(l,idx,n,score,sp)
                #print "better",l, n, score, penalty, spenalty, epenalty, sp, ds, de
        else:
            bestmum_by_n[n]=(l,idx,n,score,sp)
            #print "initial",l, n, score, penalty, spenalty, epenalty, sp, ds, de
    
    for key in sorted(bestmum_by_n,reverse=True):
        if bestmum_by_n[key][3]>=minscore:
            bestmum=bestmum_by_n[key]
            break
    
    #print "bestmum",bestmum
    #print index.nsamples, bestmum_by_n, bestmum
    return bestmum

def multimumpicker(multimums,idx):
    bestmum=None
    for multimum in multimums:
        l,n,sp=multimum
        if l<minlength:
            continue
        
        if bestmum!=None:
            if n<bestmum[2]:
                continue
            if n==bestmum[2] and l<=bestmum[0]:
                continue
                
        ds=[start-ts[start].pop()[0] for start in sp]
        ads=sum(ds)/n
        spenalty=sum([abs(p-ads) for p in ds])
        de=[ts[start].pop()[1]-(start+l) for start in sp]
        ade=sum(de)/n
        epenalty=sum([abs(p-ade) for p in de])
        penalty=min([spenalty,epenalty])
        score=(l*n)-penalty
        
        if score<minscore:
            continue
        
        bestmum=(l,idx,n,score,sp)
    
    return bestmum


def cluster(multimums,idx):
    T=idx.T
    for node in idx.nodes:
        print T[node[0]:node[1]]
    print multimums
    printSA(idx)
    
    S=[[0 for x in range(len(ts))] for x in range(len(ts))]
    for multimum in multimums:
        l,n,sp=multimum
        for i in range(len(sp)):
            for j in range(len(sp)):
                sid1=interval2sampleid[ ts[sp[i]].pop() ]
                sid2=interval2sampleid[ ts[sp[j]].pop() ]
                S[sid1][sid2]+=l
    from sklearn.decomposition import PCA
    from matplotlib import pyplot as plt
    import numpy as np
    S=np.array(S)
    pca = PCA(n_components=2)
    X=pca.fit(S).transform(S)
    for row in X:
        plt.plot(row[0],row[1],'r.')
    plt.show()    
    
def printSA(index,maxline=100,start=0,end=200):
    sa=index.SA
    lcp=index.LCP
    t=index.T
    #so=index.SO
    print len(sa), len(lcp)
    assert(len(sa)==len(lcp))
    for s,l in zip(sa[start:end],lcp[start:end]):
        print str(s).zfill(8), str(l).zfill(6), t[s:s+l].ljust(maxline) if lcp<=maxline else t[s:s+maxline].ljust(maxline)#, so[s]
