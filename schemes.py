# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 17:59:26 2015

@author: jasperlinthorst
"""

from intervaltree import IntervalTree, Interval
import networkx as nx
import sys

global minlength, minscore, pcutoff
minlength=20
minscore=0
pcutoff=1e-3
ts=IntervalTree()
G=nx.DiGraph()

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
    try:
        bestmum=None
        trace=False
        #if idx.left!=None and idx.right!=None:
        #    leftoffsets=G.node[idx.left]['offsets']
        #    rightoffsets=G.node[idx.right]['offsets']
        #    if 'CTRI_2.fasta' in leftoffsets and 'CTRI_2.fasta' in rightoffsets:
        #        if leftoffsets['CTRI_2.fasta']==341197 and rightoffsets['CTRI_2.fasta']==341395:
        #            trace=True
        
        for multimum in multimums:
            if trace:
                print multimum
            
            l,n,sp=multimum
            if l<minlength:
                continue
            if n<minn:
                continue
            
            if bestmum!=None:
                if n<bestmum[2]:
                    if trace:
                        print "n < bestmum"
                    continue
                if n==bestmum[2] and l<=bestmum[0]:
                    if trace:
                        print "n equal and l <= bestmum"
                    continue

            if idx.nsamples==len(idx.nodes):
                ds=[start-ts[start].pop()[0] for start in sp]
                ads=sum(ds)/n
                spenalty=sum([abs(p-ads) for p in ds])
            
                de=[ts[start].pop()[1]-(start+l) for start in sp]
                ade=sum(de)/n
                epenalty=sum([abs(p-ade) for p in de])
            
                penalty=min([spenalty,epenalty])
                
                if trace:
                    print "penalty",penalty
                
                score=(l*n)-penalty
                
                if trace:
                    print "score", score
                
            else:
                if trace:
                    "print nsamples != len(nodes)"
                score=minscore #in case of multi aligning graph (havent tried yet) we cant penalize gaps this easily..
            
            if score<minscore:
                if trace:
                    print "score too low"
                continue
            
            if isinstance(sp,tuple):
                sp=list(sp)
            bestmum=(l,idx,n,score,sp,penalty)

            if trace:
                print bestmum
        
        if trace:
            print bestmum
       
        return bestmum
    except:
        print "MULITMUMPICKER ERROR", sys.exc_info()[0]
        return None

def graphmumpicker(mums,idx,penalize=True):
    try:
        bestmum=None
        bestn=2
        bestscore=None
        nleft=None
        nright=None
        trace=False
        
        #if (8306560, 8306677) in set(idx.nodes):#==set([(6742328, 6742408), (2221063, 2221143)]):
        #    trace=True
        #    print idx.left, idx.right, idx.n
        #    print "nodes",idx.nodes
        
        if idx.left!=None:
            nlefttup=idx.left
            nleft=G.node[Interval(nlefttup[0],nlefttup[1])]
        
        if idx.right!=None:
            nrighttup=idx.right
            nright=G.node[Interval(nrighttup[0],nrighttup[1])]
        
        for mum in mums:
            if trace:
                print mum

            l,n,sp=mum
            
            if l<minlength:
                continue
            
            n1=iter(ts[sp[0]]).next()
            n2=iter(ts[sp[1]]).next()
            n1data=G.node[n1]
            n2data=G.node[n2]
            n1samples=set(n1data['offsets'].keys())
            n2samples=set(n2data['offsets'].keys())
            n=len(n1samples)+len(n2samples) #make sure we intersect only the paths that we're following
            
            if n<minn:
                continue
            
            if bestmum!=None:
                if n<bestn:
                    continue
                if n==bestn and l<=bestmum[0]:
                    continue
            
            if penalize:
                if nleft!=None:
                    n1distfromlmapoints=dict()
                    for sample in n1samples:
                        n1distfromlmapoints[sample]=(n1data['offsets'][sample]+(sp[0]-n1[0]))-(nleft['offsets'][sample]+(nlefttup[1]-nlefttup[0]))
                    
                    n2distfromlmapoints=dict()
                    for sample in n2samples:
                        n2distfromlmapoints[sample]=(n2data['offsets'][sample]+(sp[1]-n2[0]))-(nleft['offsets'][sample]+(nlefttup[1]-nlefttup[0]))
                    
                    a=n1distfromlmapoints.values()+n2distfromlmapoints.values()
                    ads=sum(a)/n
                    spenalty=sum([abs(p-ads) for p in a])
                else:
                    spenalty=0
                
                if nright!=None:
                    n1distfromrmapoints=dict()
                    for sample in n1samples:
                        n1distfromrmapoints[sample]=nright['offsets'][sample] - ( n1data['offsets'][sample] + (sp[0]-n1[0]) + l )
                    n2distfromrmapoints=dict()
                    for sample in n2samples:
                        n2distfromrmapoints[sample]=nright['offsets'][sample] - ( n2data['offsets'][sample] + (sp[1]-n2[0]) + l )
                    a=n1distfromrmapoints.values()+n2distfromrmapoints.values()
                    ads=sum(a)/n
                    epenalty=sum([abs(p-ads) for p in a])
                else:
                    epenalty=0
                
                penalty=min([spenalty,epenalty]) #TODO: calculate indel penalty based on distances from graph
            else:
                penalty=0
            
            score=(l*n)-penalty
            
            if trace:
                print l, penalty, score, sp
            
            if score<minscore:
                continue
            
            bestmum=(l,idx,2,score,[sp[0],sp[1]],penalty)
            bestn=n
            if trace:
                print bestmum
            
        return bestmum
    except:
        print "graphmumpicker",n1data['offsets'],len(n1data['offsets']),n2data['offsets'],len(n2data['offsets']),len(nleft['offsets']),len(nright['offsets'])
        print "GRAPHMUMPICKER ERROR", sys.exc_info()[0]
        return None 

def printSA(index,maxline=100,start=0,end=200):
    sa=index.SA
    lcp=index.LCP
    t=index.T
    #so=index.SO
    print len(sa), len(lcp)
    assert(len(sa)==len(lcp))
    for s,l in zip(sa[start:end],lcp[start:end]):
        print str(s).zfill(8), str(l).zfill(6), t[s:s+l].ljust(maxline) if lcp<=maxline else t[s:s+maxline].ljust(maxline)#, so[s]
