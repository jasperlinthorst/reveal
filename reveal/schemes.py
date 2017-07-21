# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 17:59:26 2015

@author: jasperlinthorst
"""

from intervaltree import IntervalTree, Interval
import networkx as nx
import sys
import math
import logging

def multimumpicker(multimums,idx):

    idxn=idx.nsamples
    if idxn==0:
        return None
    
    if minn==None:
        lminn=idxn
    else:
        lminn=minn
    
    try:
        logging.debug("Selecting best out of %d MUMs in multimumpicker (wpen=%s, wscore=%s, exp=%s, minn=%s ,samples=%s)."%(len(multimums),wpen,wscore,exp,lminn,idxn))
        bestmum=None
        
        for multimum in multimums:
            l,n,sp=multimum
            if l<minlength:
                continue
            if n<lminn:
                continue
            if bestmum!=None:
                if n<bestmum[2]:
                    continue
                if n==bestmum[2] and l<=bestmum[0]:
                    continue
            
            if idx.depth==0: #no other anchors yet, dont penalize first anchor
                score=int(wscore*(l*(n**exp)))
                penalty=0
            else:
                
                if wpen>0:
                    ds=[start-ts[start].pop()[0] for start in sp]
                    ads=sum(ds)/n
                    de=[ts[start].pop()[1]-(start+l) for start in sp]
                    ade=sum(de)/n
                    spenalty=sum([abs(p-ads) for p in ds])
                    epenalty=sum([abs(p-ade) for p in de])
                    penalty=min([spenalty,epenalty])
                    penalty=int((wpen*penalty)**(.5))
                else:
                    penalty=0
                
                score=int(wscore*(l*(n**exp)))
                
                if minscore!=None:
                    if penalty>score:
                        score=minscore
                    else:
                        score=score-penalty
            
            if minscore!=None:
                if score<minscore:
                    continue
            
            if isinstance(sp,tuple):
                sp=list(sp)
            
            if bestmum!=None:
                if score>bestmum[3]:
                    bestmum=(l,idx,n,score,sp,penalty)
            else:
                bestmum=(l,idx,n,score,sp,penalty)
        
        if bestmum!=None:
            if bestmum[3] > sys.maxint:
                bestmum=(bestmum[0],bestmum[1],bestmum[2],sys.maxint,bestmum[4],bestmum[5])
            if bestmum[2]!=idxn:
                logging.debug("Branching of group of %d samples out of %s based on match of length %d."%(bestmum[2],idxn,bestmum[0]))
        
        logging.debug("Selected %s"%str(bestmum))
        
        return bestmum
    except Exception as e:
        print "MULTIMUMPICKER ERROR", e, sys.exc_info()[0]
        return None

def graphmumpicker(mums,idx,penalize=True):
    idxn=idx.nsamples
    if idxn==0:
        return None
     
    if minn==None:
        lminn=idxn
    else:
        lminn=minn
    
    try:
        logging.debug("Selecting best out of %d MUMs in graphmumpicker (score weight: %d and penalty weight: %d)."%(len(mums),wscore,wpen))
        if not penalize:
            logging.debug("Not penalizing MUMs.")
        
        bestmum=None
        bestn=2
        bestscore=None
        nleft=None
        nright=None
        trace=False
        
        if idx.left!=None:
            nlefttup=idx.left
            nleft=G.node[Interval(nlefttup[0],nlefttup[1])]

        if idx.right!=None:
            nrighttup=idx.right
            nright=G.node[Interval(nrighttup[0],nrighttup[1])]

        if idx.left!=None and idx.right!=None:
            leftoffsets=G.node[idx.left]['offsets']
            rightoffsets=G.node[idx.right]['offsets']
        
        for mum in mums:
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
            
            if n<lminn:
                continue
            
            if bestmum!=None:
                if n<bestn:
                    continue
                if n==bestn and l<=bestmum[0]:
                    continue
            
            if penalize:
                spenalty=None
                epenalty=None
                
                if nleft!=None:
                    n1distfromlmapoints=dict()
                    for sample in n1samples:
                        n1distfromlmapoints[sample]=(n1data['offsets'][sample]+(sp[0]-n1[0]))-(nleft['offsets'][sample]+(nlefttup[1]-nlefttup[0]))
                    n2distfromlmapoints=dict()
                    for sample in n2samples:
                        n2distfromlmapoints[sample]=(n2data['offsets'][sample]+(sp[1]-n2[0]))-(nleft['offsets'][sample]+(nlefttup[1]-nlefttup[0]))
                    spenalty=mindist(n1distfromlmapoints.values(),n2distfromlmapoints.values())
                    
                if nright!=None:
                    n1distfromrmapoints=dict()
                    for sample in n1samples:
                        n1distfromrmapoints[sample]=nright['offsets'][sample] - ( n1data['offsets'][sample] + (sp[0]-n1[0]) + l )
                    n2distfromrmapoints=dict()
                    for sample in n2samples:
                        n2distfromrmapoints[sample]=nright['offsets'][sample] - ( n2data['offsets'][sample] + (sp[1]-n2[0]) + l )
                    epenalty=mindist(n1distfromrmapoints.values(),n2distfromrmapoints.values())

                if epenalty==None and spenalty==None:
                    penalty=0
                elif epenalty==None:
                    penalty=spenalty
                elif spenalty==None:
                    penalty=epenalty
                else:
                    penalty=min([spenalty,epenalty]) #TODO: calculate indel penalty based on distances from graph
            else:
                penalty=0
            
            score=(wscore*(l*(n**exp)))-((wpen*penalty)**(.5))
            #score=(wscore*(l*n))-(wpen*penalty)
            #score=int(math.log(l-(wpen*penalty))*n)
            
            if trace:
                print l, penalty, score
            
            if minscore!=None:
                if score<minscore:
                    continue
            
            bestmum=(l,idx,2,score,[sp[0],sp[1]],penalty)
            
            bestn=n
            
            if trace:
                print bestmum
    
        logging.debug("Selected: %s\n"%str(bestmum))
        return bestmum
    except Exception as e:
        #print "graphmumpicker",n1data['offsets'],len(n1data['offsets']),n2data['offsets'],len(n2data['offsets']),len(nleft['offsets']),len(nright['offsets'])
        print "GRAPHMUMPICKER ERROR", e, sys.exc_info()[0]
        return None 

def mindist(x,y):
    x=sorted(list(set(x)))
    y=sorted(list(set(y)))
    i=0
    j=0
    mind=abs(x[i]-y[j])
    while True:
        if i<len(x)-1:
            xdif=abs(x[i+1] - y[j])
        else:
            xdif=abs(x[i] - y[j])

        if j<len(y)-1:
            ydif=abs(x[i] - y[j+1])
        else:
            ydif=abs(x[i] - y[j])
        
        if xdif < ydif:
            if xdif<mind:
                mind=xdif
            i+=1
            if i==len(x):
                break
        else:
            if ydif<mind:
                mind=ydif
            j+=1
            if j==len(y):
                break
        
        if mind==0:
            break
    return mind

def printSA(index,maxline=100,start=0,end=200):
    sa=index.SA
    lcp=index.LCP
    t=index.T
    #so=index.SO
    print len(sa), len(lcp)
    assert(len(sa)==len(lcp))
    for s,l in zip(sa[start:end],lcp[start:end]):
        print str(s).zfill(8), str(l).zfill(6), t[s:s+l].ljust(maxline) if lcp<=maxline else t[s:s+maxline].ljust(maxline)#, so[s]
