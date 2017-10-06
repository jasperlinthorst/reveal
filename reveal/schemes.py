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
import utils

# def multimumpicker(multimums,idx):

#     idxn=idx.nsamples
#     if idxn==0:
#         return None
    
#     if minn==None:
#         lminn=idxn
#     else:
#         lminn=minn
    
#     try:
#         logging.debug("Selecting best out of %d MUMs in multimumpicker (wpen=%s, wscore=%s, exp=%s, minn=%s ,samples=%s)."%(len(multimums),wpen,wscore,exp,lminn,idxn))
#         bestmum=None
        
#         for multimum in multimums:
#             l,n,sp=multimum
#             if l<minlength:
#                 continue
#             if n<lminn:
#                 continue
#             if bestmum!=None:
#                 if n<bestmum[2]:
#                     continue
#                 if n==bestmum[2] and l<=bestmum[0]:
#                     continue
            
#             if idx.depth==0: #no other anchors yet, dont penalize first anchor
#                 if exp==None:
#                     score=int(wscore*(l*(n**n)))
#                 else:
#                     score=int(wscore*(l*(n**exp)))
#                 penalty=0
#             else:
                
#                 if wpen>0:
#                     ds=[start-ts[start].pop()[0] for start in sp]
#                     ads=sum(ds)/n
#                     de=[ts[start].pop()[1]-(start+l) for start in sp]
#                     ade=sum(de)/n
#                     spenalty=sum([abs(p-ads) for p in ds])
#                     epenalty=sum([abs(p-ade) for p in de])
#                     penalty=min([spenalty,epenalty])
#                     penalty=int((wpen*penalty)**(.5))
#                 else:
#                     penalty=0
                
#                 if exp==None:
#                     score=int(wscore*(l*(n**n)))
#                 else:
#                     score=int(wscore*(l*(n**exp)))
                
#                 if minscore!=None:
#                     if penalty>score:
#                         score=minscore
#                     else:
#                         score=score-penalty
            
#             if minscore!=None:
#                 if score<minscore:
#                     continue
            
#             if isinstance(sp,tuple):
#                 sp=list(sp)
            
#             if bestmum!=None:
#                 if score>bestmum[3]:
#                     bestmum=(l,idx,n,score,sp,penalty)
#             else:
#                 bestmum=(l,idx,n,score,sp,penalty)
        
#         if bestmum!=None:
#             if bestmum[3] > sys.maxint:
#                 bestmum=(bestmum[0],bestmum[1],bestmum[2],sys.maxint,bestmum[4],bestmum[5])
#             if bestmum[2]!=idxn:
#                 logging.debug("Branching of group of %d samples out of %s based on match of length %d."%(bestmum[2],idxn,bestmum[0]))
        
#         logging.debug("Selected %s"%str(bestmum))
        
#         return bestmum
#     except Exception as e:
#         print "MULTIMUMPICKER ERROR", e, sys.exc_info()[0]
#         return None

def graphmumpicker(mums,idx):
    try:
        idxn=idx.nsamples
        if idxn==0:
            return None
        
        if minn==None:
            lminn=idxn
        else:
            lminn=minn
        
        logging.debug("Selecting best out of %d MUMs in graphmumpicker (score weight: %d and penalty weight: %d)."%(len(mums),wscore,wpen))
        logging.debug("Selection criteria:")
        if minscore!=None:
            logging.debug("minscore=%d"%minscore)
        logging.debug("minn=%d"%lminn)
        logging.debug("minlength=%d"%minlength)
        
        bestmum=None
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
            l,mmn,sp=mum
            
            if l<minlength:
                continue
            
            if mmn<lminn:
                continue

            n=0
            qlpoint=dict()
            qrpoint=dict()
            for pos in sp:
                node=iter(ts[pos]).next()
                ndata=G.node[node]
                nsamples=set([o for o in ndata['offsets'].keys() if not G.graph['id2sample'][o].startswith("*")])
                n+=len(nsamples)
                rel=node[0]-pos
                for k in ndata['offsets']:
                    v=ndata['offsets'][k]+rel
                    qlpoint[k]=v
                    qrpoint[k]=v+l

            rlpoint=dict()
            lpenalty=0
            if nleft!=None:
                for k in qlpoint:
                    rlpoint[k]=nleft['offsets'][k]+(nlefttup[1]-nlefttup[0])
                qp=[]
                for k in sorted(qlpoint):
                    qp.append(qlpoint[k])
                rp=[]
                for k in sorted(qlpoint):
                    rp.append(rlpoint[k])
                lpenalty=utils.sumofpairs(rp,qp)

            rrpoint=dict()
            rpenalty=0
            if nright!=None:
                for k in qrpoint:
                    rrpoint[k]=nright['offsets'][k]
                qp=[]
                for k in sorted(qrpoint):
                    qp.append(qrpoint[k])
                rp=[]
                for k in sorted(qrpoint):
                    rp.append(rrpoint[k])
                rpenalty=utils.sumofpairs(qp,rp)

            penalty=lpenalty+rpenalty
            
            score=(wscore*(l*(n**n)))-((wpen*penalty)**(.5))
            
            if minscore!=None:
                if score<minscore:
                    continue
            
            if bestmum==None:
                bestmum=(l,idx,mmn,score,sp,penalty)
            else:
                if bestmum[3]<score:
                    bestmum=(l,idx,mmn,score,sp,penalty)
        
        if bestmum!=None:
            if bestmum[3] > sys.maxint: #cap score, so we dont have problems with overflows in c
                bestmum=(bestmum[0],bestmum[1],bestmum[2],sys.maxint,bestmum[4],bestmum[5])
        
        logging.debug("Selected: %s\n"%str(bestmum))
        return bestmum
    except Exception as e:
        print "GRAPHMUMPICKER ERROR", e, sys.exc_info()[0]
        return None 

# def mindist(x,y):
#     x=sorted(list(set(x)))
#     y=sorted(list(set(y)))
#     i=0
#     j=0
#     mind=abs(x[i]-y[j])
#     while True:
#         if i<len(x)-1:
#             xdif=abs(x[i+1] - y[j])
#         else:
#             xdif=abs(x[i] - y[j])

#         if j<len(y)-1:
#             ydif=abs(x[i] - y[j+1])
#         else:
#             ydif=abs(x[i] - y[j])
        
#         if xdif < ydif:
#             if xdif<mind:
#                 mind=xdif
#             i+=1
#             if i==len(x):
#                 break
#         else:
#             if ydif<mind:
#                 mind=ydif
#             j+=1
#             if j==len(y):
#                 break
        
#         if mind==0:
#             break
#     return mind

def printSA(index,maxline=100,start=0,end=200):
    sa=index.SA
    lcp=index.LCP
    t=index.T
    #so=index.SO
    print len(sa), len(lcp)
    assert(len(sa)==len(lcp))
    for s,l in zip(sa[start:end],lcp[start:end]):
        print str(s).zfill(8), str(l).zfill(6), t[s:s+l].ljust(maxline) if lcp<=maxline else t[s:s+maxline].ljust(maxline)#, so[s]
