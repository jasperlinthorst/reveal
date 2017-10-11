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

def graphmumpicker(mums,idx):
    try:
        idxn=idx.nsamples
        if idxn==0:
            return None
        
        if minn==None:
            lminn=idxn
        else:
            lminn=minn
        
        logging.debug("Selecting best out of %d MUMs in graphmumpicker (gap-penalty model: %s, score weight: %d and penalty weight: %d)."%(len(mums),gcmodel,wscore,wpen))
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
                rel=pos-node[0]
                for k in ndata['offsets']:
                    v=ndata['offsets'][k]+rel
                    qlpoint[k]=v
                    qrpoint[k]=v+l

            rlpoint=dict()
            lpenalty=0
            if nleft!=None:
                for k in qlpoint:
                    rlpoint[k]=nleft['offsets'][k]+(nlefttup[1]-nlefttup[0])
                qlp=[]
                for k in sorted(qlpoint):
                    qlp.append(qlpoint[k])
                rlp=[]
                for k in sorted(qlpoint):
                    rlp.append(rlpoint[k])
                lpenalty=utils.gapcost(rlp,qlp,model=gcmodel)

            rrpoint=dict()
            rpenalty=0
            if nright!=None:
                for k in qrpoint:
                    rrpoint[k]=nright['offsets'][k]
                qrp=[]
                for k in sorted(qrpoint):
                    qrp.append(qrpoint[k])
                rrp=[]
                for k in sorted(qrpoint):
                    rrp.append(rrpoint[k])
                rpenalty=utils.gapcost(qrp,rrp,model=gcmodel)

            penalty=lpenalty+rpenalty
            
            score=(wscore*(l*(n**exp)))-(wpen*penalty)
            
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

def printSA(index,maxline=100,start=0,end=200):
    sa=index.SA
    lcp=index.LCP
    t=index.T
    #so=index.SO
    print len(sa), len(lcp)
    assert(len(sa)==len(lcp))
    for s,l in zip(sa[start:end],lcp[start:end]):
        print str(s).zfill(8), str(l).zfill(6), t[s:s+l].ljust(maxline) if lcp<=maxline else t[s:s+maxline].ljust(maxline)#, so[s]
