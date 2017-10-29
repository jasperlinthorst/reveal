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
import traceback

def chain(mums,left,right,gcmodel="sumofpairs"):
    if len(mums)==0:
        return []
    logging.debug("Number of anchors before chaining: %d",len(mums))
    
    #use one coordinate system for sorting
    ref=mums[0][2].keys()[0]
    logging.trace("Ref is %s"%ref)
    mums.append(right)
    mums.sort(key=lambda mum: mum[2][ref]) #sort by reference dimension

    sp2mum=dict()
    for mum in mums:
        sp2mum[mum[2][ref]]=mum

    #sps=[mum[2] for mum in mums]
    #mumsptree=utils.kdtree(sps,2)
    #eps=[tuple([sp+mum[0] for sp in mum[2]]) for mum in mums]
    
    #eps=[]
    #for mum in mums:
    #    t=[]
    #    for sp in mum[2]:
    #        t.append(sp+mum[0])
    #    eps.append(tuple(t))
    #mumeptree=utils.kdtree(eps,2)

    minscore=-1*utils.gapcost([left[2][k] for k in right[2]],[right[2][k] for k in right[2]])
    logging.trace("Initial cost is: %d"%minscore)

    start=left[2][ref]
    end=right[2][ref]

    link=dict()
    score=dict({left[2][ref]:0})
    
    active=[left]
    processed=[]
    
    for mum in mums:
        #active=[ep2mum[ep] for ep in utils.range_search(mumeptree,(0,0),[sp-1 for sp in mum[2]])].sort(key=lambda x: score[x], reverse=True)
        remove=[]
        for pmum in processed:
            pendpoint1=pmum[2][ref]+pmum[0]
            if pendpoint1<=mum[2][ref]:
                active.append(pmum)
                remove.append(pmum)

        for r in remove:
            processed.remove(r)

        #subset active points
        #aeps=[ tuple([sp+mum[0] for sp in mum[2]]) for mum in active]
        #amumeptree=utils.kdtree(aeps,2)
        #active=[aep2mum[ep] for ep in utils.range_search(amumeptree,(0,0),[sp-1 for sp in mum[2]])]#.sort(key=lambda x: score[x], reverse=True) #subset active set!
        
        active.sort(key=lambda x: score[x[2][ref]], reverse=True) #sort active by score decreasing, kind of priority queue
        
        # subactive=[]
        # for amum in active:
        #     for psp,sp in zip(amum[2],mum[2]):
        #         if psp+amum[0]>=sp:
        #             break
        #     else: #pmum does not overlap in any dimension
        #         subactive.append(amum)
        #subactive.sort(key=lambda x: score[x], reverse=True) #sort active by score decreasing, kind of priority queue
        w=None
        for amum in active:
            s=score[amum[2][ref]]+(wscore*(mum[0]*((mum[1]*(mum[1]-1))/2) ))

            if w!=None:
                if w > s: #as input is sorted by score
                    break
            
            overlap=False
            for crd in mum[2]:
                if amum[2][crd]+amum[0]>mum[2][crd]:
                    overlap=True
            if overlap:
                continue

            penalty=utils.gapcost([amum[2][k] for k in mum[2]],[mum[2][k] for k in mum[2]],model=gcmodel)

            assert(penalty>=0)

            tmpw=score[amum[2][ref]]+(wscore*(mum[0]*((mum[1]*(mum[1]-1))/2)))-(wpen*penalty)

            if tmpw>w or w==None:
                w=tmpw
                best=amum

        link[mum[2][ref]]=best[2][ref]
        score[mum[2][ref]]=w

        processed.append(mum)

    assert(score[end]>=minscore)

    #backtrack
    path=[]
    while end!=start:
        path.append(sp2mum[end])
        end=link[end]

    return path[1:]

#determine a subset of genomes for which (length * n) is largest
def segment(mums):
    d=dict()
    for mum in mums:
        k=tuple(sorted(mum[2].keys()))
        if k in d:
            d[k].append(mum)
        else:
            d[k]=[mum]
    best=0
    for part in d:
        z=sum([m[0] for m in d[part]])*len(part)
        if z>best:
            best=part
    return d[best]

def lookup(mum):
    l,mmn,spd=mum
    sp=spd.values()
    n=0
    qlpoint=dict()
    qrpoint=dict()
    for pos in sp:
        node=iter(ts[pos]).next()
        ndata=G.node[node]
        nsamples=set([o for o in ndata['offsets'].keys() if not G.graph['id2sample'][o].startswith("*")])
        n+=len(nsamples)
        rel=pos-node[0]
        for k in nsamples:
            v=ndata['offsets'][k]+rel
            qlpoint[k]=v
            qrpoint[k]=v+l
    return (l,n,qlpoint)

def maptooffsets(mums):
    mapping=dict()
    relmums=[]
    for mum in mums:
        relmum=lookup(mum)
        relmums.append(relmum)
        mapping[tuple(relmum[2].values())]=mum
    return relmums,mapping

def graphmumpicker(mums,idx,precomputed=False):
    try:
        if len(mums)==0:
            return

        if not precomputed:
            logging.debug("Selecting input multimums for %d samples: %d"%(idx.nsamples, len(mums)))
            mmums=[mum for mum in mums if len(mum[2])==idx.nsamples] #subset only those mums that apply to all indexed genomes/graphs

            if len(mmums)==0:
                logging.debug("No MUMS that span all input genomes, segment genomes.")
                mmums=segment(mums)
                logging.debug("Segmented genomes/graphs into %s, now %d MUMS for chaining."%(mmums[0][2].keys(),len(mmums)))

            if len(mmums)>maxmums and maxmums!=0:
                logging.debug("Number of MUMs exceeds cap (%d), taking largest %d"%(len(mums),maxmums))
                mmums.sort(key=lambda mum: mum[0]) #sort by size
                mmums=mmums[-maxmums:] #cap to topn mums

            logging.debug("Mapping indexed positions to relative postions within genomes.")
            relmums,mapping=maptooffsets(mmums)

            logging.debug("Subset to max available number of samples in set")
            relmums.sort(key=lambda m: m[1]) #sort by n
            
            relmums=[mum for mum in relmums if mum[2].keys()==relmums[-1][2].keys()] #subset to only those mums that apply to all genomes in the graph

            if len(relmums)==0:
                logging.debug("No MUMS that span all genomes that are in the graph, segment genomes.")
                relmums=segment(relmums)
                logging.debug("Segmented genomes/graphs into %s, now %d MUMS for chaining."%(relmums[0][2].keys(),len(relmums)))

            logging.debug("Left with %d mums"%len(relmums))

            if idx.leftnode!=None:
                spd=dict()
                for k in relmums[-1][2].keys():
                    spd[k]=G.node[idx.leftnode]['offsets'][k]+(idx.leftnode[1]-idx.leftnode[0])-1
                left=(0,0,spd)
            else:
                spd=dict()
                for sid in relmums[-1][2].keys():
                    spd[sid]=-1
                left=(0,0,spd)

            if idx.rightnode!=None:
                spd=dict()
                for k in relmums[-1][2].keys():
                    spd[k]=G.node[idx.rightnode]['offsets'][k]
                right=(0,0,spd)
            else:
                spd=dict()
                for sid in relmums[-1][2].keys():
                    spd[sid]=G.graph['id2end'][sid]
                right=(0,0,spd)

            logging.debug("Chaining %d mums"%len(mmums))
            chainedmums=chain(relmums,left,right,gcmodel=gcmodel)[::-1]
            logging.debug("Selected chain of %d mums"%len(chainedmums))
            if len(chainedmums)==0:
                return
            
            optsplit=None
            for mum in chainedmums: #determine optimal split in chain
                lseq=0
                rseq=0
                for crd in mum[2]:
                    lseq=mum[2][crd]
                    assert(lseq>=0)
                    rseq=right[2][crd]-mum[2][crd]+mum[0]
                    assert(rseq>=0)
                if optsplit==None or abs(lseq-rseq)<optsplit:
                    optsplit=abs(lseq-rseq)
                    splitmum=mum

        else:
            logging.debug("Selecting MUM from precomputed chain")
            assert(len(mums)>0)
            chainedmums,mapping=maptooffsets(mums)
            splitmum=chainedmums[len(chainedmums)/2] #take the middle
        
        logging.debug("Best MUM from chain: %s"%str(splitmum))
        skipleft=[]
        skipright=[]
        if seedsize>0:
            t=skipleft
            for mum in chainedmums:
                if mum==splitmum:
                    t=skipright
                    continue
                t.append(mapping[tuple(mum[2].values())])
            skipleft=[mum for mum in skipleft if mum[0]>=seedsize]
            skipright=[mum for mum in skipright if mum[0]>=seedsize]

        logging.debug("Mapped back to index space: %s"%str(splitmum))
        splitmum=mapping[tuple(splitmum[2].values())]
        logging.debug("Skipleft: %d"%len(skipleft))
        logging.debug("Skipright: %d"%len(skipright))
        return splitmum,skipleft,skipright

    except Exception:
        print traceback.format_exc()
        return

def printSA(index,maxline=100,start=0,end=None,fn="sa.txt"):
    sa=index.SA
    lcp=index.LCP
    t=index.T
    #so=index.SO
    if end==None:
        end=len(sa)
    with open(fn,'w') as f:
        f.write("%d\t%d\n"%(len(sa), len(lcp)))
        assert(len(sa)==len(lcp))
        for i in xrange(len(sa)):
            s=sa[i]
            lcpi=lcp[i]

            if i>0 and i<len(sa)-1:
                l1=lcp[i]
                l2=lcp[i+1]
            elif i==len(sa)-1:
                l1=max([lcp[i-1],lcp[i]])
                l2=0
            else:
                l1=0
                l2=lcp[i+1]

            if i>=start and i<=end:
                #f.write("%s\t%s\t%s\n"%(str(s).zfill(8), str(lcpi).zfill(6), t[s:s+maxline].ljust(maxline) if l1<=maxline else t[s:s+maxline]+"..."+t[s+l1-40:s+l1].ljust(maxline) ) )
                f.write("%s\t%s\t%s\t%s\t%s\n"%(str(s).zfill(8), str(lcpi).zfill(6), t[s:s+maxline] ,t[s+l1-maxline:s+l1], t[s+l2-maxline:s+l2] ) )




