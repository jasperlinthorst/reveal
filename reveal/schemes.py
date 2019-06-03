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
from utils import mem2mums
import math

# from matplotlib import pyplot as plt

def chain(mums,left,right,gcmodel="sumofpairs"):
    if len(mums)==0:
        return []

    logging.debug("Number of anchors before chaining: %d",len(mums))
    
    #use one coordinate system for sorting
    ref=mums[0][2].keys()[0]
    # logging.trace("Ref is %s"%ref)
    mums.append(right)
    mums.sort(key=lambda mum: mum[2][ref]) #sort by reference dimension

    sp2mum=dict()
    for mum in mums:
        sp2mum[mum[2][ref]]=mum

    minscore=-1*utils.gapcost([left[2][k] for k in right[2]],[right[2][k] for k in right[2]])
    logging.debug("Initial cost is: %d"%minscore)

    start=left[2][ref]
    end=right[2][ref]

    link=dict()
    score=dict({left[2][ref]:0})
    
    active=[left]
    processed=[]

    for mum in mums:
        trace=False
        #active=[ep2mum[ep] for ep in utils.range_search(mumeptree,(0,0),[sp-1 for sp in mum[2]])].sort(key=lambda x: score[x], reverse=True)
        remove=[]
        for pmum in processed:
            for crd in pmum[2]:
                if pmum[2][crd]+pmum[0]>mum[2][crd]:
                    break
            else:
                active.append(pmum)
                remove.append(pmum)

        for r in remove:
            processed.remove(r)

        active.sort(key=lambda x: score[x[2][ref]], reverse=True) #sort active by score decreasing, kind of priority queue
        
        w=None
        for amum in active:
            for crd in amum[2]:
                if amum[2][crd]+amum[0]>mum[2][crd]:
                    break
            else:
                s=score[amum[2][ref]] + (args.wscore*(mum[0]*((mum[1]*(mum[1]-1))/2)))

                if w!=None:
                    if w > s: #as input is sorted by score
                        break

                penalty=utils.gapcost([amum[2][k]+amum[0] for k in mum[2]],[mum[2][k] for k in mum[2]],model=gcmodel)

                assert(penalty>=0)

                # tmpw=score[amum[2][ref]] + (args.wscore*(mum[0]*((mum[1]*(mum[1]-1))/2))) - (args.wpen*penalty)
                tmpw=s - (args.wpen*penalty)

                if tmpw>w or w==None:
                    logging.trace("mum: %s --> %s = penalty: %d and score at amum: %d, score at mum: %d"%(str(mum),str(amum),penalty,s,tmpw))
                    w=tmpw
                    best=amum

        link[mum[2][ref]]=best[2][ref]

        score[mum[2][ref]]=w

        processed.append(mum)

    logging.debug("Best score is: %d"%score[end])
    logging.trace("Min score is: %d"%minscore)

    #backtrack
    path=[]
    while end!=start:
        path.append((sp2mum[end],score[end]))
        end=link[end]

    return path[1:]

#determine a subset of genomes for which (length * n) is largest
def segment(mums):
    d=dict()
    for mum in mums:
        k=tuple(sorted([gid for gid,sp in mum[2]]))
        if k in d:
            d[k].append(mum)
        else:
            d[k]=[mum]

    best=0
    for part in d:
        z=sum([m[0] for m in d[part]])*len(part)
        if z>best:
            best=z
            partition=part

    logging.debug("Splitting input genomes: %s"%str(partition))
    return d[partition]

def lookup(mum):
    l,mmn,spd=mum
    if isinstance(spd,dict):
        sp=spd.values()
    elif isinstance(spd,tuple):
        sp=[sp for gid,sp in spd]
    else:
        logging.fatal("Unknown format: %s"%str(spd))

    n=0
    qlpoint=dict()
    qrpoint=dict()
    for pos in sp:
        t=ts[pos]
        assert(len(t)==1)
        node=iter(t).next()
        ndata=G.node[node]
        nsamples=set([o for o in ndata['offsets'].keys() if not G.graph['id2path'][o].startswith("*")])
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

def trim_overlap(mums):
    coords=mums[0][2]
    for coord in range(len(coords)):
        if len(mums)<=1: #by definition no more overlaps
            break

        mums.sort(key=lambda m: (m[2][coord][1],-m[0])) #sort by start position, then -1*size

        #filter the partial matches that are now contained
        mums=[mum for i,mum in enumerate(mums) if (i==0 and mums[i+1][2][coord][1]+mums[i+1][0] > mum[2][coord][1]+mum[0] ) or mums[i-1][2][coord][1]+mums[i-1][0]<mum[2][coord][1]+mum[0]]

        if len(mums)<=1: #by definition no more overlaps
            break

        trimmed=[mums[0]]
        for mum in mums[1:]:
            pmum=trimmed[-1]
            overlap = (pmum[2][coord][1]+pmum[0]) - mum[2][coord][1]
            if overlap>0:
                if pmum[0]-overlap>0:
                    trimmed[-1] = (pmum[0]-overlap, pmum[1], pmum[2])
                else:
                    del trimmed[-1]
                if mum[0]-overlap>0:
                    trimmed.append( (mum[0]-overlap, mum[1], tuple((k,v+overlap) for k,v in mum[2]) ))
            else:
                trimmed.append(mum)

        mums=trimmed

    return mums

args=None
splitchain="largest"
maxdepth=None #stop recursion when max depth is reached

def graphmumpicker(mums,idx,precomputed=False,minlength=0):
    try:
        if len(mums)==0:
            return ()
        
        if not precomputed:
            if maxdepth!=None:
                if idx.depth>maxdepth:
                    return ()

            if args.maxsize!=None:
                rpaths=[p for p in G.graph['paths'] if not p.startswith('*')]

                if idx.leftnode==None:
                    lo={G.graph['path2id'][p]: 0 for p in rpaths}
                else:
                    lo={k: G.node[idx.leftnode]['offsets'][k]+(idx.leftnode[1]-idx.leftnode[0]) for k in G.node[idx.leftnode]['offsets']}
                
                if idx.rightnode==None:
                    ro={G.graph['path2id'][p]: G.graph['id2end'][G.graph['path2id'][p]] for p in rpaths}
                else:
                    ro=G.node[idx.rightnode]['offsets']

                for k in set(lo.keys()) & set(ro.keys()):
                    if ro[k]-lo[k]>args.maxsize:
                        break
                else:
                    return () #no break, so all fragments in bubbles are smaller than maxsize

            logging.debug("Selecting input multimums (for %d samples) out of: %d mums"%(idx.nsamples, len(mums)))
            mmums=[mum for mum in mums if mum[1]==idx.nsamples] #subset only those mums that apply to all indexed genomes/graphs
            
            if len(mmums)==0 and idx.nsamples>2:
                logging.debug("No MUMS that span all input genomes, segment genomes.")
                mmums=segment(mums)
                logging.debug("Segmented genomes/graphs into %s, now %d MUMS for chaining."%(mmums[0][2],len(mmums)))

            if args.trim:
                logging.debug("Trimming overlap between mums.")
                mmums=trim_overlap(mmums)
                if len(mmums)==0:
                    return ()

            mmums.sort(key=lambda mum: mum[0], reverse=True) #sort by size

            logging.debug("Mapping indexed positions to relative postions within genomes.")

            relmums,mapping=maptooffsets(mmums) #and convert tuple to dict for fast lookup in chaining

            logging.debug("Subset to same group of samples")
            relmums.sort(key=lambda m: (m[1],m[0])) #sort by n, than l
            relmums=[mum for mum in relmums if mum[2].keys()==relmums[-1][2].keys()] #subset to only those mums that apply to the same set
            
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

            # if minlength==0: #autodetermine significant subset
                # relmums=[mum for mum in relmums if 1-((1-((.25**(mum[1]-1))**mum[0]))**o)<pcutoff] #subset to only significant mums
            
            if len(relmums)==0:
                logging.debug("No more significant MUMs.")
                return ()

            skipleft=[]
            skipright=[]

            if len(relmums)==1:
                splitmum=relmums[0]
            else:
                if len(relmums)>args.maxmums:
                    logging.debug("Number of MUMs exceeds cap (%d), taking largest %d"%(len(mmums),args.maxmums))
                    relmums=relmums[-args.maxmums:]

                logging.debug("Chaining %d mums"%len(relmums))
                chainedmums=chain(relmums,left,right,gcmodel=args.gcmodel)[::-1]

                logging.debug("Selected chain of %d mums"%len(chainedmums))
                if len(chainedmums)==0:
                    return ()

                if splitchain=="balanced":
                    logging.debug("Selecting MUM from chain on position within chain.")
                    optsplit=None
                    for mum,score in chainedmums: #determine optimal split in chain
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
                elif splitchain=="largest":
                    logging.debug("Selecting MUM from chain based on size.")
                    splitmum=sorted(chainedmums,key=lambda m:m[0][0])[-1][0]
                else: #select at random
                    logging.debug("Selecting MUM from chain at random.")
                    splitmum=chainedmums[random.randint(0,len(chainedmums)-1)][0]

                if args.seedsize>0:
                    t=skipleft
                    scoreatsplit=0
                    for mum,score in chainedmums:
                        if mum==splitmum:
                            scoreatsplit=score
                            t=skipright
                            continue
                        t.append( (mapping[tuple(mum[2].values())], score-scoreatsplit) )
                        # t.append( mapping[tuple(mum[2].values())] )
                    skipleft=[(mum,score) for mum,score in skipleft if mum[0]>=args.seedsize]
                    skipright=[(mum,score) for mum,score in skipright if mum[0]>=args.seedsize]

            splitmum=mapping[tuple(splitmum[2].values())]

            if minlength==0: #experimental, use significance to determine valid anchor length when minlength is set to 0
                o=1
                for p in left[2]:
                    o=o*(right[2][p]-left[2][p])
                l=splitmum[0]
                n=splitmum[1]
                p=((.25**(n-1)))**l #probability of observing this match by random chance
                if p>0:
                    p=1-math.exp(math.log(1-p) * o) #correct for the number of tests we actually did
                if p>args.pcutoff:
                    logging.info("P-value for: %s (n=%d l=%d o=%d) is %.4g"%(str(splitmum),n,l,o,p))
                    return ()
        else:
            logging.debug("Selecting MUM from precomputed chain")
            chainedmums=mums
            splitmum=chainedmums[len(chainedmums)/2][0]
            skipleft=chainedmums[:len(chainedmums)/2]
            skipright=chainedmums[(len(chainedmums)/2)+1:]
        
        logging.debug("Best MUM has length: %d"%splitmum[0])
        
        logging.debug("Skipleft: %d"%len(skipleft))
        logging.debug("Skipright: %d"%len(skipright))

        return splitmum,skipleft,skipright

    except Exception:
        logging.fatal(traceback.format_exc())

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




