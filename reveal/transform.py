import reveallib
import reveallib64
from utils import *
from multiprocessing.pool import Pool
import signal
import os
import math
import argparse
import logging
import numpy as np
import intervaltree
import sortedcontainers

try:
    from matplotlib import pyplot as plt
    from matplotlib import patches as patches
except:
    pass

def sparse_mums(r,q):
    mums=dict()
    # os.system("essamummer -qthreads 3 -F -r -c -l 20 -n %s %s -mum > mums"%(args.reference,args.contigs))
    os.system("essamummer -qthreads 3 -F -c -b -l 20 -n %s %s -mum > mums"%(r,q))
    # os.system("essamummer -qthreads 3 -F -c -l 20 -b -n %s %s -maxmatch > mums"%(args.reference,args.contigs))

    with open("mums") as mumfile:
        for line in mumfile:
            if line.startswith(">"):
                line=line.rstrip()
                if line.endswith("Reverse"):
                    qry=line.split()[1]
                    revcomp=1
                else:
                    qry=line.split()[1]
                    revcomp=0
            else:
                v=line.strip().split()
                mums[v[0]][qry][revcomp].append((int(v[1]),int(v[2]),int(v[3])))
    return mums

def plot(anchors,sep,wait=True,nc='r',rc='g',color="black",lines=False):
    
    if len(anchors)==0:
        return
    # elif len(anchors)>5000:
    #     logging.warn("More than 5000 anchors to display, plot only show the 5000 largest.")
    #     anchors=sorted(anchors,key=lambda a: a[5],reverse=True)[:5000]

    if len(anchors[0])==2: #unaligned blocks
        for start,stop in anchors:
            ax = plt.axes()
            if start<sep: #ref
                ax.add_patch(
                        patches.Rectangle(
                            (start, 0), #bottom left
                            stop-start, #width
                            sep, #height #should be qry length!
                            alpha=.25,
                            color="blue"
                        )
                    )
            else:
                ax.add_patch(
                        patches.Rectangle(
                            (0, start-sep), #bottom left
                            sep, #width
                            stop-start, #height
                            alpha=.25,
                            color="grey"
                        )
                    )
    elif len(anchors[0])==3: #mums
        for l,sps,revcomp in anchors:
            if revcomp:
                plt.plot( (sps[0],sps[0]+l), ((sps[1]-sep)+l, (sps[1]-sep)),'%s-'%rc)
            else:
                plt.plot( (sps[0],sps[0]+l), ((sps[1]-sep), (sps[1]-sep)+l),'%s-'%nc)
    elif len(anchors[0])==4: #synteny blocks, without orientation
        for anchor in anchors:
            s1,e1,s2,e2=anchor
            ax = plt.axes()
            ax.add_patch(
                    patches.Rectangle(
                        (s1, s2-sep), #bottom left
                        e1-s1, #width
                        e2-s2, #height
                        alpha=.5,
                        color=color
                    )
                )
    elif len(anchors[0])==5: #synteny blocks with orientation
        for anchor in anchors:
            s1,e1,s2,e2,revcomp=anchor
            ax = plt.axes()
            ax.add_patch(
                    patches.Rectangle(
                        (s1, s2-sep), #bottom left
                        e1-s1, #width
                        e2-s2, #height
                        alpha=.25,
                        color="green" if revcomp else "red"
                    )
                )
    elif len(anchors[0])==8: #synteny blocks with score and ctg
        
        if lines:
            rcxpoints,xpoints=[],[]
            rcypoints,ypoints=[],[]

            for anchor in anchors:
                s1,e1,s2,e2,revcomp,score,ref,ctg=anchor
                # plt.text(s1+((e1-s1)/2),(s2-sep)+(((e2-sep)-(s2-sep))/2) ,str(anchor),fontsize=6)
                if revcomp:
                    # plt.plot((s1,e1), (e2-sep,s2-sep),'g-')
                    rcxpoints.append(s1)
                    rcxpoints.append(e1)
                    rcxpoints.append(None)
                    rcypoints.append(e2-sep)
                    rcypoints.append(s2-sep)
                    rcypoints.append(None)
                else:
                    # plt.plot((s1,e1), (s2-sep,e2-sep),'r-')
                    xpoints.append(s1)
                    xpoints.append(e1)
                    xpoints.append(None)
                    ypoints.append(s2-sep)
                    ypoints.append(e2-sep)
                    ypoints.append(None)
            
            plt.plot(xpoints,ypoints,'r-')
            plt.plot(rcxpoints,rcypoints,'g-')

        else:
            for anchor in anchors:
                s1,e1,s2,e2,revcomp,score,ref,ctg=anchor
                ax = plt.axes()
                ax.add_patch(
                        patches.Rectangle(
                            (s1, s2-sep), #bottom left
                            e1-s1, #width
                            e2-s2, #height
                            alpha=.25,
                            color="green" if revcomp else "red"
                        )
                    )

    if wait:
        plt.show()
    else:
        plt.draw()

def addctginfo(mums,ctg2range):
    logging.debug("Augment contig information.")
    #add ref information to mums
    mums.sort(key=lambda m: m[1][0]) #sort mums by ref domain
    intvidx=0
    for i in range(len(mums)):
        while mums[i][1][0]>ctg2range[intvidx][1]:
            intvidx+=1
        mums[i]=mums[i]+(intvidx,)
    
    #add contig information to mums
    mums.sort(key=lambda m: m[1][1]) #sort mums by query domain
    intvidx=0
    for i in range(len(mums)):
        while mums[i][1][1]>ctg2range[intvidx][1]:
            intvidx+=1
        mums[i]=mums[i]+(intvidx,)
    logging.debug("Done.")
    return mums

def transform_cmd(args):
    for qry in args.contigs:
        logging.info("Running transform for %s"%qry)
        transform(args,qry)
        logging.info("Done")

def transform(args,qry):

    if args.output==None:
        prefix=os.path.splitext(os.path.basename(qry))[0]
    else:
        prefix=args.output
    
    refnames=[]
    ctgnames=[]

    if args.sa64:
        idx=reveallib64.index()
    else:
        idx=reveallib.index()

    ctg2range=[]
    for sample in [args.reference[0],qry]:
        idx.addsample(os.path.basename(sample))

        for name,seq in fasta_reader(sample, cutN=args.cutn):
            if len(seq)<args.minctglength:
                logging.debug("Skip transform for contig: %s"%name)
                continue

            intv=idx.addsequence(seq)
            ctg2range.append(intv)

            if sample==args.reference[0]:
                refnames.append(name)
            else:
                ctgnames.append(name)
    T=idx.T

    logging.info("Compute mums.")
    idx.construct(rc=False)
    mums=addctginfo(idx.getmums(args.minlength),ctg2range)
    logging.info("Done, %d mums."%len(mums))
    
    if args.cluster:
        logging.info("Cluster mums by diagonals.")
        blocks=clustermumsbydiagonal(mums,maxdist=args.maxdist,minclustsize=args.mincluster,rcmums=False)
        logging.info("Done, %d clusters."%len(blocks))
    else:
        blocks=[(mum[1][0], mum[1][0]+mum[0], mum[1][1], mum[1][1]+mum[0], mum[2], mum[0], mum[3], mum[4]) for mum in mums]
    
    # rcidx=idx.copy()
    # rcidx.construct(rc=True)
    # mums+=rcidx.getmums(args.minlength)

    logging.info("Compute RC mums.")
    idx.construct(rc=True)
    rcmums=addctginfo(idx.getmums(args.minlength),ctg2range)
    logging.info("Done, %d rc mums in total."%len(rcmums))


    sep=idx.nsep[0]
    idxn=idx.n

    rlength=idx.nsep[0]
    qlength=idxn-idx.nsep[0]

    del idx

    
    if args.cluster:
        logging.info("Cluster rc mums by anti-diagonals.")
        rcblocks=clustermumsbydiagonal(rcmums,maxdist=args.maxdist,minclustsize=args.mincluster,rcmums=True)
        logging.info("Done, %d rc clusters."%len(rcblocks))
    else:
        rcblocks=[(mum[1][0], mum[1][0]+mum[0], mum[1][1], mum[1][1]+mum[0], mum[2], mum[0], mum[3], mum[4]) for mum in rcmums]
    
    blocks+=rcblocks

    logging.info("Remove anchors that are contained in other clusters.")
    
    syntenyblocks=remove_contained_blocks(blocks)
    logging.info("Done, %d anchors remain."%len(syntenyblocks))

    if args.greedy:
        logging.info("Assign overlap between MUMs in a greedy manner.")
        syntenyblocks=remove_overlap_greedy_blocks(syntenyblocks)
    else:
        logging.info("Assign overlap between MUMs in a conservative manner.")
        syntenyblocks=remove_overlap_conservative_blocks(syntenyblocks)

    logging.info("Done.")

    logging.info("Remove anchors that after overlap removal are smaller then initial mum size.")
    syntenyblocks=[b for b in syntenyblocks if (b[1]-b[0]) >= args.minlength]
    logging.info("Done, %d blocks remain."%len(syntenyblocks))




    #encode mums as blocks
    # syntenyblocks=[(mum[1][0], mum[1][0]+mum[0], mum[1][1], mum[1][1]+mum[0], mum[2], mum[0], mum[3], mum[4]) for mum in mums]
    # syntenyblocks=blocks

    logging.info("Start glocal chaining for filtering anchors.")
    syntenyblocks=glocalchain(syntenyblocks,rlength,qlength,rearrangecost=args.rearrangecost,
                                                            inversioncost=args.inversioncost, 
                                                            useheap=args.useheap, 
                                                            lastn=args.lastn,
                                                            axis=0)
    # G=localcolinearchains(syntenyblocks,rlength,qlength,rearrangecost=rearrangecost,inversioncost=inversioncost)
    # chain,rcchain=colinearchains(syntenyblocks,rlength,qlength)
    logging.info("%d anchors remain after glocal chaining."%len(syntenyblocks))


    if args.plot:
        plot(syntenyblocks,sep,wait=False,lines=True)

    logging.info("Merge consecutive blocks.")
    syntenyblocks=merge_consecutive(syntenyblocks)
    logging.info("%d blocks after merging consecutive blocks."%len(syntenyblocks))

    logging.info("Remove all blocks that are shorter than minchainsum (%d)."%args.minchainsum)
    syntenyblocks=[b for b in syntenyblocks if b[5] >= args.minchainsum]
    logging.info("%d blocks after filtering for minchainsum."%len(syntenyblocks))
    
    weight,cost=chainscore(syntenyblocks, rearrangecost=args.rearrangecost,inversioncost=args.inversioncost) #determine the actual cost of the glocal chain 
    score=weight-cost
    syntenyblocks=merge_consecutive(syntenyblocks)

    if args.optimise and len(syntenyblocks)>1:
        iteration=1
        
        while True:
            iteration+=1
            logging.info("Optimise chain, iteration %d."%iteration)
            tsyntenyblocks,weight,cost=optimise(syntenyblocks,ctg2range,rearrangecost=args.rearrangecost,inversioncost=args.inversioncost)
            nscore=weight-cost
            if nscore<=score:
                break
            else:
                score=nscore
                syntenyblocks=tsyntenyblocks
                syntenyblocks=merge_consecutive(syntenyblocks)
        logging.info("Done. %d blocks after optimisation."%len(syntenyblocks))

    logging.debug("Extend %d blocks to query borders."%len(syntenyblocks))
    extendblocks(syntenyblocks,ctg2range)
    logging.debug("Done.")

    if args.plot:
        for start,end in ctg2range:
            if start<sep:
                plt.axvline(x=start, ymin=0, ymax=idxn-sep, linewidth=.1, linestyle='solid')
            else:
                plt.axhline(y=start-sep, xmin=0, xmax=sep, linewidth=.1, linestyle='solid')

        plot(syntenyblocks,sep,wait=False)
        plt.xlim(0,rlength)
        plt.ylim(0,qlength)

        if args.interactive:
            plt.show()
        else:
            plt.savefig("%s.png"%(prefix))

        plt.clf()

    #determine the subset of mappable contigs from ref and qry
    mappablectgs=set()
    for s1,e1,s2,e2,o,score,refid,ctgid in syntenyblocks:
        mappablectgs.add(ctgid)
        mappablectgs.add(refid)

    logging.info("Write breakpoint graph to: %s.gfa"%prefix)
    write_breakpointgraph(syntenyblocks,T,refnames,ctgnames,mappablectgs,prefix)


def clustermumsbydiagonal(mums,maxdist=90,minclustsize=65,rcmums=True):
    
    logging.debug("Sorting anchors by diagonals...")
    if rcmums:
        mums.sort(key=lambda m: (m[1][0]+(m[1][1]+m[0]), m[1][0]-(m[1][1]+m[0])) ) #sort mums by anti-diagonal, then diagonal
    else:
        mums.sort(key=lambda m: (m[1][0]-m[1][1], m[1][0]+m[1][1])) #sort mums by diagonal, then anti-diagonal
    logging.debug("Done.")

    l,sps,rc,ctg,ref=mums[0]
    clusters=[(sps[0],sps[0]+l,sps[1],sps[1]+l,rc,l,ctg,ref)]

    update_progress(0,len(mums))
    for i in xrange(1,len(mums)):
        update_progress(i,len(mums))

        l,sps,rc,ctg,ref=mums[i]
        s1,e1,s2,e2,prc,score,pctg,pref=clusters[-1]

        if rcmums:
            d=mums[i][1][0]+(mums[i][1][1]+mums[i][0])
            pd=e1+s2
        else:
            d=mums[i][1][0]-mums[i][1][1]
            pd=s1-s2

        if d==pd and pctg==ctg and pref==ref: #same diagonal and same contigs
            dist=mums[i][1][0]-e1
            assert(dist>=0)
            if dist < maxdist:
                if rc==0:
                    clusters[-1]=(s1,sps[0]+l,s2,sps[1]+l,rc,score+l,ctg,ref)
                else:
                    clusters[-1]=(s1,sps[0]+l,sps[1],e2,rc,score+l,ctg,ref)
            else:
                clusters.append((sps[0],sps[0]+l,sps[1],sps[1]+l,rc,l,ctg,ref))
        else:
            clusters.append((sps[0],sps[0]+l,sps[1],sps[1]+l,rc,l,ctg,ref))

    return [c for c in clusters if c[5]>=minclustsize]


def write_breakpointgraph(syntenyblocks,T,refnames,ctgnames,mappablectgs,outputprefix):
    #build a breakpoint graph, that we can write to GFA
    G=nx.MultiDiGraph()
    start=uuid.uuid4().hex
    end=uuid.uuid4().hex
    G.graph['startnodes']=[start]
    G.graph['endnodes']=[end]
    G.graph['paths']=[]
    G.graph['path2id']={}
    G.graph['id2path']={}

    G.add_node(start,offsets=dict())
    G.add_node(end,offsets=dict())
    
    pid=0
    for name in refnames:
        if pid in mappablectgs:
            # name=os.path.splitext(os.path.basename(reference))[0]+"_"+name
            name=outputprefix+"_"+name
            G.graph['paths'].append(name)
            G.graph['path2id'][name]=pid
            G.graph['id2path'][pid]=name
            G.node[start]['offsets'][pid]=0
        else:
            logging.info("No contigs were mapped to: %s"%name)
        pid+=1

    for name in ctgnames:
        if pid in mappablectgs:
            name="*"+name #prefix so we can recognise the two paths afterwards
            G.graph['paths'].append(name)
            G.graph['path2id'][name]=pid
            G.graph['id2path'][pid]=name
            G.node[start]['offsets'][pid]=0
        else:
            logging.info("Contig: %s could not be uniquely placed on the reference"%name)
        pid+=1

    #write the reference layout of the query sequences
    syntenyblocks.sort(key=lambda b: b[0]) #TODO: check if not already the case..
    prefid=None
    pnid=None
    l=0

    mapping=dict()
    nid=0
    for i,block in enumerate(syntenyblocks):
        s1,e1,s2,e2,o,score,refid,ctgid=block
        mapping[(s2,e2)]=nid

        if refid!=prefid:
            if prefid!=None:
                G.add_edge(pnid,end,paths=set([prefid]),ofrom="+", oto="+")
            pnid=start
            l=0

        if o==0:
            G.add_node(nid,seq=T[s2:e2],offsets={refid:l})
        else:
            G.add_node(nid,seq=rc(T[s2:e2]),offsets={refid:l})
        
        G.add_edge(pnid,nid,paths=set([refid]),ofrom="+", oto="+")
        prefid=refid
        pnid=nid
        nid+=1
        l+=e2-s2
        
        if i!=len(syntenyblocks)-1: #add gap node, so we later know which bubbles are caused by gaps in the assembly
            gapsize=1 #TODO: if specified use reference to add a gap
            G.add_node(nid,seq="N"*gapsize,offsets={refid:l})
            l+=gapsize
            G.add_edge(pnid,nid,paths=set([refid]),ofrom="+", oto="+")
            pnid=nid
            nid+=1

    G.add_edge(pnid,end,paths=set([refid]),ofrom="+", oto="+")

    writeorg=True
    if writeorg: #write the original layout of the query sequences, so we can reconstruct the input afterwards
        syntenyblocks.sort(key=lambda b: b[2])
        pctgid=None
        pnid=None
        
        l=0
        for nid,block in enumerate(syntenyblocks):
            s1,e1,s2,e2,o,score,refid,ctgid=block
            nid=mapping[(s2,e2)]
            if ctgid!=pctgid:
                if pctgid!=None:
                    G.add_edge(pnid,end,paths=set([pctgid]),ofrom="+" if o==0 else "-", oto="+")
                pnid=start
                l=0
                po=0
            
            G.node[nid]['offsets'][ctgid]=l
            
            l+=e2-s2
            G.add_edge(pnid,nid,paths=set([ctgid]),ofrom="+" if po==0 else "-", oto="+" if o==0 else "-")
            po=o
            pctgid=ctgid
            pnid=nid

        G.add_edge(pnid,end,paths=set([ctgid]),ofrom="+" if o==0 else "-", oto="+")

    write_gfa(G,None,outputfile=outputprefix+".gfa")


def merge_consecutive(syntenyblocks):
    if len(syntenyblocks)<2:
        return syntenyblocks
    #first merge consecutive blocks in the chain
    syntenyblocks.sort(key=lambda s: s[0]) #order by ref position
    qryorder = sorted(xrange(len(syntenyblocks)), key= lambda i: syntenyblocks[i][2]) #qry order
    qryorder_inv = sorted(xrange(len(syntenyblocks)), key=qryorder.__getitem__) #inverse qry order
    head=0
    
    for ri in xrange(1,len(syntenyblocks)):
        pblock=syntenyblocks[ri-1]
        block=syntenyblocks[ri]
        pqi=qryorder_inv[ri-1] #index within the qryorder of pblock
        qi=qryorder_inv[ri] #index within the qryorder of block
        ps1,pe2,ps2,pe2,po,pscore,prefid,pctgid=pblock #previous block on reference
        s1,e1,s2,e2,o,score,refid,ctgid=block
        es1,ee1,es2,ee2,eo,escore,erefid,ectgid=syntenyblocks[head]
        if ctgid==pctgid:
            if pqi+1==qi and o==po==0:
                syntenyblocks[head]=(es1,e1,es2,e2,eo,escore+score,erefid,ectgid)
            elif pqi-1==qi and o==po==1:
                syntenyblocks[head]=(es1,e1,s2,ee2,eo,escore+score,erefid,ectgid)
            else:
                head+=1
                syntenyblocks[head]=block
        else:
            head+=1
            syntenyblocks[head]=block
    
    while head!=ri:#len(syntenyblocks)-1:
        syntenyblocks.pop()
        head+=1

    return syntenyblocks


def extendblocks(syntenyblocks,ctg2range):

    syntenyblocks.sort(key=lambda s: s[0]) #order by reference position
    
    for i in xrange(len(syntenyblocks)):
        s1,e1,s2,e2,o,score,ref,ctg=syntenyblocks[i]
        
        if i==0: #first
            s1=ctg2range[ref][0]
        else:
            ps1,pe1,ps2,pe2,po,pscore,pref,pctg=syntenyblocks[i-1]
            if pref==ref:
                s1=pe1
            else:
                s1=ctg2range[ref][0]
        
        if i==len(syntenyblocks)-1: #last
            e1=ctg2range[ref][1]
        else:
            ns1,ne1,ns2,ne2,no,nscore,nref,nctg=syntenyblocks[i+1]
            if nref==ref:
                e1+=((ns1-e1)/2)
            else:
                e1=ctg2range[ref][1]

        assert(s1<e1)
        syntenyblocks[i]=(s1,e1,s2,e2,o,score,ref,ctg)

    syntenyblocks.sort(key=lambda s: s[2]) #order by qry position

    for i in xrange(len(syntenyblocks)):
        s1,e1,s2,e2,o,score,ref,ctg=syntenyblocks[i]
        if i==0: #first
            s2=ctg2range[ctg][0]
        else:
            ps1,pe1,ps2,pe2,po,pscore,pref,pctg=syntenyblocks[i-1]
            if pctg==ctg:
                s2=pe2
            else:
                s2=ctg2range[ctg][0]

        if i==len(syntenyblocks)-1: #last
            e2=ctg2range[ctg][1]
        else:
            ns1,ne1,ns2,ne2,no,nscore,nref,nctg=syntenyblocks[i+1]
            if nctg==ctg:
                e2+=((ns2-e2)/2)
            else:
                e2=ctg2range[ctg][1]
        
        assert(s2<e2)
        syntenyblocks[i]=(s1,e1,s2,e2,o,score,ref,ctg)

def optimise(syntenyblocks,ctg2range,rearrangecost=1000,inversioncost=1):

    orgchain=sorted(syntenyblocks,key=lambda c: c[5])
    maxchain=syntenyblocks
    maxchain_weight,maxchain_cost=chainscore(maxchain, rearrangecost=rearrangecost,inversioncost=inversioncost)
    maxchainscore=maxchain_weight-maxchain_cost

    stack=[]
    loglevel=logging.getLogger().getEffectiveLevel()

    if loglevel>logging.DEBUG:
        update_progress(0,len(orgchain))
    
    for i in xrange(len(orgchain)):
        if loglevel>logging.DEBUG:
            update_progress(i,len(orgchain))

        tmp=list(stack+orgchain[i+1:])
        weight,cost=chainscore(tmp, rearrangecost=rearrangecost,inversioncost=inversioncost)
        tmpchainscore=weight-cost

        if tmpchainscore<maxchainscore:
            stack.append(orgchain[i])
        else:
            logging.debug("Dropped block %s, gain: %d"%(orgchain[i],tmpchainscore-maxchainscore))
            maxchainscore=tmpchainscore
            maxchain=tmp
            maxchain_cost=cost
            maxchain_weight=weight

    return maxchain,maxchain_weight,maxchain_cost

def chainscore(chain, rearrangecost=1000, inversioncost=1):
    chain.sort(key=lambda s: s[0]) #order by reference position
    qryorder = sorted(xrange(len(chain)), key= lambda i: chain[i][2]) #qry order
    qryorder_inv = sorted(xrange(len(chain)), key=qryorder.__getitem__) #inverse qry order
    
    #count out of order traversals
    rearrangements=0
    inversions=0
    
    cost=0
    weight=chain[0][5]
    
    # costs=[0]
    for ri in xrange(1,len(chain)):

        pblock=chain[ri-1]
        block=chain[ri]

        ps1,pe1,ps2,pe2,po,pscore,pref,pctg=pblock
        s1,e1,s2,e2,o,score,ref2,ctg=block
        weight+=score

        # xgap=0#s1-pe1

        pqi=qryorder_inv[ri-1] #index within the qryorder of pblock
        qi=qryorder_inv[ri] #index within the qryorder of block

        if pctg==ctg and pref==ref2:
            
            if (pqi==qi-1) or (pqi==qi+1): #check if the two blocks are colinear

                gc=gapcost(pblock,block)
                cost+=gc

            else: #all other options use rearrangement penalty
                rearrangements+=1
                cost+=rearrangecost
        
        elif pctg!=ctg: #cross contigs
            
            if o==0:
                if qi>0:
                    pqs1,pqe1,pqs2,pqe2,pqo,pqscore,pq_ref,pq_ctg=chain[qryorder[qi-1]]
                else:
                    pq_ctg='start'
            else:
                if qi<len(qryorder)-1:
                    pqs1,pqe1,pqs2,pqe2,pqo,pqscore,pq_ref,pq_ctg=chain[qryorder[qi+1]]
                else:
                    pq_ctg='end'

            if po==0:
                if pqi<len(qryorder)-1:
                    nqs1,nqe1,nqs2,nqe2,nqo,nqscore,nq_ref,nq_ctg=chain[qryorder[pqi+1]]
                else:
                    nq_ctg='end'
            else:
                if pqi>0:
                    nqs1,nqe1,nqs2,nqe2,nqo,nqscore,nq_ref,nq_ctg=chain[qryorder[pqi-1]]
                else:
                    nq_ctg='start'

            if pq_ctg==ctg or nq_ctg==pctg: #there exists another block on this query contig before changing contigs, so has to be rearranged
                rearrangements+=1
                cost+=rearrangecost

    return weight*2,cost


def update_progress(i,n):
    fullbar=100
    if (i+1) % (n/fullbar if n>fullbar else 1)==0 or i+1==n:
        done=int(fullbar*((i+1)/float(n)))
        todo=fullbar-done
        sys.stdout.write('\r[%s%s]'%("#"*done," "*todo))
        if i+1==n:
            sys.stdout.write('\n')
        sys.stdout.flush()


def glocalchain(syntenyblocks,rlength,qlength, rearrangecost=1000, inversioncost=1, lastn=50, useheap=False, axis=0):

    sep=rlength
    
    start=(0,0,rlength,rlength,0,0,0,0)
    startrc=(0,0,rlength+qlength,rlength+qlength,1,0,0,0)

    end=(rlength-2,rlength-2,rlength+qlength,rlength+qlength,0,0,0,0)
    endrc=(rlength-2,rlength-2,rlength,rlength,1,0,0,0)

    syntenyblocks.append(end)
    syntenyblocks.append(endrc)

    if axis==0: #sort by ref
        c1,c2=0,2
    else: #sort by qry
        c1,c2=2,0

    syntenyblocks.sort(key=lambda s: (s[c1],s[5]) ) #order by reference position, then score

    # qryorder = sorted(xrange(len(syntenyblocks)), key= lambda i: syntenyblocks[c2][2]) #qry order
    # qryorder_inv = sorted(xrange(len(syntenyblocks)), key=qryorder.__getitem__) #inverse qry order

    if useheap:
        heap=sortedcontainers.SortedList()
        heap.add((0,start))
        heap.add((0,startrc))
    else:
        heap=[None]*(len(syntenyblocks)+2)
        heap[0]=(0,start)
        heap[1]=(0,startrc)

    G=dict()

    maxscore=None

    n=len(syntenyblocks)

    update_progress(0,n)

    for ri,block in enumerate(syntenyblocks):

        update_progress(ri,n)

        s1,e2,s2,e2,o,score,refid,ctgid=block
        score=score*2

        assert(score>=0)

        bestscore=None
        bestblock=None
        bestcost=0

        if useheap:
            iterator=heap[::-1]#[:lastn]
        else:
            s=((ri+2)-lastn) if ((ri+2)-lastn)>0 else 0
            assert(s<(ri+2))
            iterator=heap[s:ri+2][::-1]

        for cscore,pblock in iterator:
            
            ps1,pe2,ps2,pe2,po,pscore,prefid,pctgid=pblock

            if useheap:
                if bestscore!=None:
                    if cscore<=bestscore:
                        break
            
            if pctgid==ctgid:
                c=min((rearrangecost,gapcost(pblock,block,inversioncost=inversioncost)))
            else:
                c=rearrangecost
            
            assert(c>=0)

            if bestscore==None or cscore-c > bestscore:
                bestscore=cscore-c
                bestblock=pblock
                bestcost=c

        cscore=bestscore+score

        if useheap:
            heap.add((cscore,block))
        else:
            heap[ri+2]=(cscore,block)
        
        if maxscore==None or maxscore<cscore:
            maxscore=cscore
            maxnode=block

        G[block]=(bestblock,bestscore)

    node=G[end][0] if G[end][1]>G[endrc][1] else G[endrc][0]

    chain=[]
    while node!=start and node!=startrc:
        s1,e1,s2,e2,o,score,refid,ctgid=node
        nnode,score=G[node]
        chain.append(nnode)
        node=nnode

    chain.pop() #remove startnode from chain

    return chain

#unused
def colinearchains(syntenyblocks,rlength,qlength):

    sep=rlength
    heap=sortedcontainers.SortedList()

    start=(0,0,rlength,rlength,0,0,0,0)
    startrc=(0,0,rlength+qlength,rlength+qlength,1,0,0,0)

    end=(rlength-2,rlength-2,rlength+qlength,rlength+qlength,0,0,0,0)
    endrc=(rlength-2,rlength-2,rlength,rlength,1,0,0,0)

    syntenyblocks.append(start)
    syntenyblocks.append(startrc)
    syntenyblocks.append(end)
    syntenyblocks.append(endrc)

    syntenyblocks.sort(key=lambda s: (s[0],s[5])) #order by reference position
    qryorder = sorted(xrange(len(syntenyblocks)), key= lambda i: syntenyblocks[i][2]) #qry order
    qryorder_inv = sorted(xrange(len(syntenyblocks)), key=qryorder.__getitem__) #inverse qry order

    G=dict()

    maxscore=None

    for ri,block in enumerate(syntenyblocks):
        s1,e2,s2,e2,o,score,refid,ctgid=block
        qi=qryorder_inv[ri] #index within the qryorder of blocks
        assert(score>=0)

        if block==start or block==startrc:
            heap.add((score,block,ri,qi))
            continue

        bestscore=None
        bestblock=None
        bestcost=0

        for cscore,pblock,pri,pqi in heap[::-1]:
            ps1,pe2,ps2,pe2,po,pscore,prefid,pctgid=pblock

            if bestscore!=None:
                if cscore<=bestscore:
                    break
            
            if block[4]==pblock[4]==1: #both rc
                if pqi>qi:
                    c=gapcost(pblock,block)
            elif block[4]==pblock[4]==0: #both normal
                if pqi<qi:
                    c=gapcost(pblock,block)
            else: #different orientation
                c=-1

            if c>=0:
                if cscore-c > bestscore:
                    bestscore=cscore-c
                    bestblock=pblock
                    bestcost=c

        cscore=bestscore+score
        heap.add((cscore,block,ri,qi))
        
        if maxscore==None or maxscore<cscore:
            maxscore=cscore
            maxnode=block

        G[block]=(bestblock,bestscore,bestcost)

        # if bestblock[4]==0 and block[4]==0:
        #     plt.plot((bestblock[1],block[0]),(bestblock[3]-sep,block[2]-sep),'y--')
        # elif bestblock[4]==0 and block[4]==1:
        #     plt.plot((bestblock[1],block[0]),(bestblock[3]-sep,block[3]-sep),'y--')
        # elif bestblock[4]==1 and block[4]==0:
        #     plt.plot((bestblock[1],block[0]),(bestblock[2]-sep,block[2]-sep),'y--')
        # else:
        #     plt.plot((bestblock[1],block[0]),(bestblock[2]-sep,block[3]-sep),'y--')

    node=G[end][0]

    chain=[]
    while node!=start:
        s1,e1,s2,e2,o,score,refid,ctgid=node
        
        assert(o==0)

        if o==0:
            plt.plot((s1,e1),(s2-sep,e2-sep),'r-')
        else:
            plt.plot((s1,e1),(e2-sep,s2-sep),'r-')

        nnode,score,cost=G[node]

        # if nnode[4]==0 and node[4]==0:
        #     plt.plot((nnode[1],node[0]),(nnode[3]-sep,node[2]-sep),'r--')
        #     # plt.text(nnode[1],nnode[3]-sep,str(cost),fontsize=8)
        # elif nnode[4]==0 and node[4]==1:
        #     plt.plot((nnode[1],node[0]),(nnode[3]-sep,node[3]-sep),'r--')
        #     # plt.text(nnode[1],nnode[3]-sep,str(cost),fontsize=8)
        # elif nnode[4]==1 and node[4]==0:
        #     plt.plot((nnode[1],node[0]),(nnode[2]-sep,node[2]-sep),'r--')
        #     # plt.text(nnode[1],nnode[2]-sep,str(cost),fontsize=8)
        # else:
        #     plt.plot((nnode[1],node[0]),(nnode[2]-sep,node[3]-sep),'r--')
        #     # plt.text(nnode[1],nnode[2]-sep,str(cost),fontsize=8)

        chain.append(nnode)
        node=nnode

    chain.pop() #remove startnode from chain


    node=G[endrc][0]


    rcchain=[]
    while node!=startrc:
        s1,e1,s2,e2,o,score,refid,ctgid=node
        
        assert(o==1)

        if o==0:
            plt.plot((s1,e1),(s2-sep,e2-sep),'g-')
        else:
            plt.plot((s1,e1),(e2-sep,s2-sep),'g-')

        nnode,score,cost=G[node]

        # if nnode[4]==0 and node[4]==0:
        #     plt.plot((nnode[1],node[0]),(nnode[3]-sep,node[2]-sep),'r--')
        #     # plt.text(nnode[1],nnode[3]-sep,str(cost),fontsize=8)
        # elif nnode[4]==0 and node[4]==1:
        #     plt.plot((nnode[1],node[0]),(nnode[3]-sep,node[3]-sep),'r--')
        #     # plt.text(nnode[1],nnode[3]-sep,str(cost),fontsize=8)
        # elif nnode[4]==1 and node[4]==0:
        #     plt.plot((nnode[1],node[0]),(nnode[2]-sep,node[2]-sep),'r--')
        #     # plt.text(nnode[1],nnode[2]-sep,str(cost),fontsize=8)
        # else:
        #     plt.plot((nnode[1],node[0]),(nnode[2]-sep,node[3]-sep),'r--')
        #     # plt.text(nnode[1],nnode[2]-sep,str(cost),fontsize=8)

        rcchain.append(nnode)
        node=nnode

    rcchain.pop() #remove startnode from chain

    return chain,rcchain

#unused
def localcolinearchains(syntenyblocks,rlength,qlength):
    heap=sortedcontainers.SortedList()
    sep=rlength
    syntenyblocks.sort(key=lambda s: (s[0],s[5])) #order by reference position
    qryorder = sorted(xrange(len(syntenyblocks)), key= lambda i: syntenyblocks[i][2]) #qry order
    qryorder_inv = sorted(xrange(len(syntenyblocks)), key=qryorder.__getitem__) #inverse qry order

    G=dict()

    maxscore=None

    for ri,block in enumerate(syntenyblocks):
        s1,e2,s2,e2,o,score,refid,ctgid=block
        qi=qryorder_inv[ri] #index within the qryorder of blocks
        assert(score>=0)

        bestscore=None
        bestblock=None
        bestcost=0

        for cscore,pblock,pri,pqi in heap[::-1]:
            
            if bestscore!=None:
                if cscore<=bestscore:
                    break
            c=-1
            if block[4]==pblock[4]==1: #both rc
                if pqi>qi:
                    c=gapcost(pblock,block)
            elif block[4]==pblock[4]==0: #both normal
                if pqi<qi:
                    c=gapcost(pblock,block)
            else: #different orientation
                c=-1
            
            if c>=0: #can be chained in a valid way
                if (cscore-c)>score and (bestscore==None or (cscore-c)>bestscore):
                    bestscore=cscore-c
                    bestblock=pblock
                    bestcost=c

            pcscore=cscore

        if bestblock==None: #can't be chained to another block
            heap.add((score,block,ri,qi))
        else:
            cscore=bestscore+score
            heap.add((cscore,block,ri,qi))
            G[block]=(bestblock,bestscore,bestcost)
            # if bestblock[4]==0 and block[4]==0:
            #     plt.plot((bestblock[1],block[0]),(bestblock[3]-sep,block[2]-sep),'y--')
            # elif bestblock[4]==0 and block[4]==1:
            #     plt.plot((bestblock[1],block[0]),(bestblock[3]-sep,block[3]-sep),'y--')
            # elif bestblock[4]==1 and block[4]==0:
            #     plt.plot((bestblock[1],block[0]),(bestblock[2]-sep,block[2]-sep),'y--')
            # else:
            #     plt.plot((bestblock[1],block[0]),(bestblock[2]-sep,block[3]-sep),'y--')

    return G


def gapcost(block1,block2,inversioncost=0):
    if block1[4]==block2[4]==0: #both normal orientation
        if block2[2]>block1[2]:
            return abs((block2[0]-block1[1])-(block2[2]-block1[3]))
        else:
            return abs((block2[0]-block1[1])-(block1[2]-block2[3]))+inversioncost
    elif block1[4]==block2[4]==1: #both reverse comp orientation
        if block2[2]>block1[2]:
            return abs((block2[0]-block1[1])-(block2[2]-block1[3]))+inversioncost
        else:
            return abs((block2[0]-block1[1])-(block1[2]-block2[3]))
    elif block1[4]==1 and block2[4]==0: #balanced inversion
        if block2[2]>block1[2]:
            return abs((block2[0]-block1[1])-(block2[2]-block1[3]))+inversioncost
        else:
            return abs((block2[0]-block1[1])-(block1[2]-block2[3]))+inversioncost
    else: #balanced inversion
        assert(block1[4]==0 and block2[4]==1)
        if block2[2]>block1[2]:
            return abs((block2[0]-block1[1])-(block2[2]-block1[3]))+inversioncost
        else:
            return abs((block2[0]-block1[1])-(block1[2]-block2[3]))+inversioncost


def printSA(index,maxline=100,start=0,end=None,fn="sa.txt"):
    sa=index.SA
    lcp=index.LCP
    t=index.T
    #so=index.SO
    if end==None:
        end=len(sa)
    
    # with open(fn,'w') as f:
    sys.stdout.write("%d\t%d\n"%(len(sa), len(lcp)))
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
            sys.stdout.write("%s\t%s\t%s\t%s\t%s\n"%(str(s).zfill(8), str(lcpi).zfill(6), t[s:s+maxline] ,t[s+l1-maxline:s+l1], t[s+l2-maxline:s+l2] ) )




def remove_overlap_conservative_blocks(anchors):
    
    for coord in [0,2]:

        if len(anchors)<=1: #by definition no containment
            return anchors

        anchors.sort(key=lambda m: (m[coord], (m[coord+1]-m[coord])*-1)) #sort by start position, then -1*size
        
        _anchors=[anchors[0]]
        last=anchors[0]
        for anchor in anchors[1:]:
            if anchor[coord] < last[coord+1]: #overlap
                if anchor[coord+1]<=last[coord+1]: #contained
                    continue
            _anchors.append(anchor)
            last=anchor
        anchors=_anchors

        _anchors=[anchors[0]]
        for anchor in anchors[1:]:
            s1,e1,s2,e2,o,score,refid,ctgid=anchor
            ps1,pe1,ps2,pe2,po,pscore,prefid,pctgid=_anchors[-1]

            overlap=(_anchors[-1][coord+1]) - anchor[coord]
            pl=pe1-ps1

            if overlap > 0: #overlap
                
                if score<=overlap:
                    continue

                assert(score-overlap >= 0)

                if o==0:
                    anchor=(s1+overlap,e1,s2+overlap,e2,o,score-overlap if overlap<score else 0,refid,ctgid)
                else:
                    if coord==0:
                        anchor=(s1+overlap,e1,s2,e2-overlap,o,score-overlap if overlap<score else 0,refid,ctgid)
                    else:
                        anchor=(s1,e1-overlap,s2+overlap,e2,o,score-overlap if overlap<score else 0,refid,ctgid)

                assert(anchor[coord+1]>_anchors[-1][coord+1])

                while pl<=overlap or pscore<=overlap:
                    _anchors.pop()
                    ps1,pe1,ps2,pe2,po,pscore,prefid,pctgid=_anchors[-1]
                    overlap=(_anchors[-1][coord+1]) - anchor[coord]
                    if overlap<0:
                        break
                    pl=pe1-ps1

                if overlap>0:                    
                    assert(pscore-overlap >= 0)
                    if po==0:
                        _anchors[-1]=(ps1,pe1-overlap,ps2,pe2-overlap,po,pscore-overlap if overlap<pscore else 0,prefid,pctgid)
                    else:
                        if coord==0:
                            _anchors[-1]=(ps1,pe1-overlap, ps2+overlap,pe2, po,pscore-overlap if overlap<pscore else 0, prefid,pctgid)
                        else:
                            _anchors[-1]=(ps1+overlap,pe1,ps2,pe2-overlap,po,pscore-overlap if overlap<pscore else 0, prefid,pctgid)
            
            _anchors.append(anchor)

        anchors=_anchors

    return anchors

def remove_overlap_greedy_blocks(anchors):
    
    #TODO: remove duplicates!

    for coord in [0,2]:
        if len(anchors)<=1: #by definition no containment
            return anchors
        
        update_progress(0,len(anchors))

        anchors.sort(key=lambda m: (m[coord], (m[coord+1]-m[coord])*-1)) #sort by start position, then -1*size
        
        _anchors=[anchors[0]]
        last=anchors[0]
        for anchor in anchors[1:]:
            if anchor[coord] < last[coord+1]: #overlap
                if anchor[coord+1]<=last[coord+1]: #contained
                    continue
            _anchors.append(anchor)
            last=anchor
        anchors=_anchors

        _anchors=[anchors[0]]
        # for anchor in anchors[1:]:

        for i in xrange(1,len(anchors)):
            anchor=anchors[i]

            update_progress(i,len(anchors))

            s1,e1,s2,e2,o,score,refid,ctgid=anchor
            ps1,pe1,ps2,pe2,po,pscore,prefid,pctgid=_anchors[-1]
            pl=pe1-ps1

            overlap=(_anchors[-1][coord+1]) - anchor[coord]

            if overlap > 0: #overlap

                if pscore > score: #update current anchor

                    if score<=overlap:
                        continue

                    assert(score-overlap >= 0)

                    if o==0:
                        anchor=(s1+overlap,e1,s2+overlap,e2,o,score-overlap if overlap<score else 0,refid,ctgid)
                    else:
                        if coord==0:
                            anchor=(s1+overlap,e1,s2,e2-overlap,o,score-overlap if overlap<score else 0,refid,ctgid)
                        else:
                            anchor=(s1,e1-overlap,s2+overlap,e2,o,score-overlap if overlap<score else 0,refid,ctgid)

                    _anchors.append(anchor)
                else:

                    while pl<=overlap or pscore<=overlap:
                        _anchors.pop()
                        ps1,pe1,ps2,pe2,po,pscore,prefid,pctgid=_anchors[-1]
                        overlap=(_anchors[-1][coord+1]) - anchor[coord]
                        if overlap<0:
                            break
                        pl=pe1-ps1

                    if overlap>0:
                        
                        assert(pl>overlap)
                        assert(pscore>overlap)

                        assert(pscore-overlap >= 0)

                        if po==0:
                            _anchors[-1]=(ps1,pe1-overlap,ps2,pe2-overlap,po,pscore-overlap if overlap<pscore else 0,prefid,pctgid)
                        else:
                            if coord==0:
                                _anchors[-1]=(ps1,pe1-overlap, ps2+overlap,pe2, po,pscore-overlap if overlap<pscore else 0,prefid,pctgid)
                            else:
                                _anchors[-1]=(ps1+overlap,pe1,ps2,pe2-overlap,po,pscore-overlap if overlap<pscore else 0,prefid,pctgid)                    
                    _anchors.append(anchor)
            else:
                _anchors.append(anchor)

        anchors=_anchors
    return anchors

def remove_contained_blocks(anchors):
    #remove duplicates!

    for coord in [0,2]:
        logging.info("Remove overlap in %s dimension."%("first" if coord==0 else "second"))

        if len(anchors)<=1: #by definition no containment
            return anchors

        anchors.sort(key=lambda m: (m[coord], (m[coord+1]-m[coord])*-1) ) #sort by start position, then -1*size
        
        _anchors=[anchors[0]]
        last=anchors[0]

        update_progress(0,len(anchors))

        # for anchor in anchors[1:]:
        for i in xrange(1,len(anchors)):
            anchor=anchors[i]
            update_progress(i,len(anchors))

            if anchor[coord] < last[coord+1]: #overlap
                if anchor[coord+1]<=last[coord+1]: #contained
                    continue
            _anchors.append(anchor)
            last=anchor
        anchors=_anchors

    return anchors



#unused
def remove_overlap_greedy_mums(anchors):
    
    #remove duplicates!
    n=2

    for coord in range(n):
        if len(anchors)<=1: #by definition no containment
            return anchors

        anchors.sort(key=lambda m: (m[1][coord], m[0]*-1)) #sort by start position, then -1*size
        
        _anchors=[anchors[0]]
        last=anchors[0]
        for anchor in anchors[1:]:
            if anchor[1][coord] < last[1][coord]+last[0]: #overlap
                if anchor[1][coord]+anchor[0]<=last[1][coord]+last[0]: #contained
                    continue
            _anchors.append(anchor)
            last=anchor
        anchors=_anchors

        _anchors=[anchors[0]]
        for anchor in anchors[1:]:
            overlap=(_anchors[-1][1][coord]+_anchors[-1][0]) - anchor[1][coord]

            if overlap > 0: #overlap

                if _anchors[-1][0] > anchor[0]:
                    if anchor[2]==0:
                        anchor=(anchor[0]-overlap, (anchor[1][0]+overlap, anchor[1][1]+overlap), anchor[2])
                    else:
                        if coord==0:
                            anchor=(anchor[0]-overlap, (anchor[1][0]+overlap, anchor[1][1]), anchor[2])
                        else:
                            anchor=(anchor[0]-overlap, (anchor[1][0], anchor[1][1]+overlap), anchor[2])

                    _anchors.append(anchor)
                else:

                    while _anchors[-1][0]<=overlap and overlap>0:
                        _anchors.pop()
                        overlap=(_anchors[-1][1][coord]+_anchors[-1][0]) - anchor[1][coord]

                    if overlap>0:
                    
                        if _anchors[-1][2]==0:
                            _anchors[-1]=(_anchors[-1][0]-overlap,_anchors[-1][1],_anchors[-1][2]) #update stack
                        else:
                            if coord==0:
                                _anchors[-1]=(_anchors[-1][0]-overlap,_anchors[-1][1],_anchors[-1][2])
                            else:
                                _anchors[-1]=(_anchors[-1][0]-overlap,(_anchors[-1][1][0]+overlap, _anchors[-1][1][1]),_anchors[-1][2])
                    
                    _anchors.append(anchor)
            else:
                _anchors.append(anchor)

        anchors=_anchors
    return anchors

#unused
def remove_contained_mums(anchors):
    #remove duplicates!

    for coord in range(2):
        if len(anchors)<=1: #by definition no containment
            return anchors

        anchors.sort(key=lambda m: (m[1][coord], m[0]*-1)) #sort by start position, then -1*size
        
        _anchors=[anchors[0]]
        last=anchors[0]
        for anchor in anchors[1:]:
            if anchor[1][coord] < last[1][coord]+last[0]: #overlap
                if anchor[1][coord]+anchor[0]<=last[1][coord]+last[0]: #contained
                    continue
            _anchors.append(anchor)
            last=anchor
        anchors=_anchors

    return anchors

#unused
def remove_overlap_conservative_mums(anchors):
    
    #remove duplicates!
    n=2

    for coord in range(n):
        if len(anchors)<=1: #by definition no containment
            return anchors

        anchors.sort(key=lambda m: (m[1][coord], m[0]*-1)) #sort by start position, then -1*size
        
        _anchors=[anchors[0]]
        last=anchors[0]
        for anchor in anchors[1:]:
            if anchor[1][coord] < last[1][coord]+last[0]: #overlap
                if anchor[1][coord]+anchor[0]<=last[1][coord]+last[0]: #contained
                    continue
            _anchors.append(anchor)
            last=anchor
        anchors=_anchors

        _anchors=[anchors[0]]
        last=anchors[0]
        for anchor in anchors[1:]:
            if anchor[1][coord] < last[1][coord]+last[0]: #overlap
                
                assert(anchor[1][coord]+anchor[0] > last[1][coord]+last[0]) #may not be contained, as we filtered these out already
                
                overlap=(last[1][coord]+last[0])-anchor[1][coord]

                assert(overlap>=0)

                assert(anchor[0]>overlap)
                
                if anchor[2]==0:
                    anchor=(anchor[0]-overlap, (anchor[1][0]+overlap, anchor[1][1]+overlap), anchor[2])
                else:
                    if coord==0:
                        anchor=(anchor[0]-overlap, (anchor[1][0]+overlap, anchor[1][1]), anchor[2])
                    else:
                        anchor=(anchor[0]-overlap, (anchor[1][0], anchor[1][1]+overlap), anchor[2])

                # assert(last[0]>overlap)

                if last[2]==0:
                    _anchors[-1]=(last[0]-overlap,last[1],last[2]) #update last
                else:
                    if coord==0:
                        _anchors[-1]=(last[0]-overlap,(last[1][0], last[1][1]+overlap),last[2])
                    else:
                        _anchors[-1]=(last[0]-overlap,(last[1][0]+overlap, last[1][1]),last[2])

            if _anchors[-1][0]<=0:
                _anchors[-1]=anchor
            else:
                _anchors.append(anchor)

            last=anchor

        anchors=_anchors

    return anchors
