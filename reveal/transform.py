import reveallib
import reveallib64
from utils import *
from multiprocessing.pool import Pool
import signal
import os
import math
import argparse
import logging
import intervaltree
import sortedcontainers
import time

def plot(plt,anchors,sep,wait=True,nc='r',rc='g',color=None,edges=False,lines=False,alpha=1,args=None):
    
    if len(anchors)==0:
        return

    if len(anchors[0])==2: #unaligned blocks
        for start,stop in anchors:
            ax = plt.axes()
            if start<sep: #ref
                ax.add_patch(
                        matplotlib.patches.Rectangle(
                            (start, 0), #bottom left
                            stop-start, #width
                            sep, #height #should be qry length!
                            alpha=.25,
                            color="blue"
                        )
                    )
            else:
                ax.add_patch(
                        matplotlib.patches.Rectangle(
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
                plt.plot( (sps[0],sps[0]+l), ((sps[1]-sep)+l, (sps[1]-sep)),'%s-'%rc,alpha=alpha)
            else:
                plt.plot( (sps[0],sps[0]+l), ((sps[1]-sep), (sps[1]-sep)+l),'%s-'%nc,alpha=alpha)
    elif len(anchors[0])==4: #synteny blocks, without orientation
        for anchor in anchors:
            s1,e1,s2,e2=anchor
            ax = plt.axes()
            ax.add_patch(
                    matplotlib.patches.Rectangle(
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
                    matplotlib.patches.Rectangle(
                        (s1, s2-sep), #bottom left
                        e1-s1, #width
                        e2-s2, #height
                        alpha=.25,
                        color="green" if revcomp else "red"
                    )
                )
    elif len(anchors[0])==8: #synteny blocks with score and ctg
        
        if edges:
            for c in [0,2]:
                anchors.sort(key=lambda a:a[c])

                xedges,yedges=[],[]
                
                panchor=None
                for anchor in anchors:

                    s1,e1,s2,e2,revcomp,score,ref,ctg=anchor
                    
                    if panchor!=None:
                        ps1,pe1,ps2,pe2,prevcomp,pscore,pref,pctg=panchor

                        if pctg!=ctg and pref!=ref:
                            panchor=anchor
                            continue

                        if c==0:
                            xedges.append(pe1)
                            xedges.append(s1)
                            xedges.append(None)
                            
                            if prevcomp:
                                yedges.append(ps2-sep)
                            else:
                                yedges.append(pe2-sep)
                            
                            if revcomp:
                                yedges.append(e2-sep)
                            else:
                                yedges.append(s2-sep)
                            yedges.append(None)

                        else:

                            if prevcomp:
                                xedges.append(ps1)
                            else:
                                xedges.append(pe1)

                            if revcomp:
                                xedges.append(e1)
                            else:
                                xedges.append(s1)

                            xedges.append(None)
                            
                            yedges.append(pe2-sep)
                            yedges.append(s2-sep)
                            yedges.append(None)

                    panchor=anchor

                if c==0:
                    plt.plot(xedges,yedges,'b--',alpha=alpha)
                else:
                    plt.plot(xedges,yedges,'y--',alpha=alpha)

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
            
            plt.plot(xpoints,ypoints,'r-' if color==None else '%s-'%color,alpha=alpha)
            plt.plot(rcxpoints,rcypoints,'g-' if color==None else '%s-'%color,alpha=alpha)

        else: #plot squares
            for anchor in anchors:
                s1,e1,s2,e2,revcomp,score,ref,ctg=anchor
                ax = plt.axes()
                ax.add_patch(
                        matplotlib.patches.Rectangle(
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

    if not args.interactive:
        matplotlib.use("agg")

    if args.plot:
        from matplotlib import pyplot as plt

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
    logging.info("Done, %d rc mums."%len(rcmums))

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

    if args.plot:
        plot(plt,blocks,sep,wait=False,lines=True,alpha=0.2,args=args)

    # if args.plot:
    #     plot(blocks,sep,wait=False,lines=True)

    logging.info("Start glocal chaining for filtering anchors (reference).")
    
    # blocks.sort(key=lambda b: b[1]-b[0])
    # logging.info("Largest ref block: %s"%str(blocks[-1]))
    # minbacktrack=blocks[-1][1]-blocks[-1][0]

    nbefore=len(blocks)
    syntenyblocks=blocks
    nafter=None
    refiteration=0
    # maxiter=1
    while nbefore!=nafter:# or refiteration==maxiter:
        logging.info("Glocal chain iteration %d"%refiteration)
        nbefore=len(syntenyblocks)
        syntenyblocks=glocalchain(syntenyblocks,rlength,qlength,ctg2range,rearrangecost=args.rearrangecost,
                                                            inversioncost=args.inversioncost,
                                                            _lambda=args._lambda,
                                                            eps=args.eps,
                                                            useheap=args.useheap, 
                                                            lastn=args.lastn,
                                                            lastbp=args.lastbp,
                                                            alfa=args.alfa,
                                                            gapopen=args.gapopen,
                                                            axis=0)
        nafter=len(syntenyblocks)
        logging.info("Anchor before chaining: %s"%nbefore)
        logging.info("Anchor after chaining: %s"%nafter)
        refiteration+=1
        
        if args.plot:
            plot(plt,syntenyblocks,sep,wait=False,lines=True,color='k',alpha=.7)
    
    logging.info("%d anchors remain after glocal chaining (reference)."%len(syntenyblocks))

    logging.info("Start glocal chaining for filtering anchors (query).")

    nbefore=len(syntenyblocks)
    nafter=None
    qryiteration=0
    while nbefore!=nafter:# or qryiteration==maxiter:
        logging.info("Glocal chain iteration %d"%qryiteration)
        nbefore=len(syntenyblocks)
        syntenyblocks=glocalchain(syntenyblocks,rlength,qlength,ctg2range,rearrangecost=args.rearrangecost,
                                                                inversioncost=args.inversioncost,
                                                                _lambda=args._lambda,
                                                                eps=args.eps,
                                                                useheap=args.useheap, 
                                                                lastn=args.lastn,
                                                                lastbp=args.lastbp,
                                                                alfa=args.alfa,
                                                                gapopen=args.gapopen,
                                                                axis=1)
        nafter=len(syntenyblocks)
        logging.info("Anchor before chaining: %s"%nbefore)
        logging.info("Anchor after chaining: %s"%nafter)
        qryiteration+=1

    # G=localcolinearchains(syntenyblocks,rlength,qlength,rearrangecost=rearrangecost,inversioncost=inversioncost)
    # chain,rcchain=colinearchains(syntenyblocks,rlength,qlength)
    logging.info("%d anchors remain after glocal chaining (query)."%len(syntenyblocks))

    if args.plot:
        plot(plt,syntenyblocks,sep,wait=False,lines=True,color='b',alpha=.7)

    #take the intersection of both the chains
    # logging.info("Determine intersection between the chains...")
    # syntenyblocks=list(set(rsyntenyblocks) & set(qsyntenyblocks))
    # logging.info("Done. %d chains remain."%len(qsyntenyblocks))

    # logging.info("Remove anchors that are contained in other clusters."
    # syntenyblocks=remove_contained_blocks(blocks)
    # logging.info("Done, %d anchors remain."%len(syntenyblocks))
    # logging.info("Done.")

    
    logging.info("Merge consecutive blocks.")
    syntenyblocks=merge_consecutive(syntenyblocks)
    logging.info("%d blocks after merging consecutive blocks."%len(syntenyblocks))

    # if args.plot:
    #     plot(syntenyblocks,sep,wait=True,lines=True,color='b')

    # logging.info("Merge consecutive blocks.")
    # syntenyblocks=merge_consecutive(syntenyblocks)
    # logging.info("%d blocks after merging consecutive blocks."%len(syntenyblocks))

    if args.greedy:
        logging.info("Assign overlap between MUMs in a greedy manner.")
        syntenyblocks=remove_overlap_greedy_blocks(syntenyblocks)
        logging.info("Done.")
    else:
        logging.info("Assign overlap between MUMs in a conservative manner.")
        syntenyblocks=remove_overlap_conservative_blocks(syntenyblocks)
        logging.info("Done.")
    
    logging.info("Remove all blocks that are shorter than minchainsum (%d)."%args.minchainsum)
    syntenyblocks=[b for b in syntenyblocks if b[5] >= args.minchainsum]
    logging.info("%d blocks after filtering for minchainsum."%len(syntenyblocks))

    logging.info("Merge consecutive blocks.")
    syntenyblocks=merge_consecutive(syntenyblocks)
    logging.info("%d blocks after merging consecutive blocks."%len(syntenyblocks))

    if args.optimise and len(syntenyblocks)>1:

        weight,cost,edgecosts=chainscore(syntenyblocks, rlength, qlength, ctg2range,rearrangecost=args.rearrangecost,inversioncost=args.inversioncost,_lambda=args._lambda,eps=args.eps,alfa=args.alfa,gapopen=args.gapopen) #determine the actual cost of the glocal chain 
        score=weight-cost

        assert(len(edgecosts) == len(syntenyblocks)+1)

        iteration=0
        
        while True:
            iteration+=1
            logging.info("Optimise chain, iteration %d."%iteration)
            tsyntenyblocks,tweight,tcost,tedgecosts=optimise(syntenyblocks,rlength, qlength, ctg2range,rearrangecost=args.rearrangecost,inversioncost=args.inversioncost,_lambda=args._lambda,eps=args.eps,alfa=args.alfa,gapopen=args.gapopen)
            nscore=tweight-tcost
            
            if nscore<=score:
                break
            else:
                score=nscore
                syntenyblocks=tsyntenyblocks
                weight=tweight
                cost=tcost
                edgecosts=tedgecosts
                syntenyblocks=merge_consecutive(syntenyblocks)

        logging.info("Done. %d blocks after optimisation."%len(syntenyblocks))

    syntenyblocks=merge_consecutive(syntenyblocks)
    weight,cost,edgecosts=chainscore(syntenyblocks, rlength, qlength, ctg2range,rearrangecost=args.rearrangecost,inversioncost=args.inversioncost,_lambda=args._lambda,eps=args.eps,alfa=args.alfa,gapopen=args.gapopen) #determine the actual cost of the glocal chain 
    score=weight-cost

    assert(len(edgecosts) == len(syntenyblocks)+1)

    if args.outputbed: #before extending to the edges of the contig, output the breakpoint regions

        logging.info("Write bedfile with contig mappings on reference to: %s.bed"%prefix)
        with open(prefix+".bed",'w') as bedout:

            block2ctgidx=dict()
            pctgid=None

            ctgid2lastblock=dict()
            ci=0

            syntenyblocks.sort(key=lambda b: b[2]) #sort by query
            for i,block in enumerate(syntenyblocks): #sorted by query
                s1,e1,s2,e2,o,score,refid,ctgid=block
                if ctgid!=pctgid:
                    if pctgid!=None:
                        ctgid2lastblock[pctgid]=ci
                    ci=0
                else:
                    ci+=1
                block2ctgidx[block]=ci
                pctgid=ctgid

            ctgid2lastblock[pctgid]=ci

            syntenyblocks.sort(key=lambda b: b[0]) #sort by reference
            bedout.write("#reference\trefbegin\trefend\tcontig:segmentidx:lastsegmentidx:begin:end\tscore:cost\torientation\taln-start\taln-end\n")

            pblock=None

            for i,block in enumerate(syntenyblocks): #sorted by reference
                s1,e1,s2,e2,o,score,refid,ctgid=block
                
                if i>0:
                    ps1,pe1,ps2,pe2,po,pscore,prefid,pctgid=pblock
                else:
                    pblock=None

                cost=edgecosts[i] #cost to connect to pblock to block

                if i<len(syntenyblocks)-2:
                    nblock=syntenyblocks[i+1]
                    ns1,ne1,ns2,ne2,no,nscore,nrefid,nctgid=nblock
                else:
                    nblock=None

                ctgoffsets=ctg2range[ctgid]
                refoffsets=ctg2range[refid]

                if pblock!=None and prefid==refid:
                    start=(s1-refoffsets[0])-((s1-pe1)/2)
                else:
                    start=s1-refoffsets[0]

                if nblock!=None and nrefid==refid:
                    end=(e1-refoffsets[0])+((ns1-e1)/2)
                else:
                    end=e1-refoffsets[0]

                qstart=s2-ctgoffsets[0]
                qend=e2-ctgoffsets[0]

                chromname=refnames[refid].split()[0]

                qi=block2ctgidx[block]
                bedout.write("%s\t%d\t%d\t%s:%d:%d:%d:%d\t%d:%d\t%s\t%d\t%d\n"%(chromname, #chrom
                                                                start, #start
                                                                end, #end
                                                                ctgnames[ctgid-len(refnames)].split()[0], #name, make sure there's no whitespace to comply with bed 'format'
                                                                qi,
                                                                ctgid2lastblock[ctgid],
                                                                qstart,
                                                                qend,
                                                                score,
                                                                cost,
                                                                '+' if o==False else '-', #strand
                                                                s1-refoffsets[0], #thick start
                                                                e1-refoffsets[0]) #thick end
                                                                #itemRgb
                                                                #blockCount
                                                                #blockSizes
                                                                #blockStarts
                                                            )

                #bedout.write("%s\t%d\t%d\t%s\t%s\t%s\t%s\n"%(refnames[refid], pe1-refoffsets[0], s1-refoffsets[0], ctgnames[ctgid-len(refnames)], ctgnames[pctgid-len(refnames)], 'n' if po==False else 'r', 'n' if o==False else 'r'))

                pblock=block

    if args.plot:
        plot(plt,syntenyblocks,sep,wait=False,args=args)

    logging.debug("Extend %d blocks to query borders."%len(syntenyblocks))
    extendblocks(syntenyblocks,ctg2range)
    logging.debug("Done.")

    if args.plot:
        for start,end in ctg2range:
            if start<sep:
                plt.axvline(x=start, ymin=0, ymax=idxn-sep, linewidth=.1, linestyle='solid')
            else:
                plt.axhline(y=start-sep, xmin=0, xmax=sep, linewidth=.1, linestyle='solid')

        plot(plt,syntenyblocks,sep,wait=False,edges=False,args=args)
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

    if len(mappablectgs)!=0:
        logging.info("Write breakpoint graph to: %s.gfa"%prefix)
        write_breakpointgraph(syntenyblocks,T,refnames,ctgnames,mappablectgs,prefix)
    else:
        logging.info("No mappable contigs.")

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
            name=os.path.basename(outputprefix+"_"+name)
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

    write_gfa(G,None,outputfile=outputprefix if outputprefix.endswith(".gfa") else outputprefix+".gfa")

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

def optimise(syntenyblocks,rlength, qlength, ctg2range,rearrangecost=1000,inversioncost=1,_lambda=5,eps=1,alfa=1,gapopen=10):

    orgchain=sorted(syntenyblocks,key=lambda c: c[5])
    maxchain=syntenyblocks
    maxchain_weight,maxchain_cost,maxchain_edgecosts=chainscore(maxchain, rlength, qlength, ctg2range, rearrangecost=rearrangecost,inversioncost=inversioncost,_lambda=_lambda,eps=eps,alfa=alfa,gapopen=gapopen)
    maxchainscore=maxchain_weight-maxchain_cost

    stack=[]
    loglevel=logging.getLogger().getEffectiveLevel()

    if loglevel>logging.DEBUG:
        update_progress(0,len(orgchain))
    
    for i in xrange(len(orgchain)):
        if loglevel>logging.DEBUG:
            update_progress(i,len(orgchain))

        tmp=list(stack+orgchain[i+1:])
        weight,cost,edgecosts=chainscore(tmp, rlength, qlength, ctg2range, rearrangecost=rearrangecost,inversioncost=inversioncost,_lambda=_lambda,eps=eps,alfa=alfa,gapopen=gapopen)
        tmpchainscore=weight-cost

        if tmpchainscore<maxchainscore:
            stack.append(orgchain[i]) #keep it
        else:
            logging.debug("Dropped block %s, gain: %d"%(orgchain[i],tmpchainscore-maxchainscore))
            maxchainscore=tmpchainscore
            maxchain=tmp
            maxchain_cost=cost
            maxchain_weight=weight
            maxchain_edgecosts=edgecosts

    logging.debug("Optimal chain has %d blocks and scores: %d"%(len(maxchain),maxchainscore))

    return maxchain,maxchain_weight,maxchain_cost,maxchain_edgecosts

def chainscore(chain, rlength, qlength, ctg2range, rearrangecost=1000, inversioncost=1, _lambda=5, eps=1, alfa=1, gapopen=10):
    
    # logging.debug("rearrangecost=%d, inversioncost=%d, _lambda=%d, eps=%d, alfa=%d, gapopen=%d"%(rearrangecost, inversioncost, _lambda, eps, alfa, gapopen))

    if len(chain)==0:
        start=(0,0,rlength,rlength,0,0,0,0)
        end=(rlength,rlength,rlength+qlength,rlength+qlength,0,0,0,0)
        cost=gapcost(start,end,rearrangecost=rearrangecost,inversioncost=inversioncost,_lambda=_lambda,eps=eps,gapopen=gapopen,axis=0)
        return 0,cost,[cost]

    chain.sort(key=lambda s: s[0]) #order by reference position
    qryorder = sorted(xrange(len(chain)), key= lambda i: chain[i][2]) #qry order
    qryorder_inv = sorted(xrange(len(chain)), key=qryorder.__getitem__) #inverse qry order

    lastqstart,lastqend=ctg2range[chain[-1][7]]
    if chain[0][4]==0:
        end=(rlength,rlength,lastqend,lastqend,chain[0][4])
    else:
        end=(rlength,rlength,lastqstart,lastqstart,chain[0][4])

    firstqstart,firstqend=ctg2range[chain[0][7]]
    if chain[0][4]==0:
        start=(0,0,firstqstart,firstqstart,chain[0][4])
    else:
        start=(0,0,firstqend,firstqend,chain[0][4])
    
    #count out of order traversals
    rearrangements=0
    inversions=0

    startcost=gapcost(start,chain[0],rearrangecost=rearrangecost,inversioncost=inversioncost,_lambda=_lambda,eps=eps,gapopen=gapopen,axis=0)

    cost=startcost
    edgecosts=[startcost]

    weight=alfa*chain[0][5]
    
    for ri in xrange(1,len(chain)):

        pblock=chain[ri-1]
        block=chain[ri]

        ps1,pe1,ps2,pe2,po,pscore,pref,pctg=pblock
        s1,e1,s2,e2,o,score,ref2,ctg=block
        weight+=(alfa*score)

        # xgap=0#s1-pe1
        
        pqi=qryorder_inv[ri-1] #index within the qryorder of pblock
        qi=qryorder_inv[ri] #index within the qryorder of block

        if pctg==ctg and pref==ref2:

            if (pqi==qi-1) or (pqi==qi+1): #check if the two blocks are colinear
                gc=gapcost(pblock,block,rearrangecost=rearrangecost,inversioncost=inversioncost,_lambda=_lambda,eps=eps,gapopen=gapopen,axis=0)
                cost+=gc
                edgecosts.append(gc)
            else: #all other options use rearrangement penalty
                rearrangements+=1
                cost+=(gapopen+rearrangecost)
                edgecosts.append(gapopen+rearrangecost)
        
        else: #cross contigs
            
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
                cost+=(gapopen+rearrangecost)
                edgecosts.append((gapopen+rearrangecost))
            else:
                edgecosts.append(gapopen) #simple traversal between two contigs

    endcost=gapcost(chain[-1],end,rearrangecost=rearrangecost,inversioncost=inversioncost,_lambda=_lambda,eps=eps,gapopen=gapopen,axis=0)

    cost+=endcost

    edgecosts.append(endcost)

    return weight,cost,edgecosts

def update_progress(i,n):
    fullbar=100
    if (i+1) % (n/fullbar if n>fullbar else 1)==0 or i+1==n:
        done=int(fullbar*((i+1)/float(n)))
        todo=fullbar-done
        sys.stdout.write('\r[%s%s]'%("#"*done," "*todo))
        if i+1==n:
            sys.stdout.write('\n')
        sys.stdout.flush()

def glocalchain(syntenyblocks, rlength, qlength, ctg2range, rearrangecost=1000, inversioncost=1, lastn=50, lastbp=10000, useheap=False, axis=0, _lambda=5, eps=1, alfa=1, gapopen=10):

    sep=rlength
    
    #add some dummy blocks for the contig start/ends
    if axis==0:
        for refid,(refstart,refend) in enumerate(ctg2range):
            if refstart>=sep:
                break
            if refid==0:
                start=(refstart,refstart,None,None,0,0,None,None)
            syntenyblocks.append((refend,refend,None,None,0,0,None,None))
        end=syntenyblocks[-1]

    if axis==1:
        first=True
        for ctgid,(ctgstart,ctgend) in enumerate(ctg2range):
            if ctgstart<sep:
                continue
            if first:
                start=(None,None,ctgstart,ctgstart,0,0,None,None)
                first=False
            syntenyblocks.append((None,None,ctgend,ctgend,0,0,None,None))
        end=syntenyblocks[-1]

    if axis==0: #sort by ref
        c1,c2=0,2
    else: #sort by qry
        c1,c2=2,0

    syntenyblocks.sort(key=lambda s: (s[c1],-s[5]) ) #order by reference position, then score

    if useheap:
        heap=sortedcontainers.SortedList()
        heap.add((0,start))
    else:
        heap=[(0,start)]+[None]*(len(syntenyblocks))

    G={b:None for b in syntenyblocks}

    maxscore=None

    n=len(syntenyblocks)
    
    bt=range(n+1)

    update_progress(0,n)

    pri=0
    t0=time.time()

    deepest=0
    # best=None

    for ri in xrange(n):
        block=syntenyblocks[ri]

        while syntenyblocks[deepest][c1+1]<block[c1]:
            deepest+=1

        if ri%1000==0:
            t1=time.time()
            sec=t1-t0
            bd=ri-pri
            logging.debug("Blocks per sec: %d"%(bd/sec))
            t0=t1
            pri=ri

        update_progress(ri,n)

        s1,e1,s2,e2,o,score,refid,ctgid=block

        trace=False
        # starttrace=105637436
        # endtrace=starttrace+10
        # if s1>=starttrace and s1<endtrace: # and refid==ctgtrace:
        # #     # if block==(4499237, 4502780, 9008394, 9011937, 0, 3543, 0, 1) or block==end:
        #     logging.info("BLOCK: %s"%str(block))
        #     print "deepest",syntenyblocks[deepest], syntenyblocks[deepest][c1+1]
        #     trace=True

        bestscore=None
        bestblock=None
        bestcost=0
        # checkedbest=False

        l=0
        
        for j in bt: #back track on the heap

            if useheap:
                if j>=len(heap):
                    break
                cscore,pblock=heap[-j]
            else:
                i=(ri+1)-j-1
                if i<0:
                    break
                cscore,pblock=heap[i]
            
            # if best==None or cscore==best:
                # checkedbest=True

            ps1,pe1,ps2,pe2,po,pscore,prefid,pctgid=pblock

            if (pblock[c1]==block[c1] and prefid!=None and refid!=None) or (pblock[c1+1]>=block[c1+1] and prefid!=None and refid!=None):
                continue

            if (pblock[c2]>=block[c2] and prefid!=None and refid!=None) and (pblock[c2+1]<=block[c2+1] and prefid!=None and refid!=None):
                continue            

            l+=1

            if bestscore!=None:
                if cscore<=bestscore:
                    if useheap:
                        break
                    else:
                        if block[c1]-pblock[c1]>lastbp and l>=lastn and pblock[c1]<syntenyblocks[deepest][c1]:
                            break
                        else:
                            continue

            #if block is a dummy block, make it relative to pblock, if possible
            if block[6]==None and pblock[6]!=None: #update current block to be relative to pblock
                if axis==0:
                    _block=(s1, e1, pe2 if po==0 else ps2, pe2 if po==0 else ps2, po, 0, prefid, pctgid)
                else:
                    _block=(pe1 if po==0 else ps1, pe1 if po==0 else ps1, s2, e2, po, 0, prefid, pctgid)
            else:
                _block=block

            #if pblock is a dummy block, make it relative to block, if possible
            if pblock[6]==None and block[6]!=None:
                if axis==0:
                    _pblock=(ps1, pe1, s2 if o==0 else e2, s2 if o==0 else e2, o, 0, refid, ctgid)
                else:
                    _pblock=(s1 if o==0 else e1, s1 if o==0 else e1, ps2, pe2, o, 0, refid, ctgid)
            else:
                _pblock=pblock

            #if blocks come from same query contig and reference contig, compute gapcost, else introduce rearrangement cost
            if _pblock[6]==_block[6]!=None and _pblock[7]==_block[7]!=None:
                c=gapcost(_pblock,_block,rearrangecost=rearrangecost,inversioncost=inversioncost,eps=eps,_lambda=_lambda,gapopen=gapopen,axis=axis)
            elif _pblock[6]==_block[6]==None and _pblock[7]==_block[7]==None: #connect two dummy blocks
                c=gapopen+(abs(block[c1]-(pblock[c1+1]))*eps)
            else: #blocks cross contigs or ref without passing a dummy block, introduce rearrangement cost

                pblockctgstart,pblockctgend=ctg2range[_pblock[7]]
                blockctgstart,blockctgend=ctg2range[_block[7]]
                pblockrefstart,pblockrefend=ctg2range[_pblock[6]]
                blockrefstart,blockrefend=ctg2range[_block[6]]

                if _pblock[6]==_block[6] and axis==0:
                    if _pblock[4]==0:
                        cp=abs( pblockctgend-_pblock[3])
                    else:
                        cp=abs( _pblock[2]-pblockctgstart)

                    if _block[4]==0:
                        cb=abs( blockctgend-_block[3] )
                    else:
                        cb=abs( _block[2]-blockctgstart )

                    c=gapopen+min((rearrangecost,((cp+cb)*eps)))

                elif _pblock[7]==_block[7] and axis==1:

                    if _pblock[4]==0:
                        cp=abs( pblockrefend-_pblock[1])
                    else:
                        cp=abs( _pblock[0]-pblockrefstart)

                    if _block[4]==0:
                        cb=abs( _block[0]-blockrefstart )
                    else:
                        cb=abs( blockrefend-_block[1])

                    c=gapopen+min((rearrangecost,((cp+cb)*eps)))
                else:
                    c=rearrangecost+gapopen+(abs(block[c1]-(pblock[c1+1]))*eps)

            assert(c>=0)

            if trace:
                logging.info("Connect to PBLOCK: %s costs %s, depth=%s, lastbp=%d, cscore,%s, cscore-c=%d, bestscore=%s"%(pblock,c,l,block[c1]-pblock[c1],cscore,cscore-c,bestscore))

            if bestscore==None or cscore-c > bestscore:
                bestscore=cscore-c
                bestblock=pblock
                bestcost=c

            if not useheap:
                if block[c1]-pblock[c1]>lastbp and l>=lastn and pblock[c1]<syntenyblocks[deepest][c1]:
                    break

        # if l>lastn:
            # logging.info("Forced deeper %d backtrack for block: %s"%(l,block))

        cscore=bestscore+(alfa*score)

        # if best==None or cscore>best:
            # best=cscore

        if useheap:
            heap.add((cscore,block))
        else:
            heap[ri+1]=(cscore,block)
        
        if maxscore==None or maxscore<cscore:
            maxscore=cscore
            maxnode=block

        if trace:
            logging.info("CONNECT TO BLOCK: %s, score=%s, cost=%s, depth=%s"%(bestblock,bestscore,bestcost,l))

        G[block]=(bestblock,bestscore)

    node,cscore=G[end]

    chain=[]
    while node!=start:# and node!=startrc:
        if node[6]!=None: #only add to the chain if it is an actual anchor, exclude contig endpoints
            chain.append(node)
        s1,e1,s2,e2,o,score,refid,ctgid=node
        nnode,score=G[node]
        if node==nnode:
            logging.fatal("Loop in chain!")
            sys.exit(1)
        node=nnode

    logging.info("Optimal glocal chain contains: %d anchors and scores %d"%(len(chain),cscore))

    return chain[::-1]

def gapcost(block1,block2,rearrangecost=10000,inversioncost=0,eps=0,_lambda=0.5,gapopen=10,axis=0):

    if axis==0: #sorted by ref
        c1,c2=0,2
    else: #sorted by qry
        c1,c2=2,0

    assert(block1[c1]<=block2[c1])

    d1=block2[c1]-block1[c1+1]

    if block1[4]==block2[4]==0: #both normal orientation
        if block2[c2]<block1[c2]:#always has to be rearranged!
            indelcost=rearrangecost
            substitutioncost=eps*(d1 if d1>0 else 0) #do not penalize if overlap
            return gapopen+indelcost+substitutioncost
        else:
            d2=block2[c2]-block1[c2+1]
            indelcost=min((rearrangecost,_lambda*abs(d1-d2)))
            substitutioncost=eps*max(((d1 if d1<d2 else d2),0))
            return gapopen+indelcost+substitutioncost

    elif block1[4]==block2[4]==1: #both reverse comp orientation
        if block2[c2]>block1[c2]: #always has to be rearranged!
            indelcost=rearrangecost
            substitutioncost=eps*(d1 if d1>0 else 0)
            return gapopen+indelcost+substitutioncost
        else:
            d2=block1[c2]-block2[c2+1]
            indelcost=min((rearrangecost,_lambda*abs(d1-d2)))
            substitutioncost=eps*max(((d1 if d1<d2 else d2),0))
            return gapopen+indelcost+substitutioncost

    elif block1[4]==1 and block2[4]==0:
        d1=max((0,d1))
        if block2[c2]>block1[c2]:
            d2=block2[c2]-block1[c2+1]
            d2=max((0,d2))
            indelcost=min((rearrangecost,_lambda*abs(d1-d2)))
            substitutioncost=eps*max(((d1 if d1<d2 else d2),0))
            return gapopen+indelcost+substitutioncost+inversioncost
        else:
            d2=block1[c2]-block2[c2+1]
            d2=max((0,d2))
            indelcost=min((rearrangecost,_lambda*abs(d1-d2)))
            substitutioncost=eps*max(((d1 if d1<d2 else d2),0))
            return gapopen+indelcost+substitutioncost+inversioncost

    else:
        # assert(block1[4]==0 and block2[4]==1)
        d1=max((0,d1))
        if block2[c2]>block1[c2]:
            d2=block2[c2]-block1[c2+1]
            d2=max((0,d2))
            indelcost=min((rearrangecost,_lambda*abs(d1-d2)))
            substitutioncost=eps*max(((d1 if d1<d2 else d2),0))
            return gapopen+indelcost+substitutioncost+inversioncost
        else:
            d2=block1[c2]-block2[c2+1]
            d2=max((0,d2))
            indelcost=min((rearrangecost,_lambda*abs(d1-d2)))
            substitutioncost=eps*max(((d1 if d1<d2 else d2),0))
            return gapopen+indelcost+substitutioncost+inversioncost

def _gapcost(block1,block2,rearrangecost=10000,inversioncost=0,eps=0,_lambda=0.5,gapopen=10,axis=0):

    if axis==0: #sorted by ref
        c1,c2=0,2
    else: #sorted by qry
        c1,c2=2,0

    assert(block1[c1]<=block2[c1])
    
    d1=abs(block2[c1]-block1[c1+1])
    
    if block1[4]==block2[4]==0: #both normal orientation
        if block2[c2]<block1[c2]:#always has to be rearranged!
            indelcost=rearrangecost
            substitutioncost=eps*d1
            return gapopen+indelcost+substitutioncost
        else:
            d2=abs(block2[c2]-block1[c2+1])
            indelcost=min((rearrangecost,_lambda*abs(d1-d2)))
            substitutioncost=eps*(d1 if d1<d2 else d2)
            return gapopen+indelcost+substitutioncost

    elif block1[4]==block2[4]==1: #both reverse comp orientation
        if block2[c2]>block1[c2]: #always has to be rearranged!
            indelcost=rearrangecost
            substitutioncost=eps*d1
            return gapopen+indelcost+substitutioncost
        else:
            d2=abs(block1[c2]-block2[c2+1])
            indelcost=min((rearrangecost,_lambda*abs(d1-d2)))
            substitutioncost=eps*(d1 if d1<d2 else d2)
            return gapopen+indelcost+substitutioncost

    elif block1[4]==1 and block2[4]==0:
        if block2[c2]>block1[c2]:
            d2=abs(block2[c2]-block1[c2+1])
            indelcost=min((rearrangecost,_lambda*abs(d1-d2)))
            substitutioncost=eps*(d1 if d1<d2 else d2)
            return gapopen+indelcost+substitutioncost+inversioncost
        else:
            d2=abs(block1[c2]-block2[c2+1])
            indelcost=min((rearrangecost,_lambda*abs(d1-d2)))
            substitutioncost=eps*(d1 if d1<d2 else d2)
            return gapopen+indelcost+substitutioncost+inversioncost

    else:
        # assert(block1[4]==0 and block2[4]==1)
        if block2[c2]>block1[c2]:
            d2=abs(block2[c2]-block1[c2+1])
            indelcost=min((rearrangecost,_lambda*abs(d1-d2)))
            substitutioncost=eps*(d1 if d1<d2 else d2)
            return gapopen+indelcost+substitutioncost+inversioncost
        else:
            d2=abs(block1[c2]-block2[c2+1])
            indelcost=min((rearrangecost,_lambda*abs(d1-d2)))
            substitutioncost=eps*(d1 if d1<d2 else d2)
            return gapopen+indelcost+substitutioncost+inversioncost


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
