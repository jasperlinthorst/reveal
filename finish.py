
import reveallib
import reveallib64
from utils import *

def finish(args):
    if args.sa64:
        idx=reveallib64.index(sa=args.sa1, lcp=args.lcp1, cache=args.cache) #enable preconstruction of first SA and LCP array
    else:
        idx=reveallib.index(sa=args.sa1, lcp=args.lcp1, cache=args.cache) #enable preconstruction of first SA and LCP array
    
    G=nx.DiGraph()
    G.graph['samples']=[]
    t=IntervalTree()
    
    reffile=os.path.basename(args.reference)
    ctgfile=os.path.basename(args.contigs)
    
    ref2length=dict()
    idx.addsample(reffile)
    if args.reference.endswith(".gfa"):
        read_gfa(args.reference,idx,t,G)
    else:
        G.graph['samples'].append(reffile)
        for name,seq in fasta_reader(args.reference):
            ref2length[name]=len(seq)
            intv=idx.addsequence(seq)
            intv=Interval(intv[0],intv[1],name)
            t.add(intv)
            G.add_node(intv,offsets={reffile:0})
    
    contig2length=dict()
    contig2seq=dict()
    
    idx.addsample(ctgfile)
    if args.contigs.endswith(".gfa"):
        read_gfa(args.contigs,idx,t,G)
    else:
        G.graph['samples'].append(ctgfile)
        for name,seq in fasta_reader(args.contigs):
            contig2length[name]=len(seq)
            contig2seq[name]=seq
            intv=idx.addsequence(seq)
            intv=Interval(intv[0],intv[1],name)
            t.add(intv)
            G.add_node(intv,offsets={ctgfile:0})
    
    #map nodes to connected components in the graph
    refnode2component=dict()
    ctgnode2component=dict()
    component2refnode=dict()
    component2ctgnode=dict()
    refcomponents=[]
    ctgcomponents=[]

    ctg2ref=dict()
    ri=0
    ci=0
    for nodes in nx.connected_components(G.to_undirected()):
        nodes=list(nodes)
        if reffile in G.node[nodes[0]]['offsets']:
            for node in nodes:
                assert(reffile in G.node[node]['offsets']) #check the graph is valid
                refnode2component[node]=ri
                component2refnode[ri]=node
            ri+=1
            refcomponents.append(nodes)
        else:
            for node in nodes:
                assert(ctgfile in G.node[node]['offsets']) #check the graph is valid
                ctgnode2component[node]=ci
                component2ctgnode[ci]=node
            ci+=1
            ctgcomponents.append(nodes)
    
    idx.construct()
    mems=[]

    for mem in idx.getmems(args.minlength):
        refstart=mem[2][0]
        ctgstart=mem[2][1]
        rnode=t[refstart].pop() #start position on match to node in graph
        cnode=t[ctgstart].pop()
        mems.append((rnode[2], refstart-rnode[0], cnode[2], ctgstart-cnode[0], mem[0], mem[1], mem[3], 0))
    
    logging.info("Indexing reverse complement...\n")
    
    ### index reverse complement
    if args.sa64:
        idx=reveallib64.index(sa=args.sa2, lcp=args.lcp2) #enable preconstruction of second SA and LCP array
    else:
        idx=reveallib.index(sa=args.sa2, lcp=args.lcp2) #enable preconstruction of second SA and LCP array
    
    rcG=nx.DiGraph()
    t=IntervalTree()
    
    idx.addsample(reffile)
    if args.reference.endswith(".gfa"):
        read_gfa(args.reference,idx,t,rcG)
    else:
        rcG.graph['samples']=set([reffile])
        for name,seq in fasta_reader(args.reference):
            intv=idx.addsequence(seq)
            intv=Interval(intv[0],intv[1],name)
            t.add(intv)
            rcG.add_node(intv,offsets={reffile:0},aligned=0)
            refseq=seq
    
    idx.addsample(ctgfile)
    if args.contigs.endswith(".gfa"):
        read_gfa(args.contigs,idx,t,rcG,revcomp=True)
    else:
        rcG.graph['samples']=set([ctgfile])
        for name,seq in fasta_reader(args.contigs):
            intv=idx.addsequence(rc(seq))
            intv=Interval(intv[0],intv[1],name)
            t.add(intv)
            rcG.add_node(intv,offsets={ctgfile:0},aligned=0)
    
    idx.construct()
    
    for mem in idx.getmems(args.minlength):
        refstart=mem[2][0]
        ctgstart=mem[2][1]
        rnode=t[refstart].pop() #start position on match to node in graph
        cnode=t[ctgstart].pop()
        l=cnode[1]-cnode[0]
        mems.append((rnode[2], refstart-rnode[0], cnode[2], l-((ctgstart-cnode[0])+mem[0]), mem[0], mem[1], mem[3], 1))
    
    #determine best scoring chain per contig
    ctg2mums=dict()
    
    logging.debug("Relating exact matches to contigs...")
    for mem in mems:
        refchrom, refstart, ctg, ctgstart, l, n, u, o = mem
        refstart=int(refstart)
        ctgstart=int(ctgstart)
        l=int(l)
        n=int(n)
        u=int(u)
        o=int(o)
        if ctg in ctg2mums:
            if refchrom in ctg2mums[ctg]:
                ctg2mums[ctg][refchrom].append((refstart,ctgstart,l,n,u,o))
            else:
                ctg2mums[ctg][refchrom]=[(refstart,ctgstart,l,n,u,o)]
        else:
            ctg2mums[ctg]=dict({refchrom : [(refstart,ctgstart,l,n,u,o)]})

    logging.info("Done.")
    
    ref2ctg=dict()
    
    logging.info("Mapping contigs to reference.") 
    #for each contig determine the best locations on the reference
    for ctg in ctg2mums:
        logging.debug("Determining best chain for: %s"%ctg)
        paths=[]

        for ref in ctg2mums[ctg]:
            logging.debug("Checking %s"%ref)
            mems=ctg2mums[ctg][ref]
            #determine optimal chain of mems
            path,score=bestmempath(mems,contig2length[ctg],rc=False)
            
            if len(path)>0:
                refstart=path[-1][0]
                refend=path[0][0]+path[0][2]
                ctgstart=path[-1][1]
                ctgend=path[0][1]+path[0][2]
                logging.debug("%s maps to %s, with chain of %d matches, score of %d, reflength %d [%d:%d] and ctglength %d [%d:%d]"%(ctg,ref,len(path),score, refend-refstart,refstart,refend,ctgend-ctgstart,ctgstart,ctgend))
                paths.append((score,ctgstart,ctgend,refstart,refend,ref,False,path))
                 
            #determine optimal chain of mems in reverse complement
            rpath,rscore=bestmempath(mems,contig2length[ctg],rc=True)
            if len(rpath)>0:
                refstart=rpath[-1][0]
                refend=rpath[0][0]+rpath[0][2]
                ctgstart=(rpath[-1][1]+rpath[-1][2])
                ctgend=rpath[0][1]
                logging.debug("Reverse complement of %s maps to %s, with chain of %d matches, score of %d, reflength %d [%d:%d] and ctglength %d [%d:%d]"%(ctg,ref,len(rpath),rscore,refend-refstart,refstart,refend,ctgstart-ctgend,ctgstart,ctgend))
                paths.append((rscore,ctgstart,ctgend,refstart,refend,ref,True,rpath))
        
        if len(paths)==0:
            continue
        
        paths.sort(key=lambda p:p[0],reverse=True) #sort chains by score in descending order, best first
        
        #take the n-best paths that do not overlap
        it=IntervalTree()
        for path in paths:
            score,ctgstart,ctgend,refstart,refend,ref,revcomp,p=path
            if revcomp:
                ctgstart,ctgend=ctgend,ctgstart
            if it[ctgstart:ctgend]==set():
                it[ctgstart:ctgend]=path
            else:
                break
        
        if len(it)>1:
            logging.warn("Contig (len=%d) %s aligns to %d distinct locations on the reference, probably misassembly!"%(contig2length[ctg],ctg[:30],len(it)))
        
        for i,p in enumerate(it):
            start,end,path=p
            score,ctgstart,ctgend,refstart,refend,ref,revcomp,chain=path
            if ref in ref2ctg:
                ref2ctg[ref].append((ctg,revcomp,chain,score,refstart,refend,ctgstart,ctgend,contig2length[ctg]))
            else:
                ref2ctg[ref]=[(ctg,revcomp,chain,score,refstart,refend,ctgstart,ctgend,contig2length[ctg])]
    
    resbase=reffile[:reffile.rfind('.')]+"_"+ctgfile[:ctgfile.rfind('.')]
    
    if not args.split:
        finished=open(resbase+".fasta",'w')
    
    #for each reference chromosome, order the assigned contigs
    for ref in ref2ctg:
        if args.split:
            finished=open(resbase+"_"+ref.replace(" ","").replace("|","").replace("/","").replace(";","").replace(":","")+".fasta",'w')
        
        logging.info("Determining contig order for: %s"%ref)
        
        ref2ctg[ref].sort(key=lambda c: c[4]) #sort by ref start position of chains
        ctgs=ref2ctg[ref]
        
        ctgs=bestctgpath(ctgs)
        
        #TODO: output statistics on number and sequence content of unplaced contigs
        
        frp=0
        if args.plot:
            plt.clf()
            plt.figure(0,figsize=(15,15))
            ax = plt.axes()
            plt.title(ref)
        
        coffset=0
        roffset=0
        yticks=[]
        yticklabels=[]
        
        finished.write(">%s_%s\n"%(ref,ctg))
        
        i=0
        pctg="start"
        
        prefend=0
        pctgend=0
        pctgbegin=0
        prevcomp=False
        
        for ctg,revcomp,path,score,refbegin,refend,ctgbegin,ctgend,ctglength in ctgs:
            
            if revcomp:
                ctgbegin,ctgend=ctgend,ctgbegin
            
            alength=ctgend-ctgbegin
            assert(alength>0)
            
            if args.fix:
                reloffset=ctgbegin
            else:
                reloffset=0
            
            if refend<=frp:
                logging.error("Contained contig should not be in best contig path! %s with alignment length %d"%(ctg,alength))
                break
            else:
                frp=refend
            
            a=refbegin-prefend
            b=coffset-pctgend
            gapsize=a-b
            
            logging.debug("coffset=%d, pctgend=%d, pctgbegin=%d, ctgbegin=%d, ctgend=%d"%(coffset,pctgend,pctgbegin,ctgbegin,ctgend))
            logging.debug("Estimated gap size between %s and %s is %d (%d,%d)."%(pctg,ctg,gapsize,a,b))
            
            if gapsize<0 or args.fixedsize:
                gapsize=100
            
            finished.write("N"*gapsize)
            
            if args.plot:
                #plt.axhline(coffset)
                ax.add_patch(
                    patches.Rectangle(
                        (0, coffset), #bottom left
                        ref2length[ref], #width
                        gapsize, #height
                        alpha=.25
                    )
                )
                
                if args.plotall:
                    for mem in ctg2mums[ctg][ref]:
                        if revcomp:
                            ax.plot([mem[0],mem[0]+mem[2]],[coffset+gapsize+mem[1]+mem[2],coffset+gapsize+mem[1]],'g:')
                        else:
                            ax.plot([mem[0],mem[0]+mem[2]],[coffset+gapsize+mem[1],coffset+gapsize+mem[1]+mem[2]],'r:')
                
                if revcomp:
                    ax.plot([refbegin,refend],[ctgend+gapsize+coffset-reloffset,ctgbegin+gapsize+coffset-reloffset],'bx')
                else:
                    ax.plot([refbegin,refend],[ctgbegin+gapsize+coffset-reloffset,ctgend+gapsize+coffset-reloffset],'bx')
                
                if revcomp:
                    for mem in path:
                        ax.plot([mem[0],mem[0]+mem[2]],[coffset+gapsize+mem[1]+mem[2]-reloffset,coffset+gapsize+mem[1]-reloffset],'g-',linewidth=2)
                else:
                    for mem in path:
                        ax.plot([mem[0],mem[0]+mem[2]],[coffset+gapsize+mem[1]-reloffset,coffset+gapsize+mem[1]+mem[2]-reloffset],'r-',linewidth=2)
            
            pctgend=ctgend-reloffset+coffset+gapsize
            prevcomp=revcomp
            prefend=refend
            
            if args.fix:
                coffset=coffset+gapsize+alength #use the length of the aligned chain
            else:
                coffset=coffset+gapsize+ctglength #use the length of the entire contig
            
            yticks.append(coffset)
            yticklabels.append(ctg[0:15])
            
            logging.info("%d - Contig (revcomp=%d): %s"%(i,revcomp,ctg))
            i+=1
            
            if args.fix:
                if revcomp:
                    finished.write(rc(contig2seq[ctg][ctgbegin:ctgend]))
                else:
                    finished.write(contig2seq[ctg][ctgbegin:ctgend])
            else:
                if revcomp:
                    finished.write(rc(contig2seq[ctg]))
                else:
                    finished.write(contig2seq[ctg])
            
            pctg=ctg
        
        finished.write("\n")
        
        if args.plot:
            ax.set_yticks(yticks)
            ax.set_yticklabels(yticklabels)
            plt.xlim(0,ref2length[ref])
            if args.interactive:
                plt.show()
            else:
                plt.savefig(resbase+"_"+ref.split()[0]+".png")
        
        if args.split:
            finished.close()

def bestctgpath(ctgs,n=10000):
    
    start=(0,0,0,0,0,0,0,0,0)

    link=dict() #key is combination of ctgname and orientation, which should be unique for each chromosome
    score=dict({(0,0):0})

    processed=[]
    active=[start]
    maxscore=0
     
    for ctg in ctgs:
        ctgname,revcomp,path,cscore,refbegin,refend,ctgbegin,ctgend,ctglength=ctg
        
        remove=[]
        for pctg in processed:
            pctgname,prevcomp,ppath,pscore,prefbegin,prefend,pctgbegin,pctgend,pctglength=pctg
            if prefbegin<refbegin:
                active.append(pctg)
                remove.append(pctg)
        
        for r in remove:
            processed.remove(r)
        
        best=None
        w=None
        
        for actg in active:
            #calculate score of connecting to active point
            actgname,arevcomp,apath,ascore,arefbegin,arefend,actgbegin,actgend,actglength=actg
            if w==None:
                w=score[(actg[0],actg[1])]+cscore
                best=actg
            else:
                if arefend>refbegin:
                    penalty=arefend-refbegin #penalize by the amount of overlap
                else:
                    penalty=0
                
                tmpw=score[(actg[0],actg[1])]+cscore-penalty
                if tmpw>w:
                    w=tmpw
                    best=actg
        
        assert(best!=None)

        link[(ctg[0],ctg[1])]=best
        score[(ctg[0],ctg[1])]=w
        
        if w>maxscore:
            maxscore=w
            end=ctg
        
        processed.append(ctg)
    
    #backtrack
    minscore=0
    path=[]
    while end[0]!=start[0]:
        path.append(end)
        end=link[(end[0],end[1])]
    
    return path[::-1]

def bestmempath(mems,ctglength,n=10000,rc=False):
    
    if len(mems)>n: #take only n largest mems
        mems.sort(key=lambda mem: mem[2]) #sort by size
        mems=mems[:n]
    
    mems=[m for m in mems if m[5]==rc]
    
    if len(mems)==0:
        return [],0

    c=sum([m[2] for m in mems])
    logging.debug("Number of anchors: %d",len(mems))
    logging.debug("Sum of anchors: %d", c)
    logging.debug("Length of contig: %d", ctglength)
    logging.debug("Cov ratio: %s"% (c/float(ctglength)) )
    
    mems.sort(key=lambda mem: mem[0]) #sort by reference position
    init=(None, None, 0, 0, 0, 0)
    link=dict()
    score=dict({init:0})
    active=[init]
    processed=[]
    start=init
    end=None 
    maxscore=0
    
    for mem in mems:
        remove=[]
        for pmem in processed:
            pendpoint=pmem[0]+pmem[2]
            if pendpoint<mem[0]:
                active.append(pmem)
                remove.append(pmem)
        
        for r in remove:
            processed.remove(r)
        
        best=None
        w=None
        
        for amem in active:
            #calculate score of connecting to active point
            if rc:
                if w==None:
                    w=score[amem]+mem[2]
                    best=amem
                elif amem[1] >= mem[1]+mem[2]:
                    penalty=abs( (amem[1]-(mem[1]+mem[2]) ) - (mem[0]-(amem[0]+amem[2])  ) )
                    tmpw=score[amem]+mem[2]-penalty
                    if tmpw>w:
                        w=tmpw
                        best=amem
            else:
                if w==None:
                    w=score[amem]+mem[2]
                    best=amem
                elif amem[1] <= mem[1]:
                    penalty=abs( ((amem[1]+amem[2])-mem[1]) - ((amem[0]+amem[2])-mem[0]) )
                    tmpw=score[amem]+mem[2]-penalty
                    if tmpw>w:
                        w=tmpw
                        best=amem
        
        assert(best!=None)
        link[mem]=best
        score[mem]=w
        
        if w>maxscore:
            maxscore=w
            end=mem
        
        processed.append(mem)
    
    #backtrack
    minscore=0
    path=[]
    while end!=start:
        path.append(end)
        if score[end]<minscore:
            minscore=score[end]
            start=end
        end=link[end]
    
    if score[start]<0:
        logging.error("Negative score at start!")
        maxscore-=score[start]
    
    return path,maxscore

