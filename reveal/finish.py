import reveallib
import reveallib64
from utils import *

def inversions(args):
    
    #obtain matches between ref and contigs
    mums,ref2length,contig2length,contig2seq = getmums(args.reference,args.contigs,sa64=args.sa64,minlength=args.minlength)
    logging.debug("MUMS in normal orientation: %d"%len(mums))
    
    #obtain matches between ref and reverse complemented contigs
    rcmums,_,_,_ = getmums(args.reference,args.contigs,revcomp=True,sa64=args.sa64,minlength=args.minlength)
    logging.debug("MUMS in reverse complemented orientation: %d"%len(rcmums))
    
    #combine matches
    mems=mums+rcmums
    
    logging.debug("MUMS in both: %d"%len(mems))
    
    #relate mums to contigs
    ctg2mums=mapmumstocontig(mems)
    
    if args.output==None:
        pref=[]
        for f in [os.path.basename(args.reference),os.path.basename(args.contigs)]:
            bn=os.path.basename(f)
            if '.' in bn:
                pref.append(bn[:bn.find('.')])
            else:
                pref.append(bn)
        args.output="_".join(pref)    
    
    #write finished assembly based on contigs or chains that map on each reference chromosome
    if not args.split:
        finished=open(args.output+".fasta",'w')
    
    excludeflank=True
    
    gap='N'*args.gapsize
    
    for ctg in ctg2mums:
        bestscore=0
        bestpath=[]
        for ref in ctg2mums[ctg]:
            logging.debug("Checking %s"%ref)
            mems=ctg2mums[ctg][ref]
            mems=filtercontainedmums(mems)
            path,score=bestmempathwithinversions(mems)
            if score>bestscore:
                bestpath=path
        
        if len(bestpath)==0:
            logging.warn("No alignment for: %s"%ctg)
        
        if args.plot:
            from matplotlib import pyplot as plt
            plt.clf()
        
        if args.split:
            finished=open(args.output+"_"+ref.replace(" ","").replace("|","").replace("/","").replace(";","").replace(":","")+".fasta",'w')
        
        finished.write(">%s_%s\n"%(ref,os.path.basename(args.contigs)))
        
        path=bestpath
        seq=contig2seq[ctg]
        
        seqoffset=0
        pm=path[0]
        d=pm[4]
        
        for m in path[1:][::-1]:
            if args.plot:
                if pm[4]==0:
                    plt.plot((pm[0],pm[0]+pm[2]),((pm[1],pm[1]+pm[2])),'r-')
                else:
                    plt.plot((pm[0],pm[0]+pm[2]),((pm[1]+pm[2],pm[1])),'g-')
            
            if m[4]!=pm[4]:
                print "Inversion event at:",pm[0]
                if m[4]==1: #to reverse complement
                    if args.plot:
                        plt.axvline(x=pm[0]+pm[2],linewidth=1,color='b',linestyle='dashed')
                    
                    if args.plot:
                        plt.axhline(y=pm[1]+pm[2],linewidth=1,color='b',linestyle='dashed')
                    
                    assert(pm[1]+pm[2]>seqoffset)
                    finished.write(seq[seqoffset:pm[1]+pm[2]]+gap)
                    seqoffset=pm[1]+pm[2]
                    d=1
                else: #from reverse complement
                    if args.plot:
                        plt.axvline(x=pm[0]+pm[2],linewidth=1,color='b',linestyle='dashed')
                    
                    if args.plot:
                        plt.axhline(y=m[1],linewidth=1,color='b',linestyle='dashed')
                    
                    assert(m[1]>seqoffset)
                    finished.write(rc(seq[seqoffset:m[1]])+gap)
                    seqoffset=m[1]
                    d=0
            pm=m
        
        finished.write(seq[seqoffset:])
        finished.write("\n")
        
        if args.split:
            finished.close()
        
        if args.plot:
            plt.show()

def finish(args):
    
    #obtain matches between ref and contigs
    mums,ref2length,contig2length,contig2seq = getmums(args.reference,args.contigs,sa64=args.sa64,minlength=args.minlength)
    logging.debug("MUMS in normal orientation: %d"%len(mums))

    #obtain matches between ref and reverse complemented contigs
    rcmums,_,_,_ = getmums(args.reference,args.contigs,revcomp=True,sa64=args.sa64,minlength=args.minlength)
    logging.debug("MUMS in reverse complemented orientation: %d"%len(rcmums))
    
    #combine matches
    mems=mums+rcmums
    
    logging.debug("MUMS in both: %d"%len(mems))

    #relate mums to contigs
    ctg2mums=mapmumstocontig(mems)
    
    if args.order=='chains':
        ref2ctg=chainstorefence(ctg2mums,contig2length,maxgapsize=args.maxgapsize,minchainlength=args.minchainlength,minchainsum=args.minchainsum)
    else:
        ref2ctg=contigstorefence(ctg2mums,contig2length)
    
    #FOR EACH REF CHROMOSOME WE NOW HAVE A ORDERED LIST OF CHAINS STORED IN THE REF2CTG DICT, NOW JUST WRITE TO FILE...
    
    if args.output==None:
        pref=[]
        for f in [os.path.basename(args.reference),os.path.basename(args.contigs)]:
            bn=os.path.basename(f)
            if '.' in bn:
                pref.append(bn[:bn.find('.')])
            else:
                pref.append(bn)
        args.output="_".join(pref)
    
    #write finished assembly based on contigs or chains that map on each reference chromosome
    if not args.split:
        finished=open(args.output+".fasta",'w')
    
    if args.plot:
        from matplotlib import pyplot as plt
        from matplotlib import patches
    
    #for each reference chromosome, order the assigned chains
    for ref in ref2ctg:
        if args.split:
            finished=open(args.output+"_"+ref.replace(" ","").replace("|","").replace("/","").replace(";","").replace(":","")+".fasta",'w')
        
        logging.info("Determining contig/chain order for: %s"%ref)
        
        ref2ctg[ref].sort(key=lambda c: c[4]) #sort by ref start position of chains
        ctgs=ref2ctg[ref]
        ctgs=bestctgpath(ctgs)
        
        #TODO: output statistics on number and sequence content of unplaced contigs
        
        frp=0
        if args.plot:
            plt.clf()
            #plt.figure(0,figsize=(5,5))
            ax = plt.axes()
            plt.title(ref)
        
        coffset=0
        roffset=0
        yticks=[]
        yticklabels=[]
        
        finished.write(">%s_%s\n"%(ref,os.path.basename(args.contigs)))
        
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
            
            if args.order=='chains':
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
            
            logging.debug("coffset=%d, pctgend=%d, pctgbegin=%d, ctgbegin=%d, ctgend=%d, alength=%d"%(coffset,pctgend,pctgbegin,ctgbegin,ctgend,alength))
            logging.debug("Estimated gap size between %s and %s is %d (%d,%d)."%(pctg[0:10],ctg[0:10],gapsize,a,b))
            
            if gapsize<0 or args.fixedsize:
                gapsize=args.gapsize
            
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
                    for mem in path:
                        ax.plot([mem[0],mem[0]+mem[2]],[coffset+gapsize+mem[1]+mem[2]-reloffset,coffset+gapsize+mem[1]-reloffset],'g-',linewidth=2)
                else:
                    ax.plot([refbegin,refend],[ctgbegin+gapsize+coffset-reloffset,ctgend+gapsize+coffset-reloffset],'bx')
                    for mem in path:
                        ax.plot([mem[0],mem[0]+mem[2]],[coffset+gapsize+mem[1]-reloffset,coffset+gapsize+mem[1]+mem[2]-reloffset],'r-',linewidth=2)
            
            pctgend=ctgend-reloffset+coffset+gapsize
            prevcomp=revcomp
            prefend=refend
            
            if args.order=='chains':
                coffset=coffset+gapsize+alength #use the length of the aligned chain
            else:
                coffset=coffset+gapsize+ctglength #use the length of the entire contig
            
            yticks.append(coffset)
            yticklabels.append(ctg[0:15])
            
            logging.info("%d - Contig (revcomp=%d,refstart=%d,refend=%d,ctgstart=%d,ctgend=%d): %s"%(i,revcomp,refbegin,refend,ctgbegin,ctgend,ctg))
            i+=1
            
            if args.order=='chains':
                if revcomp:
                    finished.write(rc(contig2seq[ctg][ctgbegin:ctgend])) #write only the part of the contig that's part of the chain
                else:
                    finished.write(contig2seq[ctg][ctgbegin:ctgend])
            else:
                if revcomp:
                    finished.write(rc(contig2seq[ctg])) #write the entire contig
                else:
                    finished.write(contig2seq[ctg])
            
            pctg=ctg
        
        finished.write("\n")
        
        if args.split:
            finished.close()
        
        if args.plot:
            ax.set_yticks(yticks)
            ax.set_yticklabels(yticklabels)
            plt.xlim(0,ref2length[ref])
            if args.interactive:
                plt.show()
            else:
                plt.savefig(args.output+"_"+ref.split()[0]+".png")

def chainstorefence(ctg2mums,contig2length,maxgapsize=20000,minchainlength=100,minchainsum=65):
    
    ref2ctg=dict()

    for ctg in ctg2mums:
        logging.debug("Determining best chain(s) for: %s"%ctg)
        paths=[]
        
        for ref in ctg2mums[ctg]:
            logging.debug("Checking %s"%ref)
            mems=ctg2mums[ctg][ref]

            logging.debug("Number of mums for contig: %d"%len(mems))
            mems=filtercontainedmums(mems)
            logging.debug("After filtering for containment: %d"%len(mems))

            candidatepaths=mempaths(mems,contig2length[ctg],revcomp=False,maxgapsize=maxgapsize,minchainlength=minchainlength,minchainsum=minchainsum)
            candidaterpaths=mempaths(mems,contig2length[ctg],revcomp=True,maxgapsize=maxgapsize,minchainlength=minchainlength,minchainsum=minchainsum)
            
            for path,score in candidatepaths:
                if len(path)>0:
                    refstart=path[-1][0]
                    refend=path[0][0]+path[0][2]
                    ctgstart=path[-1][1]
                    ctgend=path[0][1]+path[0][2]
                    #logging.debug("%s maps to %s, with chain of %d matches, score of %d, reflength %d [%d:%d] and ctglength %d [%d:%d]"%(ctg,ref,len(path),score, refend-refstart,refstart,refend,ctgend-ctgstart,ctgstart,ctgend))
                    paths.append((score,ctgstart,ctgend,refstart,refend,ref,False,path))

            for rpath,rscore in candidaterpaths:
                if len(rpath)>0:
                    refstart=rpath[-1][0]
                    refend=rpath[0][0]+rpath[0][2]
                    ctgstart=(rpath[-1][1]+rpath[-1][2])
                    ctgend=rpath[0][1]
                    #logging.debug("Reverse complement of %s maps to %s, with chain of %d matches, score of %d, reflength %d [%d:%d] and ctglength %d [%d:%d]"%(ctg,ref,len(rpath),rscore,refend-refstart,refstart,refend,ctgstart-ctgend,ctgstart,ctgend))
                    paths.append((rscore,ctgstart,ctgend,refstart,refend,ref,True,rpath))
        
        if len(paths)==0:
            continue
        
        nrefchroms=len(set([p[5] for p in paths]))
        
        logging.debug("Found a total of %d chains that map to %d different reference chromosomes."%(len(paths),nrefchroms))
        paths.sort(key=lambda p:p[0],reverse=True) #sort chains by score in descending order, best first
        
        #take the n-best paths that do not overlap on the contig
        it=IntervalTree()
        for path in paths:
            score,ctgstart,ctgend,refstart,refend,ref,revcomp,p=path
            if revcomp:
                ctgstart,ctgend=ctgend,ctgstart
            assert(ctgstart<ctgend)
            if it[ctgstart:ctgend]==set():
                it[ctgstart:ctgend]=path
        
        if len(it)>1:
            logging.warn("Contig (len=%d) aligns to %d distinct locations or orientations on the reference."%(contig2length[ctg],len(it)))
            for path in it:
                score,ctgstart,ctgend,refstart,refend,ref,revcomp,p=path[2]
                logging.debug("Actually using path from ctg:%d:%d to ref:%d:%d with orientation %d with score %d"%(ctgstart,ctgend,refstart,refend,revcomp,score))
        
        for i,p in enumerate(it):
            start,end,path=p
            score,ctgstart,ctgend,refstart,refend,ref,revcomp,chain=path
            if ref in ref2ctg:
                ref2ctg[ref].append((ctg,revcomp,chain,score,refstart,refend,ctgstart,ctgend,contig2length[ctg]))
            else:
                ref2ctg[ref]=[(ctg,revcomp,chain,score,refstart,refend,ctgstart,ctgend,contig2length[ctg])]
    
    return ref2ctg


def contigstorefence(ctg2mums,contig2length):
    
    ref2ctg=dict()

    for ctg in ctg2mums:
        logging.debug("Determining best chain(s) for: %s"%ctg)
        paths=[]
        
        for ref in ctg2mums[ctg]:
            logging.debug("Checking %s"%ref)
            mems=ctg2mums[ctg][ref]
            mems=filtercontainedmums(mems)
            
            path,score=bestmempath(mems,contig2length[ctg],revcomp=False)
            rpath,rscore=bestmempath(mems,contig2length[ctg],revcomp=True)
            
            if score>rscore:
                refstart=path[-1][0]
                refend=path[0][0]+path[0][2]
                ctgstart=path[-1][1]
                ctgend=path[0][1]+path[0][2]
                logging.debug("%s maps to %s, with chain of %d matches, score of %d, reflength %d [%d:%d] and ctglength %d [%d:%d]"%(ctg,ref,len(path),score, refend-refstart,refstart,refend,ctgend-ctgstart,ctgstart,ctgend))
                paths.append((score,ctgstart,ctgend,refstart,refend,ref,False,path))
            elif rscore>0:
                refstart=rpath[-1][0]
                refend=rpath[0][0]+rpath[0][2]
                ctgstart=(rpath[-1][1]+rpath[-1][2])
                ctgend=rpath[0][1]
                logging.debug("Reverse complement of %s maps to %s, with chain of %d matches, score of %d, reflength %d [%d:%d] and ctglength %d [%d:%d]"%(ctg,ref,len(rpath),rscore,refend-refstart,refstart,refend,ctgstart-ctgend,ctgstart,ctgend))
                paths.append((rscore,ctgstart,ctgend,refstart,refend,ref,True,rpath))
        
        if len(paths)==0:
            continue
        
        logging.debug("Found a total of %d chains that map to $d different reference chromosomes.")
        paths.sort(key=lambda p:p[0],reverse=True) #sort chains by score in descending order, best first
        
        #take the n-best paths that do not overlap on the contig
        it=IntervalTree()
        for path in paths:
            score,ctgstart,ctgend,refstart,refend,ref,revcomp,p=path
            if revcomp:
                ctgstart,ctgend=ctgend,ctgstart
            if it[ctgstart:ctgend]==set():
                it[ctgstart:ctgend]=path
        
        if len(it)>1:
            logging.warn("Contig (len=%d) aligns to %d distinct locations or orientations on the reference."%(contig2length[ctg],len(it)))
        
        for i,p in enumerate(it):
            start,end,path=p
            score,ctgstart,ctgend,refstart,refend,ref,revcomp,chain=path
            if ref in ref2ctg:
                ref2ctg[ref].append((ctg,revcomp,chain,score,refstart,refend,ctgstart,ctgend,contig2length[ctg]))
            else:
                ref2ctg[ref]=[(ctg,revcomp,chain,score,refstart,refend,ctgstart,ctgend,contig2length[ctg])]
    
    return ref2ctg


def mapmumstocontig(mems):
    ctg2mums=dict()
    for mem in mems:
        refchrom, refstart, ctg, ctgstart, l, n, o = mem
        refstart=int(refstart)
        ctgstart=int(ctgstart)
        l=int(l)
        n=int(n)
        o=int(o)
        if ctg in ctg2mums:
            if refchrom in ctg2mums[ctg]:
                ctg2mums[ctg][refchrom].append((refstart,ctgstart,l,n,o))
            else:
                ctg2mums[ctg][refchrom]=[(refstart,ctgstart,l,n,o)]
        else:
            ctg2mums[ctg]=dict({refchrom : [(refstart,ctgstart,l,n,o)]})
    return ctg2mums

def getmums(reference, query, revcomp=False, sa64=False, minlength=20):

    if sa64:
        idx=reveallib64.index()
    else:
        idx=reveallib.index()
    
    G=nx.DiGraph()
    G.graph['samples']=[]
    t=IntervalTree()
    
    reffile=os.path.basename(reference)
    ctgfile=os.path.basename(query)
    
    ref2length=dict()
    idx.addsample(reffile)
    if reference.endswith(".gfa"):
        read_gfa(reference,idx,t,G)
    else:
        G.graph['samples'].append(reffile)
        for name,seq in fasta_reader(reference):
            ref2length[name]=len(seq)
            intv=idx.addsequence(seq)
            intv=Interval(intv[0],intv[1],name)
            t.add(intv)
            G.add_node(intv,offsets={reffile:0})
    
    contig2length=dict()
    contig2seq=dict()
    
    idx.addsample(ctgfile)
    if query.endswith(".gfa"):
        read_gfa(query,idx,t,G,revcomp=revcomp)
    else:
        G.graph['samples'].append(ctgfile)
        for name,seq in fasta_reader(query):
            contig2length[name]=len(seq)

            if revcomp:
                rcseq=rc(seq)
                contig2seq[name]=rcseq
                intv=idx.addsequence(rcseq)
            else:
                contig2seq[name]=seq
                intv=idx.addsequence(seq)

            intv=Interval(intv[0],intv[1],name)
            t.add(intv)
            G.add_node(intv,offsets={ctgfile:0})
    
    idx.construct()
    
    mums=[]
    
    for mum in idx.getmums(minlength):
        refstart=mum[2][0]
        ctgstart=mum[2][1]
        rnode=t[refstart].pop() #start position on match to node in graph
        cnode=t[ctgstart].pop()
        if revcomp:
            l=cnode[1]-cnode[0]
            mums.append((rnode[2], refstart-rnode[0], cnode[2], l-((ctgstart-cnode[0])+mum[0]), mum[0], mum[1], 1))
        else:
            mums.append((rnode[2], refstart-rnode[0], cnode[2], ctgstart-cnode[0], mum[0], mum[1], 0))
    
    return mums,ref2length,contig2length,contig2seq

def filtercontainedmums(mems):
    #filter mems before trying to form chains
    mems.sort(key=lambda m: m[0]) #sort by reference position
    filteredmems=[]
    prefstart=0
    prefend=0
    for mem in mems:
        refstart=mem[0]
        refend=mem[0]+mem[2]
        if refend<prefend and prefstart<refstart: #ref contained
            continue
        else:
            filteredmems.append(mem)
            prefend=refend
            prefstart=prefstart
    
    filteredmems.sort(key=lambda m: m[1]) #sort by qry position
    mems=[]
    pqrystart=0
    pqryend=0
    for mem in filteredmems:
        qrystart=mem[1]
        qryend=mem[1]+mem[2]
        if qryend<pqryend and pqrystart<qrystart: #qry contained
            continue
        else:
            mems.append(mem)
            pqryend=qryend
            pqrystart=pqrystart
    return mems

def bestctgpath(ctgs,n=10000):
    
    #start=(0,0,0,0,0,0,0,0,0)
    start=(0,0,[],0,0,0,0,0,0)

    link=dict() #key is combination of ctgname and orientation, which should be unique for each chromosome
    score=dict({():0})
    
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
                w=score[tuple(actg[2])]+cscore
                best=actg
            else:
                if arefend>refbegin:
                    penalty=arefend-refbegin #penalize by the amount of overlap
                else:
                    penalty=0
                
                tmpw=score[tuple(actg[2])]+cscore-penalty
                if tmpw>w:
                    w=tmpw
                    best=actg
            
        assert(best!=None)

        link[tuple(path)]=best
        score[tuple(path)]=w
        
        if w>maxscore:
            maxscore=w
            end=ctg
        
        processed.append(ctg)
    
    #backtrack
    minscore=0
    path=[]
    while end[0]!=start[0]:
        path.append(end)
        end=link[tuple(end[2])]
    
    return path[::-1]


def bestmempath(mems,ctglength,n=10000,revcomp=False):
    
    if len(mems)>n: #take only n largest mems
        mems.sort(key=lambda mem: mem[2]) #sort by size
        mems=mems[:n]
    
    mems=[m for m in mems if m[4]==revcomp]
    
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
    
    active=[]
    
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
        
        best=init
        w=0
        
        #sort active by score decreasing
        active.sort(key=lambda x: score[x], reverse=True)
        
        for amem in active:
            s=score[amem]
            
            assert(s>=0)
            if w > s: #as input is sorted by score
                break
            
            #calculate score of connecting to active point
            if revcomp:
                if amem[1] >= mem[1]+mem[2]:
                    penalty=abs( (amem[1]-(mem[1]+mem[2]) ) - (mem[0]-(amem[0]+amem[2])  ) )
                    tmpw=score[amem]-penalty
                    if tmpw>w:
                        w=tmpw
                        best=amem
            else:
                if amem[1] <= mem[1]:
                    penalty=abs( ((amem[1]+amem[2])-mem[1]) - ((amem[0]+amem[2])-mem[0]) )
                    tmpw=score[amem]+mem[2]-penalty
                    if tmpw>w:
                        w=tmpw
                        best=amem
        
        assert(best!=None)
        
        link[mem]=best
        score[mem]=w+mem[2]
        
        if w+mem[2]>maxscore:
            maxscore=w+mem[2]
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
    
    #if score[start]<0:
    #    logging.error("Negative score at start!")
    #    maxscore-=score[start]
    
    return path,maxscore


def mempaths(mems,ctglength,n=10000,revcomp=False,maxgapsize=20000,minchainlength=100,minchainsum=65):
    
    if len(mems)>n: #take only n largest mems
        mems.sort(key=lambda mem: mem[2]) #sort by size
        mems=mems[:n]
    
    mems=[m for m in mems if m[4]==revcomp]
    
    if len(mems)==0:
        return []
    
    c=sum([m[2] for m in mems])
    logging.debug("Number of anchors: %d",len(mems))
    logging.debug("Sum of anchors: %d", c)
    logging.debug("Length of contig: %d", ctglength)
    logging.debug("Cov ratio: %s"% (c/float(ctglength)) )
    
    mems.sort(key=lambda mem: mem[0]) #sort by reference position
    paths=[]
    
    #extract the best path, remove mems that are part of or start/end within the range of the best path, until no more mems remain
    while len(mems)>0:
        init=(None, None, 0, 0, 0, 0)
        link=dict()
        score=dict({init:0})
        active=[]
        processed=[]
        start=init
        end=None
        maxscore=0
        
        for mem in mems:
            remove=[]
            for pmem in processed:
                pendpoint=pmem[0]+pmem[2]
                if pendpoint<mem[0]+mem[2]:
                    active.append(pmem)
                    remove.append(pmem)
            
            for r in remove:
                processed.remove(r)
            
            best=init
            w=mem[2]
            
            if revcomp:
                subactive=[amem for amem in active if amem[0]+amem[2]>mem[0]-maxgapsize and amem[1]<mem[1]+mem[2]+maxgapsize]
            else:
                subactive=[amem for amem in active if amem[0]+amem[2]>mem[0]-maxgapsize and amem[1]+amem[2]>mem[1]-maxgapsize]
            
            for amem in subactive:
                #calculate score of connecting to active point
                if revcomp:
                    if amem[1] >= mem[1]:
                        p1=(mem[0], mem[1]+mem[2])
                        p2=(amem[0]+amem[2], amem[1])
                        penalty=sumofpairs(p1,p2)
                        tmpw=score[amem]+mem[2]-penalty
                        if tmpw>w:
                            w=tmpw
                            best=amem
                else:
                    if amem[1]+amem[2] <= mem[1]+mem[2]:
                        p1=(amem[0]+amem[2], amem[1]+amem[2])
                        p2=(mem[0], mem[1])
                        penalty=sumofpairs(p1,p2)
                        tmpw=score[amem]+mem[2]-penalty
                        if tmpw>w:
                            w=tmpw
                            best=amem
            
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
        
        if len(path)!=0:
            refstart=path[-1][0]
            refend=path[0][0]+path[0][2]
            
            if revcomp:
                ctgstart=path[-1][1]+path[0][2]
                ctgend=path[0][1]
            else:
                ctgstart=path[-1][1]
                ctgend=path[0][1]+path[0][2]
            
            assert(refend>refstart)
            
            if refend-refstart>minchainlength:
                if sum([m[2] for m in path])>minchainsum:
                    paths.append((path,maxscore))
                    logging.debug("Added path from ctg:%d:%d to ref:%d:%d with orientation %d with score %d."%(refstart,refend,ctgstart,ctgend,revcomp,maxscore))
            
            if revcomp: #remove all mems that were part of or overlap the chain that was just found
                assert(ctgstart>ctgend)
                mems=[mem for mem in mems if (mem[0]+mem[2]<refstart or mem[0]>refend) and (mem[1]>ctgend or mem[1]+mem[2]<ctgstart)]
            else:
                mems=[mem for mem in mems if (mem[0]+mem[2]<refstart or mem[0]>refend) and (mem[1]+mem[2]<ctgstart or mem[1]>ctgend)]
        else:
            break
    
    logging.debug("Detected number of chains: %d."%len(paths))
    return paths


def bestmempathwithinversions(mems,n=10000):
    
    if len(mems)>n: #take only n largest mems
        mems.sort(key=lambda mem: mem[2]) #sort by size
        mems=mems[:n]
    
    if len(mems)==0:
        return [],0
    
    c=sum([m[2] for m in mems])
    logging.debug("Number of anchors: %d",len(mems))
    logging.debug("Sum of anchors: %d", c)
    
    mems.sort(key=lambda mem: mem[0]) #sort by reference position
    
    init=(None, None, 0, 0, 0, 0)
    link=dict()
    score=dict({init:0})
    active=[]
    processed=[]
    start=init
    end=None
    maxscore=0
    
    trace=False
    for mem in mems:

        remove=[]
        for pmem in processed:
            pendpoint=pmem[0]+pmem[2]
            if pendpoint<mem[0]+mem[2]:
                active.append(pmem)
                remove.append(pmem)
        
        for r in remove:
            processed.remove(r)
        
        best=init
        w=mem[2]
        
        for amem in active:
            #calculate score of connecting to active point
            tmpw=0
            
            if mem[4]==1: #reverse complement match, connected to top
                if amem[4]==1: #active mem also from reverse complement
                    if amem[1] >= mem[1]: #may overlap, but mem may not be contained
                        p1=(mem[0], mem[1]+mem[2])
                        p2=(amem[0]+amem[2], amem[1])
                        penalty=sumofpairs(p1,p2)
                        tmpw=score[amem]+mem[2]-penalty
                else: #active mem has regular orientation mem from reverse complement
                    if mem[1]+mem[2] >= amem[1]+amem[2]:
                        p1=tuple([mem[0]])#, mem[1]+mem[2])
                        p2=tuple([amem[0]])#, amem[1]+amem[2])
                        penalty=sumofpairs(p1,p2)   #need different penalty here! *********
                        tmpw=score[amem]+mem[2]-penalty
            else: #mem on regular orientation
                if amem[4]==1: #active mem from reverse complement
                    if amem[1] <= mem[1]: #then it has to be 'above' mem
                        p1=tuple([mem[0]])#, mem[1])
                        p2=tuple([amem[0]+amem[2]])#, amem[1])
                        penalty=sumofpairs(p1,p2)   #need different penalty here! ********
                        tmpw=score[amem]+mem[2]-penalty
                else: #mem and amem on same regular orientation
                    if amem[1]+amem[2] <= mem[1]+mem[2]: #may overlap, but mem may not be contained
                        p1=(mem[0], mem[1])
                        p2=(amem[0]+amem[2], amem[1]+amem[2])
                        penalty=sumofpairs(p1,p2)
                        tmpw=score[amem]+mem[2]-penalty
            
            if tmpw>w:
                w=tmpw
                best=amem
            
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
    
    #from matplotlib import pyplot as plt
    
    #for mem in path:
    #    if mem[4]==1:
    #        plt.plot((mem[0],mem[0]+mem[2]),(mem[1]+mem[2], mem[1]),'g-')
    #    else:
    #        plt.plot((mem[0],mem[0]+mem[2]),(mem[1], mem[1]+mem[2]),'r-')
    #plt.show()
    
    return path,maxscore

