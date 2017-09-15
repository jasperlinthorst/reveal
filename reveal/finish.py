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
    ctg2mums=mapmumstocontig(mems,filtercontained=args.filtercontained)
    
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
                    
                    assert(m[1]>=seqoffset)
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
    
    #from multiprocessing import Pool
    #p=Pool(2)
    #p.apply_async(getmums,(args.reference,args.contigs),dict(sa64=args.sa64,minlength=args.minlength))
    
    logging.debug("Extracting mums in normal orientation.")
    #obtain matches between ref and contigs
    mums,ref2length,contig2length,contig2seq = getmums(args.reference,args.contigs,sa64=args.sa64,minlength=args.minlength)
    logging.debug("MUMS in normal orientation: %d"%len(mums))
    
    logging.debug("Extracting mums in reverse complemented normal orientation.")
    #obtain matches between ref and reverse complemented contigs
    rcmums,_,_,_ = getmums(args.reference,args.contigs,revcomp=True,sa64=args.sa64,minlength=args.minlength)
    logging.debug("MUMS in reverse complemented orientation: %d"%len(rcmums))
    
    #combine matches
    mems=mums+rcmums

    assert(len(mems)==len(mums)+len(rcmums))
    
    logging.debug("Associating mums to contigs.")
    #relate mums to contigs
    ctg2mums=mapmumstocontig(mems,filtercontained=args.filtercontained,maxgapsize=args.maxgapsize,minchainlength=args.minchainlength)
    
    if args.order=='chains':
        ref2ctg=chainstorefence(ctg2mums,contig2length,maxn=args.maxn,maxgapsize=args.maxgapsize,minchainlength=args.minchainlength,minchainsum=args.minchainsum)
    else:
        ref2ctg=contigstorefence(ctg2mums,contig2length,maxn=args.maxn)
    
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
        unplaced=open(args.output+".unplaced.fasta",'w')
    
    if args.plot:
        from matplotlib import pyplot as plt
        from matplotlib import patches
    
    #for each reference chromosome, order the assigned chains
    for ref in ref2ctg:
        if args.split:
            finished=open(args.output+"_"+ref.replace(" ","").replace("|","").replace("/","").replace(";","").replace(":","")+".fasta",'w')
            unplaced=open(args.output+"_"+ref.replace(" ","").replace("|","").replace("/","").replace(";","").replace(":","")+"unplaced.fasta",'w')
        
        if ref=='unplaced':
            logging.info("The following %d contigs could not be placed anywhere on the reference sequence.")
            for name in ref2ctg[ref]:
                logging.info("%s (length=%d)"%(name,contig2length[name]))
                unplaced.write(">%s\n"%name)
                unplaced.write("%s\n"%contig2seq[name])
            continue
        
        logging.info("Determining %s order for: %s"%(args.order,ref))
        
        ref2ctg[ref].sort(key=lambda c: c[4]) #sort by ref start position

        ctgs=ref2ctg[ref]
        nctgsin=len(ctgs)

        b=set([(c[0],c[8]) for c in ctgs])
        
        ctgs=bestctgpath(ctgs)
        
        a=set([(c[0],c[8]) for c in ctgs])
       
        logging.debug("Selected %d out of %d %s to layout assembly with respect to %s."%(len(ctgs),nctgsin,args.order,ref))

        if args.order=='contigs':
            if len(b)-len(a)>0:
                logging.info("The following %d contigs were placed on reference sequence %s but were not used:"%(len(b)-len(a),ref))
                for name,l in b - a:
                    logging.info("%s (length=%d)"%(name,l))
                    unplaced.write(">%s\n"%name)
                    unplaced.write("%s\n"%contig2seq[name])
        
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
        o=0
        pctg=("start",ctgs[0][1],[],0,0,0,0,0,0,0)
        
        for ctg in ctgs:
            
            ctgname,revcomp,path,score,refbegin,refend,ctgbegin,ctgend,ctglength,ci=ctg
            
            if revcomp:
                ctgbegin,ctgend=ctgend,ctgbegin
            
            pctgname,prevcomp,ppath,pscore,prefbegin,prefend,pctgbegin,pctgend,pctglength,pci=pctg
            
            if prevcomp:
                pctgbegin,pctgend=pctgend,pctgbegin
            
            alength=ctgend-ctgbegin
            
            assert(alength>0)
            
            if args.order=='chains':
                reloffset=ctgbegin
            else:
                reloffset=0
            
            if refend<=prefend:
                logging.error("Contained contig should not be in best contig path! %s with alignment length %d"%(ctgname,alength))
                break
            
            gapsize=refbegin-prefend
            
            if gapsize==0: #perfect boundary, stil use one N to be able to distinguish the event
                gapsize=1
            
            logging.debug("coffset=%d, pctgend=%d, pctgbegin=%d, ctgbegin=%d, ctgend=%d, alength=%d"%(coffset,pctgend,pctgbegin,ctgbegin,ctgend,alength))
            #logging.info("Estimated gap size between %s and %s is %d (%d,%d)."%(pctgname[0:10],ctgname[0:10],gapsize,a,b))
            
            if gapsize<0 or args.fixedsize:
                gapsize=args.gapsize
            
            #pctgend_=ctgend-reloffset+coffset+gapsize
            
            logging.debug("%d (index on ctg: %d->%d) - Order %s (revcomp=%d,refstart=%d,refend=%d,ctgstart=%d,ctgend=%d)"%(i,pci,ci,args.order,revcomp,refbegin,refend,ctgbegin,ctgend))
            
            if args.order=='chains':
                
                if pctg[0]!="start":
                    #if pctg was consecutive on ref and ctg and with same orientation, write seq in between instead of gap
                    
                    if pctgname==ctgname and revcomp==prevcomp and ((pci==ci+1 and revcomp==1) or (pci==ci-1 and revcomp==0)):
                        logging.debug("Consecutive chains, no gap.")
                        
                        if (pctgend>ctgbegin and revcomp==0):
                            ctgbegin=pctgend #they overlap in normal order, so no sequence between chains
                        elif (ctgend>pctgbegin and revcomp==1):
                            ctgend=pctgbegin
                        else:
                            if revcomp:
                                assert(ctgend<=pctgbegin)
                                finished.write(rc(contig2seq[ctgname][ctgend:pctgbegin]))
                                gapsize=pctgbegin-ctgend
                            else:
                                assert(pctgend<=ctgbegin)
                                finished.write(contig2seq[ctgname][pctgend:ctgbegin])
                                gapsize=ctgbegin-pctgend
                        
                        gap=False
                    else:
                        if i==0:
                            event='start'
                        elif pctgname!=ctgname:
                            event='cross-contig'
                        elif revcomp!=prevcomp:
                            event='inversion'                        
                        elif ((pci!=ci+1 and revcomp==1) or (pci!=ci-1 and revcomp==0)):
                            event='translocation'
                        else:
                            event='error!'
                        
                        logging.info("Event of type: \'%s\' detected at: %d inserting gap of size: %d"%(event,prefend,gapsize))
                        logging.debug("Chains are not consecutive on contig: %s,%s - %s,%s - %s,%s"%(revcomp,prevcomp,pctgname,ctgname,pci,ci-1))
                        
                        gap=True
                        finished.write("N"*gapsize)
                else:
                    gap=False
                
                l=gapsize+alength
                if revcomp:
                    finished.write(rc(contig2seq[ctgname][ctgbegin:ctgend])) #write only the part of the contig that's part of the chain
                else:
                    finished.write(contig2seq[ctgname][ctgbegin:ctgend])
                
            else:
                gap=True
                finished.write("N"*gapsize)
                l=gapsize+len(contig2seq[ctgname])
                
                if revcomp:
                    finished.write(rc(contig2seq[ctgname])) #write the entire contig
                else:
                    finished.write(contig2seq[ctgname])
            
            if args.plot:
                
                if gap:
                    ax.add_patch(
                        patches.Rectangle(
                            (0, o), #bottom left
                            ref2length[ref], #width
                            gapsize, #height
                            alpha=.25
                        )
                    )
                
                if revcomp:
                    #ax.plot([refbegin,refend],[o+alength+gapsize,o+gapsize],'bx')
                    ax.plot([refbegin,refend],[o+gapsize,o+alength+gapsize],'bx')
                    #ax.plot([refbegin,refend],[o+alength+gapsize,o+gapsize],'g-')
                    ax.plot([refbegin,refend],[o+gapsize,o+alength+gapsize],'g-')
                else:
                    ax.plot([refbegin,refend],[o+gapsize,o+alength+gapsize],'bx')
                    ax.plot([refbegin,refend],[o+gapsize,o+alength+gapsize],'r-')
            i+=1
            o=o+l
            yticks.append(o)
            yticklabels.append("%s:%d"%(ctgname[0:15],ctgend))
            pctg=ctg
        
        finished.write("\n")
        
        if args.split:
            finished.close()
            unplaced.close()
        
        logging.info("Done.")
        
        if args.plot:
            ax.set_yticks(yticks)
            ax.set_yticklabels(yticklabels)
            plt.xlim(0,ref2length[ref])
            if args.interactive:
                plt.show()
            else:
                plt.savefig(args.output+"_"+ref.split()[0]+".png")
    
    if not args.split:
        finished.close()
        unplaced.close()

def chainstorefence(ctg2mums,contig2length,maxgapsize=1500,minchainlength=1500,minchainsum=65,maxn=15000):
    
    ref2ctg=dict()
    
    for ctg in ctg2mums:
        logging.debug("Determining best chain(s) for: %s"%ctg)
        paths=[]
        
        for ref in ctg2mums[ctg]:
            logging.debug("Checking %s"%ref)
            mems=ctg2mums[ctg][ref]
            candidatepaths=mempathsbothdirections(mems,contig2length[ctg],n=maxn,maxgapsize=maxgapsize,minchainlength=minchainlength,minchainsum=minchainsum)
            
            for path,score,rc,ctgstart,ctgend,refstart,refend in candidatepaths:
                if len(path)>0:
                    paths.append((score,ctgstart,ctgend,refstart,refend,ref,rc,path))
        
        if len(paths)==0:
            continue
        
        nrefchroms=len(set([p[5] for p in paths]))
        
        logging.debug("Found a total of %d chains that map to %d different reference chromosomes."%(len(paths),nrefchroms))
        
        #sort by score
        paths=sorted(paths,key=lambda c: c[0])
        
        selectedpaths=[]
        #take n-best paths that dont overlap
        it=IntervalTree()
        for path in paths:
            score,ctgstart,ctgend,refstart,refend,ref,revcomp,p=path
            
            logging.debug("Path: ctg:%d:%d - ref:%d:%d (%d) with score %d"%(ctgstart,ctgend,refstart,refend,revcomp,score))
            
            if revcomp:
                ctgend,ctgstart=ctgstart,ctgend
            
            s=it[ctgstart:ctgend]
            
            if s==set():
                it[ctgstart:ctgend]=path
                selectedpaths.append(path)
            else:
                for start,end,v in s:
                    if start<=ctgstart and end>=ctgend: #allow overlap, not containment
                        break
                else:
                    it[ctgstart:ctgend]=path
                    selectedpaths.append(path)
        
        paths=sorted(selectedpaths,key=lambda c: c[1] if c[6] else c[2])
        
        for i,path in enumerate(paths):
            score,ctgstart,ctgend,refstart,refend,ref,revcomp,chain=path
            #logging.info("%d -> ctg:%d:%d - ref:%d:%d - %d"%(i,ctgstart,ctgend,refstart,refend,revcomp))

            if ref in ref2ctg:
                ref2ctg[ref].append((ctg,revcomp,chain,score,refstart,refend,ctgstart,ctgend,contig2length[ctg],i))
            else:
                ref2ctg[ref]=[(ctg,revcomp,chain,score,refstart,refend,ctgstart,ctgend,contig2length[ctg],i)]

    return ref2ctg


def contigstorefence(ctg2mums,contig2length,maxn=15000):
    
    ref2ctg=dict()

    for ctg in ctg2mums:
        logging.debug("Determining best chain for: %s"%ctg)
        paths=[]
        
        for ref in ctg2mums[ctg]:
            logging.debug("Checking %s"%ref)
            mems=ctg2mums[ctg][ref]
            path,score=bestmempath(mems,contig2length[ctg],n=maxn,revcomp=False)
            rpath,rscore=bestmempath(mems,contig2length[ctg],n=maxn,revcomp=True)
            
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
            if 'unplaced' in ref2ctg:
                ref2ctg['unplaced'].append(ctg)
            else:
                ref2ctg['unplaced']=[ctg]
            continue
        
        paths.sort(key=lambda p:p[0],reverse=True) #sort chains by score in descending order, best first
        
        score,ctgstart,ctgend,refstart,refend,ref,revcomp,chain=paths[0] #just take the best path
        if ref in ref2ctg:
            ref2ctg[ref].append((ctg,revcomp,chain,score,refstart,refend,ctgstart,ctgend,contig2length[ctg],0))
        else:
            ref2ctg[ref]=[(ctg,revcomp,chain,score,refstart,refend,ctgstart,ctgend,contig2length[ctg],0)]
    
    return ref2ctg

def mapmumstocontig(mems,filtercontained=True,maxgapsize=1500,minchainlength=1500):
    ctg2mums=dict()
    logging.debug("Mapping %d mums to contigs."%len(mems))
    for mem in mems:
        refchrom, refstart, ctg, ctgstart, l, n, o = mem
        refstart=int(refstart)
        ctgstart=int(ctgstart)
        l=int(l)
        n=int(n)
        o=int(o)
        if ctg in ctg2mums:
            if refchrom in ctg2mums[ctg]:
                ctg2mums[ctg][refchrom].append([refstart,ctgstart,l,n,o])
            else:
                ctg2mums[ctg][refchrom]=[[refstart,ctgstart,l,n,o]]
        else:
            ctg2mums[ctg]=dict({refchrom : [[refstart,ctgstart,l,n,o]]})
    
    if filtercontained:
        for ctg in ctg2mums:
            for ref in ctg2mums[ctg]:
                #ctg2mums[ctg][ref]=filtercontainedmums(ctg2mums[ctg][ref])
                #ctg2mums[ctg][ref]=filtercontainedmumsboth(ctg2mums[ctg][ref])
                #ctg2mums[ctg][ref]=filtercontainedmumsratio(ctg2mums[ctg][ref])
                logging.debug("Filtering mums between %s and %s"%(ref,ctg))
                ctg2mums[ctg][ref]=filtermumsrange(ctg2mums[ctg][ref],maxgapsize=maxgapsize,minchainlength=minchainlength)
    
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
            #mums.append([rnode[2], refstart-rnode[0], cnode[2], l-((ctgstart-cnode[0])+mum[0]), mum[0], mum[1], 1])
        else:
            mums.append((rnode[2], refstart-rnode[0], cnode[2], ctgstart-cnode[0], mum[0], mum[1], 0))
            #mums.append([rnode[2], refstart-rnode[0], cnode[2], ctgstart-cnode[0], mum[0], mum[1], 0])
    
    return mums,ref2length,contig2length,contig2seq


def filtercontainedmumsboth(mems):
    #filter mems before trying to form chains
    before=len(mems)
    mems.sort(key=lambda m: m[0]) #sort by reference position
    filteredmems=[]
    prefstart=0
    prefend=0
    for mem in mems:
        refstart=mem[0]
        refend=mem[0]+mem[2]
        if refend<prefend and prefstart<refstart: #ref contained
            filteredmems.append(mem+[1])
            continue
        else:
            filteredmems.append(mem+[0])
            prefend=refend
            prefstart=prefstart
    
    filteredmems.sort(key=lambda m: m[1]) #sort by qry position
    mems=[]
    pqrystart=0
    pqryend=0
    for mem in filteredmems:
        qrystart=mem[1]
        qryend=mem[1]+mem[2]
        if qryend<pqryend and pqrystart<qrystart and mem[5]==1: #qry and ref contained
            continue
        else:
            mems.append(mem[:5])
            pqryend=qryend
            pqrystart=pqrystart
    
    after=len(mems)
    
    logging.info("Filtered mums from %d to %d."%(before,after))
    
    return mems

def filtercontainedmumsratio(mems,ratio=.1):
    #filter mems before trying to form chains

    mems.sort(key=lambda m: m[0]) #sort by reference position
    filteredmems=[]
    prefstart=0
    prefend=0
    for mem in mems:
        refstart=mem[0]
        refend=mem[0]+mem[2]
        if refend<prefend and prefstart<refstart: #ref contained
            if float(refend-refstart)/float(prefend-prefstart)<ratio:
                continue

        filteredmems.append(mem)
        prefend=refend
        prefstart=refstart
        pmem=mem
    
    filteredmems.sort(key=lambda m: m[1]) #sort by qry position
    mems=[]
    pqrystart=0
    pqryend=0
    for mem in filteredmems:
        qrystart=mem[1]
        qryend=mem[1]+mem[2]
        if qryend<pqryend and pqrystart<qrystart: #qry contained
            if float(qryend-qrystart)/float(pqryend-pqrystart)<ratio:
                continue
        mems.append(mem)
        pqryend=qryend
        pqrystart=qrystart
        pmem=mem
    

    return mems

def filtermumsrange(mems,maxgapsize=1500,minchainlength=1500):
    #filter mems before trying to form chains
    
    logging.debug("Number of mems before filtering: %d"%len(mems))

    mems.sort(key=lambda m: m[0]) #sort by reference position
    filteredmems=[]

    for i in xrange(len(mems)-1):
        if i==0:
            continue
        if mems[i][2]>=minchainlength:
            filteredmems.append(mems[i])
            continue
        
        prefstart=mems[i-1][0]
        prefend=mems[i-1][0]+mems[i-1][2]
        refstart=mems[i][0]
        refend=mems[i][0]+mems[i][2]
        nrefstart=mems[i+1][0]
        nrefend=mems[i+1][0]+mems[i+1][2]
        
        if abs(refstart-prefend)>maxgapsize:
            continue
        elif abs(refend-nrefstart)>maxgapsize:
            continue
        else:
            filteredmems.append(mems[i])

    filteredmems.sort(key=lambda m: m[1]) #sort by qry position
    mems=[]

    for i in xrange(len(filteredmems)-1):
        if i==0:
            continue
        if filteredmems[i][2]>=minchainlength:
            mems.append(filteredmems[i])
            continue

        pqrystart=filteredmems[i-1][1]
        pqryend=filteredmems[i-1][1]+filteredmems[i-1][2]
        qrystart=filteredmems[i][1]
        qryend=filteredmems[i][1]+filteredmems[i][2]
        nqrystart=filteredmems[i+1][1]
        nqryend=filteredmems[i+1][1]+filteredmems[i+1][2]
        
        if abs(qrystart-pqryend)>maxgapsize:
            continue
        elif abs(qryend-nqrystart)>maxgapsize:
            continue
        else:
            mems.append(filteredmems[i])
    
    logging.debug("Number of mems after filtering: %d"%len(mems))
    
    return mems


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

def bestctgpath(ctgs):
    
    start=(0,0,[],0,0,0,0,0,0,0)
    
    link=dict() #key is combination of ctgname and orientation, which should be unique for each chromosome
    score=dict({():0})
    
    processed=[]
    active=[start]
    maxscore=0
    
    for ctg in ctgs:
        ctgname,revcomp,path,cscore,refbegin,refend,ctgbegin,ctgend,ctglength,ci=ctg
        
        remove=[]
        for pctg in processed:
            pctgname,prevcomp,ppath,pscore,prefbegin,prefend,pctgbegin,pctgend,pctglength,pci=pctg

            if prefend<refend: #may overlap, may not be contained
                active.append(pctg)
                remove.append(pctg)
        
        for r in remove:
            processed.remove(r)
        
        best=None
        w=None

        for actg in active:
            #calculate score of connecting to active point
            actgname,arevcomp,apath,ascore,arefbegin,arefend,actgbegin,actgend,actglength,aci=actg
            
            if w==None:
                w=score[tuple(actg[2])]+cscore
                best=actg
            else:
                if arefend>refbegin:
                    penalty=abs(arefend-refbegin) #penalize by the amount of overlap
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


def bestmempath(mems,ctglength,n=15000,revcomp=False):
    
    if len(mems)>n: #take only n largest mems
        mems.sort(key=lambda mem: mem[2]) #sort by size
        mems=mems[:n]
    
    mems=[tuple(m) for m in mems if m[4]==revcomp]
    
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
    
    return path,maxscore


def mempaths(mems,ctglength,n=15000,revcomp=False,maxgapsize=1500,minchainlength=1500,minchainsum=65):
    
    nmums=len(mems)
    if nmums>n: #take only n largest mems
        logging.info("Too many mums (%d), taking the %d largest."%(nmums,n))
        mems.sort(key=lambda mem: mem[2],reverse=True) #sort by size
        mems=mems[:n] #take the n largest
    
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


def bestmempathwithinversions(mems,n=15000):
    nmums=len(mems)
    if nmums>n: #take only n largest mems
        logging.info("Too many mums (%d), taking the %d largest."%(nmums,n))
        mems.sort(key=lambda mem: mem[2],reverse=True) #sort by size
        mems=mems[:n] #take the n largest

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

def mempathsbothdirections(mems,ctglength,n=15000,maxgapsize=1500,minchainlength=1500,minchainsum=1000):

    nmums=len(mems)
    if nmums>n: #take only n largest mems
        logging.info("Too many mums (%d), taking the %d largest."%(nmums,n))
        mems.sort(key=lambda mem: mem[2],reverse=True) #sort by size
        mems=mems[:n] #take the n largest
    
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
            
            #make active a kd tree, so we can search amems in log(n)! Add amems dynamically to the tree...

            #or remove amems for which amem[0]+amem[2]<mem[0]-maxgapsize
            
            best=init
            w=mem[2]
            
            if mem[4]==1:
                subactive=[amem for amem in active if amem[0]+amem[2]>mem[0]-maxgapsize and amem[1]<mem[1]+mem[2]+maxgapsize and amem[4]==mem[4]]
            else:
                subactive=[amem for amem in active if amem[0]+amem[2]>mem[0]-maxgapsize and amem[1]+amem[2]>mem[1]-maxgapsize and amem[4]==mem[4]]
            
            for amem in subactive:
                
                #calculate score of connecting to active point
                if mem[4]==1:
                    if amem[1] >= mem[1]:
                        p1=(mem[0], mem[1]+mem[2])
                        p2=(amem[0]+amem[2], amem[1])
                        penalty=sumofpairs(p1,p2)
                        #tmpw=score[amem]+mem[2]-penalty
                        tmpw=score[tuple(amem)]+mem[2]-penalty
                        if tmpw>w:
                            w=tmpw
                            best=amem
                else:
                    if amem[1]+amem[2] <= mem[1]+mem[2]:
                        p1=(amem[0]+amem[2], amem[1]+amem[2])
                        p2=(mem[0], mem[1])
                        penalty=sumofpairs(p1,p2)
                        #tmpw=score[amem]+mem[2]-penalty
                        tmpw=score[tuple(amem)]+mem[2]-penalty
                        if tmpw>w:
                            w=tmpw
                            best=amem
            
            #link[mem]=best
            #score[mem]=w
            link[tuple(mem)]=tuple(best)
            score[tuple(mem)]=w
            
            if w>maxscore:
                maxscore=w
                end=mem
            
            processed.append(mem)
        
        path=[]
        o=end[4]
        while end!=start:
            assert(o==end[4])
            path.append(tuple(end))
            #end=link[end]
            end=link[tuple(end)]
        
        if len(path)!=0:
            refstart=path[-1][0]
            refend=path[0][0]+path[0][2]
            
            if o==1:
                ctgstart=path[-1][1]+path[-1][2]
                ctgend=path[0][1]
            else:
                ctgstart=path[-1][1]
                ctgend=path[0][1]+path[0][2]
            
            assert(refend>refstart)
            
            if refend-refstart>minchainlength:
                chainsum=sum([m[2] for m in path])
                if chainsum>=minchainsum:
                    paths.append((path,maxscore,o,ctgstart,ctgend,refstart,refend))
                    logging.debug("Added path consisting of %d mums, from ref:%d:%d to ctg:%d:%d with orientation %d with score %d and length %d."%(len(path),refstart,refend,ctgstart,ctgend,o,maxscore,refend-refstart))
                else:
                    logging.debug("Skipping path consisting of %d mums, from ref:%d:%d to ctg:%d:%d with orientation %d with score %d and length %d and sum %d (minchainsum=%d)."%(len(path),refstart,refend,ctgstart,ctgend,o,maxscore,refend-refstart,chainsum,minchainsum))
            else:
                logging.debug("Skipping path consisting of %d mums, from ref:%d:%d to ctg:%d:%d with orientation %d with score %d and length %d (minchainlength=%d)."%(len(path),refstart,refend,ctgstart,ctgend,o,maxscore,refend-refstart,minchainlength))
            
            if o==1:
                assert(ctgstart>ctgend)
                ctgstart,ctgend=ctgend,ctgstart #flip
            
            umems=[]
            for mem in mems:

                #filter mems that were contained in the extracted chain
                if (mem[0]>=refstart and mem[0]+mem[2]<=refend):
                    continue

                if (mem[1]>=ctgstart and mem[1]+mem[2]<=ctgend):
                    continue
                
                #update mems that were overlapping the ends of the chain, such that we obtain a non-overlapping segmentation of the query sequence
                
                #update left overlapping reference
                #on reference
                if mem[0]<refstart and mem[0]+mem[2]>refstart:
                    mem=(mem[0],mem[1],refstart-mem[0],mem[3],mem[4])
                    #check if now it is contained on ctg domain of the chain
                    if (mem[1]>=ctgstart and mem[1]+mem[2]<=ctgend):
                        continue
                    
                #on contig
                if mem[1]<ctgstart and mem[1]+mem[2]>ctgstart:
                    mem=(mem[0],mem[1],ctgstart-mem[1],mem[3],mem[4])                
                    #check if now it is contained on reference domain of the chain
                    if (mem[0]>=refstart and mem[0]+mem[2]<=refend):
                        continue
                    
                #update right overlapping
                #on reference
                if mem[0]<refend and mem[0]+mem[2]>refend:
                    t=refend-mem[0]
                    mem=(mem[0]+t,mem[1]+t,mem[2]-t,mem[3],mem[4])
                    #check if now it is contained on ctg domain of the chain
                    if (mem[1]>=ctgstart and mem[1]+mem[2]<=ctgend):
                        continue

                #on contig
                if mem[1]<ctgend and mem[1]+mem[2]>ctgend:
                    t=ctgend-mem[1]
                    mem=(mem[0]+t,mem[1]+t,mem[2]-t,mem[3],mem[4])
                    #check if now it is contained on reference domain of the chain
                    if (mem[0]>=refstart and mem[0]+mem[2]<=refend):
                        continue
                
                assert(mem[2]>0)
                
                umems.append(mem)

            mems=umems
            
        else:
            break
    
    logging.debug("Detected number of chains: %d."%len(paths))

    return paths


