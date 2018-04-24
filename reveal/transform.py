import reveallib
import reveallib64
from utils import *
from multiprocessing.pool import Pool
import signal
try:
    from matplotlib import pyplot as plt
    from matplotlib import patches as patches
except:
    pass

def plotchains(ctg2mums,ctg2ref,contig2length,ref2length):

    for ctgname in ctg2ref:
        refs=set([refname for refname,revcomp,score,refbegin,refend,ctgbegin,ctgend,ctglength,ci in ctg2ref[ctgname]])

        for refname in refs:
            plt.ioff()
            plt.clf()
            plt.title(ctgname)
            plt.ylabel(ctgname)
            plt.xlabel(refname)
            mums=sorted(ctg2mums[ctgname][refname],key=lambda m: m[2],reverse=True)[:10000]
            for s1,s2,l,revcomp in mums:
                if revcomp==1:
                    plt.plot([s1,s1+l],[s2+l,s2],'g-')
                else:
                    plt.plot([s1,s1+l],[s2,s2+l],'r-')

            ax = plt.axes()
            last=0

            for ref,revcomp,score,refbegin,refend,ctgbegin,ctgend,ctglength,ci in sorted(ctg2ref[ctgname],key=lambda c: c[6] if c[1]==0 else c[5]):
                logging.info("Plot chain: %s"%str((revcomp,refbegin,refend,ctgbegin,ctgend)))
                if refname != ref:
                    continue

                if revcomp:
                    ctgbegin,ctgend=ctgend,ctgbegin

                if last!=ctgbegin:
                    ax.add_patch(
                        patches.Rectangle(
                            (0, last), #bottom left
                            ref2length[refname], #width
                            ctgbegin-last, #height
                            alpha=.1,
                            color="grey"
                        )
                    )

                ax.add_patch(
                        patches.Rectangle(
                            (refbegin, ctgbegin), #bottom left
                            refend-refbegin, #width
                            ctgend-ctgbegin, #height
                            alpha=.1
                        )
                    )

                last=ctgend

            plt.xlim(0,ref2length[refname])
            plt.ylim(0,contig2length[ctgname])
            plt.show()
            # plt.close()

def transform(args):
    logging.debug("Extracting mums.")

    if args.output==None:
        pref=[]
        for f in [os.path.basename(args.reference),os.path.basename(args.contigs)]:
            bn=os.path.basename(f)
            if '.' in bn:
                pref.append(bn[:bn.find('.')])
            else:
                pref.append(bn)
        args.output="_".join(pref)

    if args.nproc>1:
        original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
        pool = Pool(processes=2 if args.nproc>=2 else 1)
        signal.signal(signal.SIGINT, original_sigint_handler)
        try:
            async_result1 = pool.apply_async(getmums, (args.reference,args.contigs), {'sa64':args.sa64,'minlength':args.minlength,'cutN':args.cutn})
            async_result2 = pool.apply_async(getmums, (args.reference,args.contigs), {'revcomp':True,'sa64':args.sa64,'minlength':args.minlength,'cutN':args.cutn})
        except:
            pool.terminate()
            sys.exit(1)

        pool.close()
        pool.join()

        logging.debug("Retrieving results...")
        mums = async_result1.get()
        logging.debug("Done. MUMS in normal orientation: %d."%len(mums))

        logging.debug("Retrieving RC results...")
        rcmums = async_result2.get()
        logging.debug("Done. MUMS in reverse complemented orientation: %d."%len(rcmums))
    else:
        mums=getmums(args.reference,args.contigs,sa64=args.sa64,minlength=args.minlength,cutN=args.cutn)
        rcmums=getmums(args.reference,args.contigs,revcomp=True,sa64=args.sa64,minlength=args.minlength,cutN=args.cutn)

    reffile=os.path.basename(args.reference)
    ctgfile=os.path.basename(args.contigs)
    
    ref2length=dict()
    for name,seq in fasta_reader(args.reference):
        ref2length[name]=len(seq)
    
    contig2length=dict()
    contig2seq=dict()
    
    totl=0
    for name,seq in fasta_reader(args.contigs,cutN=args.cutn):
        l=len(seq)
        contig2length[name]=l
        totl+=l
        contig2seq[name]=seq

    #combine matches
    mums=mums+rcmums

    if len(mums)==0:
        logging.error("No mums! Exit")
        sys.exit()

    #if args.minlength==None:
    #    args.minlength=1

    #sort by length
    logging.debug("Sorting %d MUMs by size..."%len(mums))
    mums=sorted(mums,key=lambda m: m[4],reverse=True)
    logging.debug("Done.")
    
    if args.minlength==0: #auto determine minlength, prevent use of too many mums
        cov=0
        for i,mem in enumerate(mums):
            cov+=mem[4]
            if cov/float(totl)>1:
                break

        if i<len(mums)-1:
            mums=mums[:i+1]
            logging.info("Over representation of MUMs, auto determined min-mum-length to %d for cov. of %f"%(mums[-1][4],cov/float(totl)))

    ld=[mem[4] for mem in mums]
    bpcovered=sum(ld)

    bpncovered=totl-bpcovered
    if bpncovered<0:
        logging.info("Over representation of MUMs, probably better to use larger -m.")
        bpncovered=1
    
    avgcov=bpcovered/float(totl)

    if args.plot:
        from matplotlib import pyplot as plt
        from matplotlib import patches
        #plot mum lengths
        # plt.hist(ld,bins=range(0,100,1))
        # plt.show()

    if args.minchainsum==None: #auto set minchainsum with 0.5x of the genome wide coverage
        args.minchainsum=int((.5*avgcov)*args.mineventsize)
        logging.info("Auto determined minchainsum to %d"%args.minchainsum)

    logging.info("Assembly consists of %d contigs."%len(contig2seq))
    
    logging.debug("Associating mums to contigs.")
    #relate mums to contigs
    ctg2mums=mapmumstocontig(mums,filtermums=args.filtermums,mineventsize=args.mineventsize)    

    logging.debug("Number of contigs that contain MUMs larger than %d: %d."%(args.minlength,len(ctg2mums)))
    
    logging.info("Using %s to layout the assembly."%args.order)
    if args.order=='chains':
        ref2ctg,ctg2ref=chainstorefence(ctg2mums,contig2length,ref2length,maxmums=args.maxmums,mineventsize=args.mineventsize,minchainsum=args.minchainsum,nproc=args.nproc,extend=args.extend)
    else:
        ref2ctg,ctg2ref=contigstorefence(ctg2mums,contig2length,maxmums=args.maxmums,mineventsize=args.mineventsize,minchainsum=args.minchainsum,nproc=args.nproc)
    
    #write finished assembly based on contigs or chains that map on each reference chromosome
    if not args.split and args.outputtype=='fasta':
        finished=open(args.output+".fasta",'w')
        unplaced=open(args.output+".unplaced.fasta",'w')
    
    totsequnplaced=0
    totseqplaced=0
    totseq=sum(contig2length.values())

    G=nx.MultiDiGraph()
    G.graph['paths']=[]
    G.graph['path2id']=dict()
    G.graph['id2path']=dict()
    G.graph['startnodes']=[]
    G.graph['endnodes']=[]
    
    gapi=0
    pathi=0
    
    defref2ctg=dict()
    unused=[]

    original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
    pool = Pool(processes=args.nproc)
    signal.signal(signal.SIGINT, original_sigint_handler)
    try:
        #multi-process contig-path computation
        for ref in ref2ctg:
            if ref=='unchained' or ref=='unplaced':
                defref2ctg[ref]=ref2ctg[ref]
                continue
            defref2ctg[ref]=pool.apply_async(bestctgpath, (ref2ctg[ref],))
    except KeyboardInterrupt:
        pool.terminate()
    else:
        pool.close()
    pool.join()


    #retrieve multi-process results
    for ref in ref2ctg:
        if ref=='unchained' or ref=='unplaced':
            continue
        
        b=set(ref2ctg[ref])
        defref2ctg[ref]=defref2ctg[ref].get()
        a=set(defref2ctg[ref])
        
        logging.debug("Selected %d out of %d %s to layout assembly with respect to %s."%(len(a),len(b),args.order,ref))
        
        if len(b)-len(a)>0:
            logging.debug("The following %d %s were placed on reference sequence %s but were not used in the layout:"%(len(b)-len(a),args.order,ref))
            if args.order=='contigs':
                for ctgname,revcomp,score,refbegin,refend,ctgbegin,ctgend,ctglength,ci in b - a:
                    logging.debug("Unused: %s (length=%d)"%(ctgname,contig2length[ctgname]))
                    ref2ctg['unplaced'].append(ctgname)
            else:
                for ctgname,revcomp,score,refbegin,refend,ctgbegin,ctgend,ctglength,ci in b - a:
                    if ctgbegin<ctgend:
                        logging.debug("Unused: (%s,%s,%s,%d) (length=%d)"%(ctgname,ctgbegin,ctgend,ci,ctgend-ctgbegin))
                        # ref2ctg['unchained'][ctgname][ctgbegin:ctgend]=0
                    else:
                        logging.debug("Unused: (%s,%s,%s,%d) (length=%d)"%(ctgname,ctgbegin,ctgend,ci,ctgbegin-ctgend))
                        # ref2ctg['unchained'][ctgname][ctgend:ctgbegin]=0
                    unused.append((ctgname,ci))

    # if args.plot:
        # logging.info("Plot chains before extend.")
        # plotchains(ctg2mums,ctg2ref,contig2length,ref2length)

    #remove unused chains from the ctg2ref mapping
    if args.order=="chains":
        defctg2ref=ctg2ref.copy()
        unused.sort(reverse=True)
        for name,i in unused:
            del defctg2ref[name][i]
            uchains=[]
            for chain in defctg2ref[name]:
                ctgname,revcomp,score,refbegin,refend,ctgbegin,ctgend,ctglength,ci=chain
                assert(ci!=i)
                if ci>i:
                    chain=ctgname,revcomp,score,refbegin,refend,ctgbegin,ctgend,ctglength,ci-1
                uchains.append(chain)
            defctg2ref[name]=uchains
        
        keys=sorted(defref2ctg)
        for ref in keys: #update the index of the chains, so that we can detect consecutive chains again
            if ref=='unchained' or ref=='unplaced':
                assert(False)
                continue
            for name,i in unused:
                ctgs=[]
                for ctg in defref2ctg[ref]:
                    ctgname,revcomp,score,refbegin,refend,ctgbegin,ctgend,ctglength,ci=ctg
                    assert(not(ctgname==name and ci==i))
                    if ctgname==name and ci>i:
                        ctg=ctgname,revcomp,score,refbegin,refend,ctgbegin,ctgend,ctglength,ci-1
                    ctgs.append(ctg)
                defref2ctg[ref]=ctgs

        if args.extend:
            defref2ctg,defctg2ref=extendandjoinchains(defref2ctg,defctg2ref,ref2length,contig2length)

        #add parts of contigs that are not part of a chain
        addunchained(defref2ctg,defctg2ref,contig2length)

    else:
        defctg2ref=ctg2ref

    if args.plot:
        logging.debug("Plot chains after extend.")
        plotchains(ctg2mums,defctg2ref,contig2length,ref2length)

    #build graph/fasta for the structural layout of the genome
    for ref in sorted(defref2ctg):
        
        pn=None

        if args.split and args.outputtype=='fasta':
            finished=open(args.output+"_"+ref.replace(" ","_").replace("|","").replace("/","").replace(";","").replace(":","")+".fasta",'w')
            unplaced=open(args.output+"_"+ref.replace(" ","_").replace("|","").replace("/","").replace(";","").replace(":","")+".unplaced.fasta",'w')
                
        if ref=='unchained' or ref=='unplaced':
            continue
        
        logging.debug("Determining %s order for: %s"%(args.order,ref))
        
        ctgs=defref2ctg[ref]
        ctgs.sort(key=lambda c: c[3]) #sort by ref start position

        if args.plot:
            plt.clf()
            #plt.figure(0,figsize=(5,5))
            ax = plt.axes()
            plt.title(args.reference+" vs. "+args.contigs)
        
        coffset=0
        roffset=0
        
        yticks=[]
        yticklabels=[]
        base=os.path.splitext(os.path.basename(args.contigs))[0]
        ctgchromname=base+"_"+ref #name for the finished pseudomolecule
        ctgchromnameorg="*"+base+"_"+ref

        refid=len(G.graph['paths'])
        G.graph['path2id'][ctgchromname]=refid
        G.graph['id2path'][refid]=ctgchromname
        G.graph['paths'].append(ctgchromname)

        startnode=uuid.uuid4().hex
        G.add_node(startnode,offsets={refid:0},endpoint=True)
        G.graph['startnodes'].append(startnode)

        endnode=uuid.uuid4().hex
        G.add_node(endnode,offsets={refid:0},endpoint=True)
        G.graph['endnodes'].append(endnode)

        for ctg in ctgs:
            p="*"+base+"_"+ctg[0] #prefix with asterisk so they're recognisable
            if p not in G.graph['path2id']:
                G.graph['path2id'][p]=len(G.graph['paths'])
                G.graph['id2path'][len(G.graph['paths'])]=p
                G.graph['paths'].append(p)

        if args.outputtype=='fasta':
            finished.write(">%s (finished using %s)\n"%(ctgchromname,ref))
        
        i=0
        o=0

        refpath=[] #path that describes the 'transformed' genome
        orgpath=[] #path that describes the 'original' genome
        
        lastrefchain=False
        lastctgchain=False

        pctg=(None,ctgs[0][1],0,0,0,0,0,0,0)

        for ctg in ctgs:
            ctgname,revcomp,score,refbegin,refend,ctgbegin,ctgend,ctglength,ci=ctg

            if ctg==ctgs[-1]: #the last chain for this chromosome
                lastrefchain=True

            if (ci==len(defctg2ref[ctgname])-1 and not revcomp) or (ci==0 and revcomp): #the last chain for this contig 
                lastctgchain=True
            else:
                lastctgchain=False

            if revcomp:
                ctgbegin,ctgend=ctgend,ctgbegin
            
            pctgname,prevcomp,pscore,prefbegin,prefend,pctgbegin,pctgend,pctglength,pci=pctg
            
            if prevcomp:
                pctgbegin,pctgend=pctgend,pctgbegin
            
            if args.order=='chains':
                reloffset=ctgbegin
            else:
                reloffset=0
            
            if refend<=prefend:
                logging.error("Contained contig should not be in best contig path! %s with alignment length %d"%(ctgname,ctgend-ctgbegin))
                sys.exit(1)

            gapsize=refbegin-prefend
            
            if gapsize<0 or args.fixedsize:
                if gapsize<0:
                    logging.debug("Chains for contigs %s and %s overlap by %d bases."%(pctgname,ctgname,abs(gapsize)))
                gapsize=args.gapsize
            
            logging.debug("%d (index on ctg: %d->%d) - Order %s (revcomp=%d,prefstart=%d,prefend=%d,refstart=%d,refend=%d,ctgstart=%d,ctgend=%d,gapsize=%d)"%(i,pci,ci,args.order,revcomp,prefbegin,prefend,refbegin,refend,ctgbegin,ctgend,gapsize))
            
            if args.order=='chains':
                event=None
                if ((ci==0 and revcomp==0) or (revcomp==1 and ci==len(defctg2ref[ctgname])-1)) and (pctgname==None or ((pci==len(ctg2ref[pctgname])-1 and prevcomp==0) or (prevcomp==1 and pci==0))): #consecutive contigs, no chains in between
                    event='contig break'
                else: #not first or last chain of contig, so has to be stuctural event
                    logging.debug("Non consecutive chains between %s [%d:%d:%d] and %s [%d:%d:%d]."%(pctgname,pctgbegin,pctgend,prevcomp,ctgname,ctgbegin,ctgend,revcomp))
                    if pctgname!=ctgname:
                        event="translocation between contigs" #between contig
                    else:
                        if revcomp!=prevcomp:
                            event="inversion"
                        else:
                            event="translocation within contig" #within contig
                    logging.info("Event of type: \'%s\' between %d and %d."%(event,prefend,refbegin))
                    logging.debug("Index within contig (%s, %d) layout: %d (of %d)"%(ctgname,revcomp,ci,len(defctg2ref[ctgname])))
                    if pctgname!=None:
                        logging.debug("Index within previous contig (%s, %d) layout: %d (of %d)"%(pctgname,prevcomp,pci,len(defctg2ref[pctgname])))
                
                logging.debug("Inserting gap of size: %d"%gapsize)
                
                gap=True
                
                if gapsize==0:
                    if args.outputtype=='fasta':
                        finished.write("N") #write at least one N so we can still distinguish events within fasta
                else:
                    if args.outputtype=='fasta':
                        finished.write("N"*gapsize)
                
                alength=ctgend-ctgbegin
                assert(alength>0)

                l=gapsize+alength
                
                if revcomp:
                    seq=rc(contig2seq[ctgname][ctgbegin:ctgend])
                else:
                    seq=contig2seq[ctgname][ctgbegin:ctgend]
                
                if args.outputtype=='fasta':
                    finished.write(seq)
                
                assert(ctgend-ctgbegin>=0)
                totseqplaced+=ctgend-ctgbegin
                assert(alength==ctgend-ctgbegin)
                
                if args.outputtype=='graph':
                    if event==None: #consecutive chains
                        G.node[pn]['seq']+=seq
                    else: #non-consecutive chains: different contig or structural variant
                        
                        if gapsize>0: #add a gap node
                            gapseq="N"*gapsize
                        else:
                            gapseq=""
                        
                        n=(ctgname,ctgbegin,ctgend,revcomp)
                        
                        G.add_node(n,seq=gapseq+seq,offsets={refid:o,G.graph['path2id']["*"+base+"_"+n[0]]:n[1]})

                        refpath.append(n)
                        
                        if pn!=None:
                            G.add_edge(pn,n,ofrom="+",oto="+",paths={refid})
                        else: #has to be first node for reference chrom
                            G.add_edge(startnode,n,ofrom="+",oto="+",paths={refid})

                        pn=n

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
                        ax.plot([refbegin,refend],[o+gapsize,o+alength+gapsize],'bx')
                        ax.plot([refbegin,refend],[o+gapsize,o+alength+gapsize],'g-')
                    else:
                        ax.plot([refbegin,refend],[o+gapsize,o+alength+gapsize],'bx')
                        ax.plot([refbegin,refend],[o+gapsize,o+alength+gapsize],'r-')
            
            else: # ordering contigs
                gap=True
                
                assert((pctglength-pctgend)>=0)

                alength=contig2length[ctgname]

                if prevcomp:
                    a_prefend=prefend+pctgbegin
                else:
                    a_prefend=prefend+(pctglength-pctgend)
                
                if revcomp:
                    a_refbegin=refbegin-(alength-ctgend)
                else:
                    a_refbegin=refbegin-ctgbegin

                gapsize=a_refbegin-a_prefend
                
                if gapsize==0: #perfect boundary, stil use one N to be able to distinguish the event
                    gapsize=1
                
                if gapsize<0 or args.fixedsize:
                    gapsize=args.gapsize
                
                if pctgname!=None:
                    logging.debug("\'%s\' follows \'%s\' inserting gap of size: %d"%(ctgname[:20],pctgname[:20],gapsize))
                    if args.outputtype=='fasta':
                        finished.write("N"*gapsize)
                
                if revcomp:
                    seq=rc(contig2seq[ctgname])
                    if args.outputtype=='fasta':
                        finished.write(seq) #write the entire contig
                else:
                    seq=contig2seq[ctgname]
                    if args.outputtype=='fasta':
                        finished.write(seq)
                
                if args.outputtype=='graph':

                    gapi+=1
                    n=(gapi)
                    G.add_node(n,seq="N"*gapsize,offsets={refid:o})
                    if pn!=None:
                        G.add_edge(pn,n,ofrom="+",oto="+",paths={refid})
                    pn=n

                    n=(ctgname,0,contig2length[ctgname],revcomp)
                    G.add_node(n,seq=seq,offsets={refid:o+gapsize,G.graph['path2id']["*"+base+"_"+n[0]]:n[1]})
                    
                    if pn!=None:
                        G.add_edge(pn,n,ofrom="+",oto="+",paths={refid})
                    pn=n
                
                #l=gapsize+len(contig2seq[ctgname])
                assert(len(seq)==contig2length[ctgname])
                assert(gapsize>0)

                l=gapsize+contig2length[ctgname]

                assert(contig2length[ctgname]>=0)
                totseqplaced+=contig2length[ctgname]
            
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
                        ax.plot([refbegin,refend],[o+gapsize+(alength-ctgend),o+(alength-ctgbegin)+gapsize],'bx')
                        ax.plot([refbegin,refend],[o+gapsize+(alength-ctgend),o+(alength-ctgbegin)+gapsize],'g-')
                    else:
                        ax.plot([refbegin,refend],[o+gapsize+ctgbegin,o+ctgend+gapsize],'bx')
                        ax.plot([refbegin,refend],[o+gapsize+ctgbegin,o+ctgend+gapsize],'r-')
            
            i+=1
            o=o+l
            yticks.append(o)
            yticklabels.append("%s:%d"%(ctgname[0:15],ctgend))
            # pctg=ctg

            if revcomp:
                pctg=ctgname,revcomp,score,refbegin,refend,ctgend,ctgbegin,ctglength,ci
            else:
                pctg=ctgname,revcomp,score,refbegin,refend,ctgbegin,ctgend,ctglength,ci

        if args.outputtype=='fasta':
            finished.write("\n")
        
        pathi+=2
        
        if args.split and args.outputtype=='fasta':
            finished.close()
            unplaced.close()
        
        logging.debug("Done.")
        
        if args.plot:
            ax.set_yticks(yticks)
            ax.set_yticklabels(yticklabels)
            plt.xlim(0,ref2length[ref])
            if args.interactive:
                plt.show()
            else:
                plt.savefig(args.output+"_"+ref.split()[0]+".png")

        G.add_edge(pn,endnode,ofrom="+",oto="+",paths={refid})

    if args.outputtype=='graph':
        ctgswithevents=set()

        if args.order=="chains":#reconnect the chains based on their layout in the draft assembly
            sortednodes=sorted([n for n in G.nodes() if type(n)!=str])
            pn=sortednodes[0]

            startnode=uuid.uuid4().hex
            G.graph['startnodes'].append(startnode)
            G.add_node(startnode,offsets={G.graph['path2id']["*"+base+"_"+pn[0]]:0},endpoint=True)
            G.add_edge(startnode,pn,ofrom="+",oto="+" if pn[3]==0 else '-',paths={G.graph['path2id']["*"+base+"_"+pn[0]]})

            if len(sortednodes)>1:
                for n in sortednodes[1:]:
                    if n[0]!=pn[0]:
                        startnode=uuid.uuid4().hex
                        G.graph['startnodes'].append(startnode)
                        G.add_node(startnode,offsets={G.graph['path2id']["*"+base+"_"+n[0]]:0},endpoint=True)
                        G.add_edge(startnode,n,ofrom="+",oto="+" if n[3]==0 else '-',paths={G.graph['path2id']["*"+base+"_"+n[0]]})

                        endnode=uuid.uuid4().hex
                        G.graph['endnodes'].append(endnode)
                        G.add_node(endnode,offsets={G.graph['path2id']["*"+base+"_"+pn[0]]:0},endpoint=True) #TODO: correct offset?
                        G.add_edge(pn,endnode,ofrom="+" if pn[3]==0 else '-',oto="+",paths={G.graph['path2id']["*"+base+"_"+pn[0]]})
                    else:
                        ctgswithevents.add("*"+base+"_"+pn[0])
                        G.add_edge(pn,n,ofrom="+" if pn[3]==0 else '-',oto="+" if n[3]==0 else '-',paths={G.graph['path2id']["*"+base+"_"+pn[0]]})
                    pn=n

            endnode=uuid.uuid4().hex
            G.graph['endnodes'].append(endnode)
            G.add_node(endnode,offsets={G.graph['path2id']["*"+base+"_"+pn[0]]:0},endpoint=True) #TODO: correct offset?
            G.add_edge(pn,endnode,ofrom="+" if pn[3]==0 else '-',oto="+",paths={G.graph['path2id']["*"+base+"_"+pn[0]]})

        if not args.allcontigs:
            G.graph['paths']=[sample for sample in G.graph['paths'] if sample in ctgswithevents or not sample.startswith("*")]

    if 'unplaced' in defref2ctg:
        for ctgname in defref2ctg['unplaced']:
            unplaced.write(">%s\n"%(ctgname))
            unplaced.write("%s\n"%contig2seq[ctgname])
    
    if 'unchained' in defref2ctg:
        if len(defref2ctg['unchained'])>0:
            logging.info("The following parts of contigs could not be placed anywhere on the reference sequence.")
            for name in defref2ctg['unchained']:
                # for start,end,i in defref2ctg['unchained'][name]:
                for start,end in defref2ctg['unchained'][name]:
                    logging.info("%s%s (start=%d,end=%d,length=%d,total-contig-length=%d)"%('*' if end-start!=contig2length[name] else '', name,start,end,end-start,contig2length[name]))
                    if args.outputtype=='fasta':
                        unplaced.write(">%s[%d:%d]\n"%(name,start,end))
                        unplaced.write("%s\n"%contig2seq[name][start:end])
                    totsequnplaced+=end-start
    
    if not args.split and args.outputtype=='fasta':
        finished.close()
        unplaced.close()
    
    if args.outputtype=='graph':
        # write_gfa(G,None,outputfile=os.path.splitext(os.path.basename(args.contigs))[0],paths=True)
        write_gfa(G,None,outputfile=args.output,paths=True)
    
    if totseqplaced==0:
        logging.info("No sequence could be placed!")
    else:
        logging.info("%.2f%% (%d out of %d) of the assembly was placed with respect to the reference."% ( (totseqplaced/float(totseq))*100, totseqplaced, totseq ))

def addunchained(defref2ctg,defctg2ref,contig2length):
    #assign unchained parts
    defref2ctg['unchained']=dict()
    for ctg in contig2length: #assign (parts of) contigs that are not part of a chain
        # defref2ctg['unchained'][ctg]=IntervalTree()
        defref2ctg['unchained'][ctg]=[]
        offset=0
        if ctg in defctg2ref:
            defctg2ref[ctg].sort(key=lambda c: c[8])
            for ref,revcomp,score,refstart,refend,ctgstart,ctgend,l,ci in defctg2ref[ctg]:
                logging.debug("Checking %s:%d:%d %d."%(ctg,ctgstart,ctgend,ci))
                if revcomp:
                    ctgstart,ctgend=ctgend,ctgstart
                if ctgstart>offset:
                    logging.debug("Marking %s:%d:%d as unchained."%(ctg,offset,ctgstart))
                    # defref2ctg['unchained'][ctg][offset:ctgstart]=0
                    defref2ctg['unchained'][ctg].append((offset,ctgstart))
                offset=ctgend
        assert(offset<=contig2length[ctg])
        if offset<contig2length[ctg]:
            # defref2ctg['unchained'][ctg][offset:contig2length[ctg]]=0
            defref2ctg['unchained'][ctg].append((offset,contig2length[ctg]))

def extendandjoinchains(ref2ctg,ctg2ref,ref2length,contig2length):
    # extref2ctg={'unchained':dict()}
    extref2ctg={}
    extctg2ref=dict()

    for ref in ref2ctg:
        if ref=="unchained":
            assert(False)
            continue

        ref2ctg[ref]=sorted(ref2ctg[ref],key=lambda c: c[4])
        extref2ctg[ref]=[]

        pchain=None
        join=[]
        for ri,chain in enumerate(ref2ctg[ref]):
            # update=False
            ctgname,revcomp,score,refstart,refend,ctgstart,ctgend,l,ci=chain
            logging.debug("Evaluate chain: %s:%d:%d - %s:%d:%d %d %d"%(ctgname[:10],ctgstart,ctgend,ref[:10],refstart,refend,revcomp,ci))

            #extend first and last chains
            if (ci==0 or ci==len(ctg2ref[ctgname])-1) and (ri==0 or ri==len(ref2ctg[ref])-1):
                logging.debug("Trying to extend: %s"%(ctgname))

                if ci==0 and ri==0 and not revcomp:
                    if ctgstart<refstart:
                        refstart-=ctgstart
                    else:
                        refstart=0
                    ctgstart=0

                if ci==len(ctg2ref[ctgname])-1 and ri==0 and revcomp:
                    if refstart>(contig2length[ctgname]-ctgstart):
                        refstart-=contig2length[ctgname]-ctgstart
                    else:
                        refstart=0
                    ctgstart=contig2length[ctgname]

                if ci==len(ctg2ref[ctgname])-1 and ri==len(ref2ctg[ref])-1 and not revcomp:
                    if ref2length[ref]>refend+(contig2length[ctgname]-ctgend):
                        refend+=contig2length[ctgname]-ctgend
                    else:
                        refend=ref2length[ref]
                    ctgend=contig2length[ctgname]

                if ci==0 and ri==len(ref2ctg[ref])-1 and revcomp:
                    if ref2length[ref]>(refend+ctgstart):
                        refend+=ctgstart
                    else:
                        refend=ref2length[ref]
                    ctgend=0

                chain=ctgname,revcomp,score,refstart,refend,ctgstart,ctgend,l,ci

            if len(extref2ctg[ref])>0:
                pctgname,prevcomp,pscore,prefstart,prefend,pctgstart,pctgend,pl,pci=extref2ctg[ref][-1]
                if pctgname==ctgname:
                    if revcomp==prevcomp:
                        if (not revcomp and ci==pci+1) or (revcomp and ci==pci-1): #consecutive chains, update boundaries
                            pctgname,prevcomp,pscore,prefstart,prefend,pctgstart,pctgend,pl,pci=extref2ctg[ref][-1]
                            logging.debug("Joining chains (%d): %d:%d - %d:%d --> %d:%d for contig: %s"%(revcomp,pctgstart,pctgend,ctgstart,ctgend,pctgstart,ctgend,ctgname))
                            prefend=refend
                            pctgend=ctgend
                            pscore+=score
                            extref2ctg[ref][-1]=(pctgname,prevcomp,pscore,prefstart,prefend,pctgstart,pctgend,pl,ci)
                            extctg2ref[ctgname][-1]=(ref,prevcomp,pscore,prefstart,prefend,pctgstart,pctgend,pl,ci)
                            continue
            
            extref2ctg[ref].append(chain)

            _,revcomp,score,refstart,refend,ctgstart,ctgend,l,ci=chain
            ctgchain=(ref,revcomp,score,refstart,refend,ctgstart,ctgend,l,ci)
            if ctgname not in extctg2ref:
                extctg2ref[ctgname]=[]
            extctg2ref[ctgname].append(ctgchain)

    return extref2ctg, extctg2ref

def decompose_contig(ctg,mums,contiglength,mineventsize=1500,minchainsum=1000,maxmums=15000):

    logging.debug("Determining best chain(s) for: %s"%ctg)
    paths=[]

    results=[]
    for ref in mums:
        rmums=mums[ref]

        candidatepaths=mempathsbothdirections(rmums,contiglength,n=maxmums,mineventsize=mineventsize,minchainsum=minchainsum)
        for path,score,rc,ctgstart,ctgend,refstart,refend in candidatepaths:
            if len(path)>0:
                paths.append((score,ctgstart,ctgend,refstart,refend,ref,rc,path))

    for r in results:
        candidatepaths=r.get()
        for path,score,rc,ctgstart,ctgend,refstart,refend in candidatepaths:
           if len(path)>0:
               paths.append((score,ctgstart,ctgend,refstart,refend,ref,rc,path))
    
    if len(paths)==0:
        return paths
    
    nrefchroms=len(set([p[5] for p in paths]))
    
    logging.debug("Found a total of %d chains for %s that map to %d different reference chromosomes."%(len(paths),ctg,nrefchroms))
     
    paths=sorted(paths,key=lambda c: c[0],reverse=True) #sort chains by alignment score
    
    selectedpaths=[]
    #take n-best paths that dont overlap
    it=IntervalTree()
    for path in paths:
        score,ctgstart,ctgend,refstart,refend,ref,revcomp,p=path
        
        logging.debug("Path: ctg:%d:%d - ref:%s:%d:%d (%d) with score %d"%(ctgstart,ctgend,ref,refstart,refend,revcomp,score))
        
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
                #update boundaries such that path does not overlap other paths
                bctgstart=ctgstart
                bctgend=ctgend
                assert(ctgstart<ctgend)
                
                for start,end,v in s:
                    assert(start<end)
                    if ctgstart<=start and ctgend>=end: #chain contains a smaller chain with better score, reduce to a point
                        ctgend=ctgstart
                        break
                    if ctgstart<=start: #left overlap, update ctgend
                        ctgend=start
                    if ctgend>=end:
                        ctgstart=end
                
                assert(bctgstart!=ctgstart or bctgend!=ctgend) #one of the two has to have been updated
                
                assert(ctgend>=ctgstart)
                if ctgend!=ctgstart:
                    #update path
                    if revcomp:
                        path=(score,ctgend,ctgstart,refstart,refend,ref,revcomp,p) #note that the reference domain is not updated
                    else:
                        path=(score,ctgstart,ctgend,refstart,refend,ref,revcomp,p)

                    if abs(ctgend-ctgstart)>=mineventsize:
                        it[ctgstart:ctgend]=path
                        selectedpaths.append(path)
    
    paths=sorted(selectedpaths,key=lambda c: c[1] if c[6] else c[2]) #sort by endposition on contig

    return paths

def chainstorefence(ctg2mums,contig2length,ref2length,mineventsize=1500,minchainsum=1000,maxmums=15000,nproc=1,extend=True):    
    # ref2ctg={'unchained':dict()}
    ref2ctg={}
    ctg2ref=dict()
    results=dict()

    original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
    pool=Pool(processes=nproc)
    signal.signal(signal.SIGINT, original_sigint_handler)
    try:
        for ctg in ctg2mums:
            results[ctg]=pool.apply_async(decompose_contig,(ctg,ctg2mums[ctg],contig2length[ctg],),{'mineventsize':mineventsize,'minchainsum':minchainsum,'maxmums':maxmums})
    except KeyboardInterrupt:
        pool.terminate()
    else:
        pool.close()

    for ctg in ctg2mums:
        # ref2ctg['unchained'][ctg]=IntervalTree()
        logging.debug("Determining best chains for: %s"%ctg)
        paths=results[ctg].get()
        
        if len(paths)==0:
            logging.info("No valid chains found for contig: %s"%ctg)
            # ref2ctg['unchained'][ctg][0:contig2length[ctg]]=0
            continue
        
        logging.info("Found %d chains for contig: %s"%(len(paths),ctg))
        offset=0
        for i,path in enumerate(paths):
            score,ctgstart,ctgend,refstart,refend,ref,revcomp,chain=path
            if revcomp:
                assert(ctgend<ctgstart)
            else:
                assert(ctgstart<ctgend)
            
            assert(offset<=ctgstart) #should not be any overlap on the contig anymore
            
            if ref in ref2ctg:
                ref2ctg[ref].append((ctg,revcomp,score,refstart,refend,ctgstart,ctgend,contig2length[ctg],i))
            else:
                ref2ctg[ref]=[(ctg,revcomp,score,refstart,refend,ctgstart,ctgend,contig2length[ctg],i)]
            
            if ctg in ctg2ref:
                ctg2ref[ctg].append((ref,revcomp,score,refstart,refend,ctgstart,ctgend,contig2length[ctg],i))
            else:
                ctg2ref[ctg]=[(ref,revcomp,score,refstart,refend,ctgstart,ctgend,contig2length[ctg],i)]
            
            if revcomp==1:
                ctgstart,ctgend=ctgend,ctgstart
            
            if offset!=ctgstart:
                logging.debug("%d:%d:%d --> unchained"%(offset,ctgstart,revcomp))
                # ref2ctg['unchained'][ctg][offset:ctgstart]=i
            logging.debug("%d:%d:%d --> %s:%d:%d"%(ctgstart,ctgend,revcomp,ref,refstart,refend))

            offset=ctgend

        if offset!=contig2length[ctg]:
            assert(offset<contig2length[ctg])
            # ref2ctg['unchained'][ctg][offset:contig2length[ctg]]=i

    return ref2ctg,ctg2ref

def map_contig(ctg,mums,contiglength,mineventsize=1500,minchainsum=1000,maxmums=15000):
    logging.debug("Determining best chain for: %s"%ctg)
    paths=[]
    for ref in mums:
        logging.debug("Checking %s"%ref)
        mpaths=mempathsbothdirections(mums[ref],contiglength,n=maxmums,all=False,mineventsize=mineventsize,minchainsum=minchainsum)
        if len(mpaths)>0:
            path,score,o,ctgstart,ctgend,refstart,refend=mpaths[0]
            paths.append((score,ctgstart,ctgend,refstart,refend,ref,o,path))
    return paths

def contigstorefence(ctg2mums,contig2length,mineventsize=1500,minchainsum=1000,maxmums=15000,nproc=1):
    
    ref2ctg={'unplaced':[]}
    #ref2ctg={}
    ctg2ref=dict()
    results=dict()

    original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
    pool=Pool(processes=nproc)
    signal.signal(signal.SIGINT, original_sigint_handler)
    try:
        for ctg in ctg2mums:
            results[ctg]=pool.apply_async(map_contig,(ctg,ctg2mums[ctg],contig2length[ctg],),{'mineventsize':mineventsize,'minchainsum':minchainsum,'maxmums':maxmums})
    except KeyboardInterrupt:
        pool.terminate()
    else:
        pool.close()
    pool.join()

    for ctg in ctg2mums:
        paths=results[ctg].get()
        if len(paths)==0:
            ref2ctg['unplaced'].append(ctg)
            continue
        
        paths.sort(key=lambda p:p[0],reverse=True) #sort chains by score in descending order, best first
        
        score,ctgstart,ctgend,refstart,refend,ref,revcomp,chain=paths[0] #just take the best path
        
        if ref in ref2ctg:
            ref2ctg[ref].append((ctg,revcomp,score,refstart,refend,ctgstart,ctgend,contig2length[ctg],0))
        else:
            ref2ctg[ref]=[(ctg,revcomp,score,refstart,refend,ctgstart,ctgend,contig2length[ctg],0)]

        if ctg in ctg2ref:
            ctg2ref[ctg].append((ref,revcomp,score,refstart,refend,ctgstart,ctgend,contig2length[ctg],0))
        else:
            ctg2ref[ctg]=[(ref,revcomp,score,refstart,refend,ctgstart,ctgend,contig2length[ctg],0)]
    
    return ref2ctg,ctg2ref

def mapmumstocontig(mums,filtermums=False,mineventsize=1500,minchainsum=1000):
    ctg2mums=dict()
    logging.debug("Mapping %d mums to contigs."%len(mums))
    for mum in mums:
        refchrom, refstart, ctg, ctgstart, l, n, o = mum
        refstart=int(refstart)
        ctgstart=int(ctgstart)
        l=int(l)
        # n=int(n)
        o=int(o)
        if ctg in ctg2mums:
            if refchrom in ctg2mums[ctg]:
                ctg2mums[ctg][refchrom].append([refstart,ctgstart,l,o])
            else:
                ctg2mums[ctg][refchrom]=[[refstart,ctgstart,l,o]]

        else:
            ctg2mums[ctg]=dict({refchrom : [[refstart,ctgstart,l,o]]})
    
    if filtermums:
        for ctg in ctg2mums:
            for ref in ctg2mums[ctg]:
                logging.info("Filtering mums between %s and %s"%(ref,ctg))
                ctg2mums[ctg][ref]=filtermumsrange(ctg2mums[ctg][ref],mineventsize=mineventsize,minchainsum=minchainsum)#,minchainlength=minchainlength)
    
    return ctg2mums

def getmums(reference, query, revcomp=False, sa64=False, minlength=20, cutN=1000):
    if sa64:
        idx=reveallib64.index()
    else:
        idx=reveallib.index()
    
    # G=nx.DiGraph() #dummy for now
    t=IntervalTree()
    
    reffile=os.path.basename(reference)
    ctgfile=os.path.basename(query)
    
    idx.addsample(reffile)
    # if reference.endswith(".gfa"): #dummy for now
    #     logging.error("Not yet supported.")
    #     read_gfa(reference,idx,t,G)
    # else:
    for name,seq in fasta_reader(reference):
        intv=idx.addsequence(seq)
        intv=Interval(intv[0],intv[1],name)
        t.add(intv)
    
    idx.addsample(ctgfile)
    # if query.endswith(".gfa"): #dummy for now
    #     logging.error("Not yet supported.")
    #     read_gfa(query,idx,t,G,revcomp=revcomp)
    # else:
    for name,seq in fasta_reader(query,cutN=cutN):
        if revcomp:
            rcseq=rc(seq)
            intv=idx.addsequence(rcseq)
        else:
            intv=idx.addsequence(seq)
        
        intv=Interval(intv[0],intv[1],name)
        t.add(intv)
    
    idx.construct()
    
    mums=[]
    
    minlength=minlength if minlength!=0 else 1

    logging.debug("Extracting all MUMs of size %d or larger."%minlength)
    
    for mum in idx.getmums(minlength):
        refstart=mum[2][0][1]
        ctgstart=mum[2][1][1]
        rnode=t[refstart].pop() #start position on match to node in graph
        cnode=t[ctgstart].pop()
        if revcomp:
            l=cnode[1]-cnode[0]
            mums.append((rnode[2], refstart-rnode[0], cnode[2], l-((ctgstart-cnode[0])+mum[0]), mum[0], mum[1], 1))
        else:
            mums.append((rnode[2], refstart-rnode[0], cnode[2], ctgstart-cnode[0], mum[0], mum[1], 0))
    
    if revcomp:
        logging.debug("Extracted %d, 3'-5' MUMs (size=%d bytes)."%(len(mums),sys.getsizeof(mums)))
    else:
        logging.debug("Extracted %d, 5'-3' MUMs (size=%d bytes)."%(len(mums),sys.getsizeof(mums)))
    
    return mums

def filtercontainedmumsboth(mums):
    #filter mums before trying to form chains
    before=len(mums)
    mums.sort(key=lambda m: m[0]) #sort by reference position
    filteredmums=[]
    prefstart=0
    prefend=0
    for mem in mums:
        refstart=mem[0]
        refend=mem[0]+mem[2]
        if refend<prefend and prefstart<refstart: #ref contained
            filteredmums.append(mem+[1])
            continue
        else:
            filteredmums.append(mem+[0])
            prefend=refend
            prefstart=prefstart
    
    filteredmums.sort(key=lambda m: m[1]) #sort by qry position
    mums=[]
    pqrystart=0
    pqryend=0
    for mem in filteredmums:
        qrystart=mem[1]
        qryend=mem[1]+mem[2]
        if qryend<pqryend and pqrystart<qrystart and mem[5]==1: #qry and ref contained
            continue
        else:
            mums.append(mem[:5])
            pqryend=qryend
            pqrystart=pqrystart
    
    after=len(mums)
    
    logging.info("Filtered mums from %d to %d."%(before,after))
    
    return mums

def filtercontainedmumsratio(mums,ratio=.1):
    #filter mums before trying to form chains

    mums.sort(key=lambda m: m[0]) #sort by reference position
    filteredmums=[]
    prefstart=0
    prefend=0
    for mem in mums:
        refstart=mem[0]
        refend=mem[0]+mem[2]
        if refend<prefend and prefstart<refstart: #ref contained
            if float(refend-refstart)/float(prefend-prefstart)<ratio:
                continue

        filteredmums.append(mem)
        prefend=refend
        prefstart=refstart
        pmem=mem
    
    filteredmums.sort(key=lambda m: m[1]) #sort by qry position
    mums=[]
    pqrystart=0
    pqryend=0
    for mem in filteredmums:
        qrystart=mem[1]
        qryend=mem[1]+mem[2]
        if qryend<pqryend and pqrystart<qrystart: #qry contained
            if float(qryend-qrystart)/float(pqryend-pqrystart)<ratio:
                continue
        mums.append(mem)
        pqryend=qryend
        pqrystart=qrystart
        pmem=mem
    
    return mums

def filtermumsrange(mums,mineventsize=1500,minchainsum=1000):#minchainlength=1500):
    #filter mums before trying to form chains
    
    logging.debug("Number of mums before filtering: %d"%len(mums))

    mums.sort(key=lambda m: m[0]) #sort by reference position

    filteredmums=[]

    for i in xrange(len(mums)):
        if mums[i][2]>=minchainsum:#minchainlength:
            filteredmums.append(mums[i])
            continue
        
        if len(mums)<=1:
            continue
        
        refstart=mums[i][0]
        refend=mums[i][0]+mums[i][2]
        
        if i==0: #first
            nrefstart=mums[i+1][0]
            nrefend=mums[i+1][0]+mums[i+1][2]
            if abs(refend-nrefstart)>mineventsize:
                continue
            else:
                filteredmums.append(mums[i])
        elif i==len(mums)-1: #last
            prefstart=mums[i-1][0]
            prefend=mums[i-1][0]+mums[i-1][2]
            if abs(refstart-prefend)>mineventsize:
                continue
            else:
                filteredmums.append(mums[i])
        else:
            prefstart=mums[i-1][0]
            prefend=mums[i-1][0]+mums[i-1][2]
            nrefstart=mums[i+1][0]
            nrefend=mums[i+1][0]+mums[i+1][2]
            if abs(refstart-prefend)>mineventsize and abs(refend-nrefstart)>mineventsize:
                continue
            else:
                filteredmums.append(mums[i])
    
    filteredmums.sort(key=lambda m: m[1]) #sort by qry position
    mums=[]

    for i in xrange(len(filteredmums)):
        if filteredmums[i][2]>=minchainsum:#minchainlength:
            mums.append(filteredmums[i])
            continue
        
        if len(filteredmums)<=1:
            continue
        
        qrystart=filteredmums[i][1]
        qryend=filteredmums[i][1]+filteredmums[i][2]
        
        if i==0: #first
            nqrystart=filteredmums[i+1][1]
            nqryend=filteredmums[i+1][1]+filteredmums[i+1][2]
            if abs(qryend-nqrystart)>mineventsize:
                continue
            else:
                mums.append(filteredmums[i])
        elif i==len(filteredmums)-1: #last
            pqrystart=filteredmums[i-1][1]
            pqryend=filteredmums[i-1][1]+filteredmums[i-1][2]
            if abs(qrystart-pqryend)>mineventsize:
                continue
            else:
                mums.append(filteredmums[i])
        else:
            pqrystart=filteredmums[i-1][1]
            pqryend=filteredmums[i-1][1]+filteredmums[i-1][2]
            nqrystart=filteredmums[i+1][1]
            nqryend=filteredmums[i+1][1]+filteredmums[i+1][2]
            if abs(qrystart-pqryend)>mineventsize and abs(qryend-nqrystart)>mineventsize:
                continue
            else:
                mums.append(filteredmums[i])
    
    logging.debug("Number of mums after filtering: %d"%len(mums))
    
    return mums

def filtercontainedmums(mums):
    #filter mums before trying to form chains

    mums.sort(key=lambda m: m[0]) #sort by reference position
    filteredmums=[]
    prefstart=0
    prefend=0
    for mem in mums:
        refstart=mem[0]
        refend=mem[0]+mem[2]
        if refend<prefend and prefstart<refstart: #ref contained
            continue
        else:
            filteredmums.append(mem)
            prefend=refend
            prefstart=prefstart
    
    filteredmums.sort(key=lambda m: m[1]) #sort by qry position
    mums=[]
    pqrystart=0
    pqryend=0
    for mem in filteredmums:
        qrystart=mem[1]
        qryend=mem[1]+mem[2]
        if qryend<pqryend and pqrystart<qrystart: #qry contained
            continue
        else:
            mums.append(mem)
            pqryend=qryend
            pqrystart=pqrystart

    return mums

def bestctgpath(ctgs):
    ctgs.sort(key=lambda c: c[3])
    start=(0,0,0,0,0,0,0,0,0)
    
    link=dict() 
    score=dict({start:0})
    
    processed=[]
    active=[start]
    maxscore=0
     
    for ctg in ctgs:
        ctgname,revcomp,cscore,refbegin,refend,ctgbegin,ctgend,ctglength,ci=ctg
        
        remove=[]
        for pctg in processed:
            pctgname,prevcomp,pscore,prefbegin,prefend,pctgbegin,pctgend,pctglength,pci=pctg

            if prefend<refend: #may overlap, may not be contained
                active.append(pctg)
                remove.append(pctg)
        
        for r in remove:
            processed.remove(r)
        
        best=start
        w=0

        for actg in active:
            #calculate score of connecting to active point
            actgname,arevcomp,ascore,arefbegin,arefend,actgbegin,actgend,actglength,aci=actg
            
            if arefend>refend:
                continue

            if arefend>refbegin:
                penalty=arefend-refbegin #penalize by the amount of overlap
                assert(penalty>0)
            else:
                penalty=0
            
            tmpw=score[actg]+cscore-penalty
            if tmpw>w:
                w=tmpw
                best=actg
        
        assert(best!=None)
        
        link[ctg]=best
        score[ctg]=w
        
        if w>maxscore:
            maxscore=w
            end=ctg
        
        processed.append(ctg)
    
    #backtrack
    minscore=0
    path=[]
    while end[0]!=start[0]:
        path.append(end)
        end=link[end]
    
    return path[::-1]

def mempathsbothdirections(mums,ctglength,n=15000,mineventsize=1500,minchainsum=1000,wscore=1,wpen=1,all=True):
    nmums=len(mums)
    if nmums>n and n!=0: #take only n largest mums
        logging.info("Too many mums (%d), taking the %d largest."%(nmums,n))
        mums.sort(key=lambda mem: mem[2],reverse=True) #sort by size
        mums=mums[:n] #take the n largest
    
    if len(mums)==0:
        return []

    c=sum([m[2] for m in mums])
    logging.debug("Number of anchors: %d",len(mums))
    logging.debug("Sum of anchors: %d", c)
    logging.debug("Length of contig: %d", ctglength)
    logging.debug("Cov ratio: %s"% (c/float(ctglength)) )
    logging.debug("Min chain sum: %d"% minchainsum)
    logging.debug("Max gap size: %d"% mineventsize)
     
    paths=[]
    
    #extract the best path, remove mums that are part of or start/end within the range of the best path, until no more mums remain
    while len(mums)>0:
        
        mums.sort(key=lambda mem: mem[0]+mem[2]) #sort by reference position
        
        init=(None, None, 0, 0, 0, 0)
        link=dict()
        score=dict({init:0})
        active=[]
        processed=[]
        start=init
        end=None
        
        endpoints=[]
        rcendpoints=[]
        pointtomem=dict()
        for mem in mums:
            # if mem[4]==0:
            if mem[3]==0:
                p=(mem[0]+mem[2],mem[1]+mem[2])
                endpoints.append(p)
                pointtomem[p]=mem
            else:
                p=(mem[0]+mem[2],mem[1])
                rcendpoints.append(p)
                pointtomem[p]=mem
        
        #build kdtrees
        memtree=kdtree(endpoints,2)
        rcmemtree=kdtree(rcendpoints,2)
        maxscore=0
        
        for mem in mums:
            best=init
            w=wscore*mem[2]
            
            # if mem[4]==1:
            if mem[3]==1:
                frompoint=(mem[0]-mineventsize, mem[1])
                topoint=(mem[0]+mem[2]-1, mem[1]+(mem[2]-1)+mineventsize)
                assert(frompoint[0]<topoint[0])
                assert(frompoint[1]<topoint[1])
                assert((topoint[1]-frompoint[1])==(topoint[0]-frompoint[0]))
                subactive=[pointtomem[p] for p in range_search(rcmemtree,frompoint,topoint)]
            else:
                frompoint=(mem[0]-mineventsize, mem[1]-mineventsize)
                topoint=(mem[0]+mem[2]-1, mem[1]+mem[2]-1)
                assert(frompoint[0]<topoint[0])
                assert(frompoint[1]<topoint[1])
                assert((topoint[1]-frompoint[1])==(topoint[0]-frompoint[0]))
                subactive=[pointtomem[p] for p in range_search(memtree,frompoint,topoint)]
            
            subactive.sort(key=lambda s: score[tuple(s)], reverse=True)

            for amem in subactive:
                if score[tuple(amem)]+(wscore*mem[2])<w:
                    break

                #calculate score of connecting to active point
                if mem[3]==1:
                    p1=(mem[0], mem[1]+mem[2])
                    p2=(amem[0]+amem[2], amem[1])
                    penalty=gapcost(p1,p2)
                    tmpw=score[tuple(amem)]+(wscore*mem[2])-(wpen*penalty)
                    if tmpw>w:
                        w=tmpw
                        best=amem
                else:
                    p1=(amem[0]+amem[2], amem[1]+amem[2])
                    p2=(mem[0], mem[1])
                    penalty=gapcost(p1,p2)
                    tmpw=score[tuple(amem)]+(wscore*mem[2])-(wpen*penalty)
                    if tmpw>w:
                        w=tmpw
                        best=amem
            
            link[tuple(mem)]=tuple(best)
            score[tuple(mem)]=w
            
            if w>maxscore:
                maxscore=w
                end=mem
        
        path=[]
        o=end[3]

        while end!=start:
            # assert(o==end[4])
            assert(o==end[3])
            path.append(tuple(end))
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
            
            if not all: #just return first best path
                logging.debug("Best path contains %d anchors (refstart=%d, refend=%d)."%(len(path),refstart,refend))
                return [(path,maxscore,o,ctgstart,ctgend,refstart,refend)]

            chainsum=sum([m[2] for m in path])
            if chainsum>=minchainsum:
                paths.append((path,maxscore,o,ctgstart,ctgend,refstart,refend))
                logging.debug("Added path consisting of %d mums, from ref:%d:%d to ctg:%d:%d with orientation %d with score %d and length %d."%(len(path),refstart,refend,ctgstart,ctgend,o,maxscore,refend-refstart))
            else:
                logging.debug("Skipping path consisting of %d mums, from ref:%d:%d to ctg:%d:%d with orientation %d with score %d and length %d and sum %d (minchainsum=%d)."%(len(path),refstart,refend,ctgstart,ctgend,o,maxscore,refend-refstart,chainsum,minchainsum))
                return paths
            
            if o==1:
                assert(ctgstart>ctgend)
                ctgstart,ctgend=ctgend,ctgstart #flip
            
            #TODO: can be sped up by looking up mums from kdtree!
            umums=[]
            for mem in mums:
                
                #filter mums that were contained in the extracted chain
                if (mem[0]>=refstart and mem[0]+mem[2]<=refend):
                    continue

                if (mem[1]>=ctgstart and mem[1]+mem[2]<=ctgend):
                    continue
                
                #update mums that were overlapping the ends of the chain, such that we obtain a non-overlapping segmentation of the query sequence
                
                #update left overlapping reference
                #on reference
                if mem[0]<refstart and mem[0]+mem[2]>refstart:
                    # mem=(mem[0],mem[1],refstart-mem[0],mem[3],mem[4])
                    # print "update left overlapping reference",mem,(mem[0],mem[1],refstart-mem[0],mem[3])
                    mem=(mem[0],mem[1],refstart-mem[0],mem[3])
                    #check if now it is contained on ctg domain of the chain
                    if (mem[1]>=ctgstart and mem[1]+mem[2]<=ctgend):
                        continue
                    
                #on contig
                if mem[1]<ctgstart and mem[1]+mem[2]>ctgstart:
                    # mem=(mem[0],mem[1],ctgstart-mem[1],mem[3],mem[4])
                    # print "update left overlapping ctg",mem,(mem[0],mem[1],ctgstart-mem[1],mem[3])
                    mem=(mem[0],mem[1],ctgstart-mem[1],mem[3])
                    #check if now it is contained on reference domain of the chain
                    if (mem[0]>=refstart and mem[0]+mem[2]<=refend):
                        continue
                    
                #update right overlapping
                #on reference
                if mem[0]<refend and mem[0]+mem[2]>refend:
                    t=refend-mem[0]
                    # mem=(mem[0]+t,mem[1]+t,mem[2]-t,mem[3],mem[4])
                    # print "update right overlapping ref",mem,(mem[0]+t,mem[1]+t,mem[2]-t,mem[3])
                    mem=(mem[0]+t,mem[1]+t,mem[2]-t,mem[3])
                    #check if now it is contained on ctg domain of the chain
                    if (mem[1]>=ctgstart and mem[1]+mem[2]<=ctgend):
                        continue

                #on contig
                if mem[1]<ctgend and mem[1]+mem[2]>ctgend:
                    t=ctgend-mem[1]
                    # mem=(mem[0]+t,mem[1]+t,mem[2]-t,mem[3],mem[4])
                    # print "update right overlapping ctg",mem,(mem[0]+t,mem[1]+t,mem[2]-t,mem[3])
                    mem=(mem[0]+t,mem[1]+t,mem[2]-t,mem[3])
                    #check if now it is contained on reference domain of the chain
                    if (mem[0]>=refstart and mem[0]+mem[2]<=refend):
                        continue
                
                assert(mem[2]>0)
                
                umums.append(mem)

            mums=umums
        else:
            break
    
    logging.debug("Detected number of chains: %d."%len(paths))
    
    return paths


