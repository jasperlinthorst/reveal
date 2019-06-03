import reveallib
import reveallib64
from utils import *

def gplot(args):
    G=nx.DiGraph()
    read_gfa(args.graph,None,None,G)
    
    if len(G.graph['paths'])<2:
        logging.error("Can't make a plot for less than two samples.")
        sys.exit(1)

    if args.x==None and args.y==None:
        args.x=G.graph['paths'][0]
        args.y=G.graph['paths'][1]
    else:
        if args.x not in G.graph['paths']:
            logging.error("%s not contained in the specified graph. Samples that are: %s."%(args.x,G.graph['paths']))
            sys.exit(1)
        if args.y not in G.graph['paths']:
            logging.error("%s not contained in the specified graph. Samples that are: %s."%(args.y,G.graph['paths']))
            sys.exit(1)
    plotgraph(G, args.x, args.y, interactive=args.interactive, region=args.region, minlength=args.minlength)

def bedplot(args):
    from matplotlib import pyplot as plt

    if len(args.fastas)==1:
        with open(args.fastas[0]) as bedfile:
            xpoints,rcxpoints=[],[]
            ypoints,rcypoints=[],[]
            pref=None
            xoffset=0
            for line in bedfile:
                if line.startswith("#"):
                    continue

                if pref!=reference:
                    off

                reference,refbegin,refend,contig,score_cost,orientation,alnstart,alnend=line.rstrip().split()
                name,ctgidx,lastsegmentidx,ctgbegin,ctgend=contig.split(":")
                refbegin,refend,alnstart,alnend,ctgbegin,ctgend=[int(v) for v in [refbegin,refend,alnstart,alnend,ctgbegin,ctgend]]


                if orientation=='-':
                    rcxpoints.append(alnstart)
                    rcxpoints.append(alnend)
                    rcxpoints.append(None)
                    rcypoints.append(ctgend)
                    rcypoints.append(ctgbegin)
                    rcypoints.append(None)
                else:
                    xpoints.append(alnstart)
                    xpoints.append(alnend)
                    xpoints.append(None)
                    ypoints.append(ctgbegin)
                    ypoints.append(ctgend)
                    ypoints.append(None)

            print len(xpoints)
            plt.plot(xpoints,ypoints,'r-')
            plt.plot(rcxpoints,rcypoints,'g-')
            plt.plot(1,1)
            plt.show()

def plot(args):

    import matplotlib

    if not args.interactive:
        matplotlib.use('Agg')

    from matplotlib import pyplot as plt
    from matplotlib import patches as patches

    vertgaps=[]
    horzgaps=[]
    vertgapsizes=[]
    horzgapsizes=[]
    ctgoffsets=[]
    refoffsets=[]
    qrylength=0
    reflength=0
    ax = plt.axes()
    
    if len(args.fastas)==2:
        if args.sa64:
            idx=reveallib64.index()
        else:
            idx=reveallib.index()
        
        ctgid=0
        
        sample=args.fastas[0]
        idx.addsample(sample)
        refoffset=0
        for name,seq in fasta_reader(sample):
            pc=None
            gapsize=None
            for i,c in enumerate(seq):
                if c=='N' and pc!='N':
                    horzgaps.append(i)
                    gapsize=1
                elif c=='N' and pc=='N':
                    gapsize+=1
                elif c!='N' and pc=='N':
                    horzgapsizes.append(gapsize)
                pc=c
            refoffset+=i+2
            reflength+=len(seq)+1
            refoffsets.append(refoffset)
            intv=idx.addsequence(seq.upper())
        
        sample=args.fastas[1]
        idx.addsample(sample)
        qryoffset=0
        for name,seq in fasta_reader(sample):
            pc=None
            gapsize=None
            for i,c in enumerate(seq):
                if c=='N' and pc!='N':
                    vertgaps.append(qryoffset+i)
                    gapsize=1
                elif c=='N' and pc=='N':
                    gapsize+=1
                elif c!='N' and pc=='N':
                    vertgapsizes.append(gapsize)
                pc=c
            qryoffset+=i+2
            qrylength+=len(seq)+1
            ctgoffsets.append(qryoffset)
            intv=idx.addsequence(seq.upper())
        
        qrylength=qrylength-1
        idx.construct()
        
        logging.info("Extracting mums...")
        mmems=idx.getmums(args.minlength)
        logging.info("Done.")
        
        sep=idx.nsep[0]

        if args.rc:
            
            #get mums for reverse orientation
            idx.construct(rc=True)
            
            logging.info("Extracting RC mums...")
            mmems+=idx.getmums(args.minlength)
            logging.info("Done.")
     
    elif len(args.fastas)==1 and args.fastas[0].endswith(".bed"):
        bedplot(args)
        return
    else:
        logging.fatal("Can only create mumplot for 2 sequences or self plot for 1 sequence.")
        return
    
    start=0
    end=sep
    qend=idx.n

    del idx
    
    if len(mmems)>args.maxmums:
        logging.info("Too many mums (%d), taking the %d largest."%(len(mmems),args.maxmums))
        mmems.sort(key=lambda mem: mem[0],reverse=True) #sort by size
        mmems=mmems[:args.maxmums] #take the n largest
    
    logging.info("Drawing %d matches."%len(mmems))
    
    xlist,rcxlist = [],[]
    ylist,rcylist = [],[]
    
    for mem in mmems:
        # sps=sorted(mem[2])
        sps=mem[1]
        l=mem[0]
        
        sp1=sps[0]
        sp2=sps[1]-(sep+1)
        ep1=sp1+l
        ep2=sp2+l
        
        if sp1>=start and ep1<=end:
            
            if mem[2]==0:
                xlist.append(sp1)
                xlist.append(ep1)
                ylist.append(sp2)
                ylist.append(ep2)
                xlist.append(None)
                ylist.append(None)    
            else:
                rcxlist.append(ep1)
                rcxlist.append(sp1)
                rcylist.append(sp2)
                rcylist.append(ep2)
                rcxlist.append(None)
                rcylist.append(None)

    plt.plot(xlist,ylist,'r-')
    plt.plot(rcxlist,rcylist,'g-')
    
    if args.endpoints:
        plt.plot(xlist,ylist,'b*')
        plt.plot(rcxlist,rcylist,'y*')

    for p in ctgoffsets:
        plt.axhline(y=p,linewidth=.5,color='black',linestyle='solid')
    
    for p in refoffsets:
        plt.axvline(x=p,linewidth=.5,color='black',linestyle='solid')
    
    if args.showgaps:
        for p,l in zip(horzgaps,horzgapsizes):
            ax.add_patch(
                patches.Rectangle(
                    (p, 0), #bottom left
                    l, #width
                    qrylength, #height
                    alpha=.1
                )
            )
        
        for p,l in zip(vertgaps,vertgapsizes):
            ax.add_patch(
                patches.Rectangle(
                    (0, p), #bottom left
                    reflength, #width
                    l, #height
                    alpha=.1
                )
            )
        
    plt.xlim(start,end)
    plt.ylim(0,qend-end)
    plt.title(" vs. ".join(args.fastas))
    if len(args.fastas)==2:
        plt.xlabel(args.fastas[0])
        plt.ylabel(args.fastas[1])
    else:
        plt.xlabel(args.fastas[0])
        plt.xlabel(args.fastas[0]+"_rc")
    plt.autoscale(enable=False)
    
    if args.xregion!=None:
        xregions=[]

        for region in args.xregion.split(","):

            if region.count("-")==1:
                rstart,rend=region.split("-") #should be rectangle with alfa here
            elif region.count(":")==1:
                rstart,rend=region.split(":") #should be rectangle with alfa here
            else:
                logging.fatal("Invalid region specification, use - : <start>-<end>")
                sys.exit(1)

            xregions.append((int(rstart),int(rend)))
            plt.axvline(x=int(rstart),linewidth=1,color='b',linestyle='dashed')
            plt.axvline(x=int(rend),linewidth=1,color='b',linestyle='dashed')

    if args.yregion!=None:
        yregions=[]

        for region in args.yregion.split(","):
            
            if region.count("-")==1:
                rstart,rend=region.split("-") #should be rectangle with alfa here
            elif region.count(":")==1:
                rstart,rend=region.split(":") #should be rectangle with alfa here
            else:
                logging.fatal("Invalid region specification, use - : <start>-<end>")
                sys.exit(1)

            yregions.append((int(rstart),int(rend)))
            plt.axhline(y=int(rstart),linewidth=1,color='b',linestyle='dashed')
            plt.axhline(y=int(rend),linewidth=1,color='b',linestyle='dashed')

    if args.interactive:
        plt.show()
    else:
        b1=os.path.basename(args.fastas[0])
        b2=os.path.basename(args.fastas[1])
        
        fn1=b1[:b1.rfind('.')] if b1.find('.')!=-1 else b1
        fn2=b2[:b2.rfind('.')] if b2.find('.')!=-1 else b2

        if args.xregion!=None and args.yregion!=None:
            assert(len(xregions)==len(yregions))
            
            if args.flanksize!=None:
                flanksizes=[int(v) for v in args.flanksize.split(",")]
            else:
                flanksizes=[0]*len(xregions)

            for xregion,yregion,flanksize in zip(xregions,yregions,flanksizes):
                plt.xlim(xregion[0]-flanksize,xregion[1]+flanksize)
                plt.ylim(yregion[0]-flanksize,yregion[1]+flanksize)
                plt.savefig(fn1+"_"+str(xregion[0])+"-"+str(xregion[1])+"_"+fn2+"_"+str(yregion[0])+"-"+str(yregion[1])+"."+args.extension)
        else:
            plt.savefig(fn1+"_"+fn2+"."+args.extension)
