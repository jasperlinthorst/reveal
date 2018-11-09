import reveallib
import reveallib64
from utils import *

try:
    from matplotlib import pyplot as plt
    from matplotlib import patches as patches
except:
    pass

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

def plot(args):
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
        
        print "Extracting mums..."
        #mmems=[(mem[0],mem[1],mem[2].values(),0) for mem in idx.getmums(args.minlength)]
        # mmems=[(mem[0],mem[1],[sp for gid,sp in mem[2]],0) for mem in idx.getmums(args.minlength)]
        # mmems=[(mem[0],2,[sp for sp in mem[1]],0) for mem in idx.getmums(args.minlength)]
        mmems=idx.getmums(args.minlength)#[(mem[0],mem[1]) for mem in idx.getmums(args.minlength)]

        sep=idx.nsep[0]

        if args.rc:
            
            #get mums for reverse orientation
            idx.construct(rc=True)

            # if args.sa64:
            #     idx=reveallib64.index()
            # else:
            #     idx=reveallib.index()
            
            # sample=args.fastas[0]
            # idx.addsample(sample)
            # for name,seq in fasta_reader(sample):
            #     idx.addsequence(seq.upper())
            
            # sample=args.fastas[1]
            # idx.addsample(sample)

            # qryintvs=[]
            # for name,seq in fasta_reader(sample):
            #     intv=idx.addsequence(rc(seq.upper()))
            #     qryintvs.append(intv)
            
            # idx.construct()
            
            print "Extracting RC mums..."
            mmems+=idx.getmums(args.minlength)
            
            # vi=iter(qryintvs)
            # v=vi.next()
            
            #tmp=[(m[0],m[1],sorted(m[2].values())) for m in tmp] #make sure start positions are sorted
            # tmp=[(m[0],m[1],sorted([sp for gid,sp in m[2]])) for m in tmp] #make sure start positions are sorted
            # tmp=[(m[0],2,sorted([sp for sp in m[1]])) for m in tmp]
            # tmp=[(m[0],m[1]) for m in tmp]
            # tmp.sort(key=lambda l: l[2][1]) #sort by query pos
            
            # nmmems=[]
            # for mem in tmp:
            #     if mem[2][1]>v[1]:
            #         v=vi.next()
            #     start,end=v
            #     newqstart=end-(mem[2][1]-start)-mem[0]
            #     ntup=(mem[0],mem[1],(mem[2][0],newqstart),1)
            #     nmmems.append(ntup)
            
            # mmems+=nmmems
            
            print "done."
     
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
    
    print "Drawing",len(mmems),"matches."
    
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
                plt.plot([sp1,ep1],[sp2,ep2],'r-')
            else:
                plt.plot([ep1,sp1],[sp2,ep2],'g-')
    
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
        for region in args.xregion.split(","):
            rstart,rend=region.split(":") #should be rectangle with alfa here
            plt.axvline(x=int(rstart),linewidth=3,color='b',linestyle='dashed')
            plt.axvline(x=int(rend),linewidth=3,color='b',linestyle='dashed')

    if args.yregion!=None:
        for region in args.yregion.split(","):
            rstart,rend=region.split(":") #should be rectangle with alfa here
            plt.axhline(y=int(rstart),linewidth=3,color='b',linestyle='dashed')
            plt.axhline(y=int(rend),linewidth=3,color='b',linestyle='dashed')

    # plt.title("")
    # plt.xticks([], [])
    # plt.yticks([], [])
    # plt.xlabel("")
    # plt.ylabel("")

    # plt.xlim(15000000,21000000)
    # plt.ylim(15000000,21000000)

    if args.interactive:
        plt.show()
    else:
        b1=os.path.basename(args.fastas[0])
        b2=os.path.basename(args.fastas[1])
        fn1=b1[0:args.fastas[0].rfind('.')] if b1.find('.')!=-1 else b1
        fn2=b2[0:args.fastas[1].rfind('.')] if b2.find('.')!=-1 else b2
        plt.savefig(fn1+"_"+fn2+"."+args.extension)
