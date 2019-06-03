import sys
import os
import uuid

def align(args):

    graphs=[args.reference[0]]

    step=0
    
    if args.transform:
        step+=1
        sys.stdout.write("#(%d) Convert draft assemblies to graphs (address rearrangements)\n"%step)
        for file in args.inputfiles:
            graph=os.path.splitext(file)[0]+'.gfa'
            sys.stdout.write("reveal transform %s %s %s -o %s\n"% ("--64" if args.sa64 else "", args.reference[0], file, graph))
            graphs.append(graph)
    else:
        graphs=args.reference+args.inputfiles

    step+=1
    sys.stdout.write("#(%d) Use REM to construct an anchor based alignment graph (brake down the problem)\n"%step)
    tmpfiles=[]

    #TODO: use tree based/progressive approach or simultaneous anchor based alignment
    if args.order=='sequential':
        level=0

        while len(graphs)>1:
            step+=1
            sys.stdout.write("#(%d) Level (%d) alignments\n"%(step,level))
            n=args.chunksize
            k,m=divmod(len(graphs),n)
            if k==0:
                chunks=[graphs]
                graphs=[]
            else:
                chunks=[graphs[i*n:i*n+n] for i in range(k)]
                if m!=0:
                    graphs=graphs[-m:]
                else:
                    graphs=[]

            for chunk in chunks:
                if len(chunks)==1 and graphs==[]: #final merge
                    sys.stdout.write("reveal rem %s %s -o %s.gfa\n"% ("--64" if args.sa64 else "", " ".join(chunk), args.output))
                    graphs.append(args.output+".gfa")
                else:
                    tmp=uuid.uuid4().hex
                    sys.stdout.write("reveal rem %s %s -o %s.gfa\n"% ("--64" if args.sa64 else "", " ".join(chunk), tmp))
                    graphs.append(tmp+".gfa")
                    tmpfiles.append(tmp+".gfa")
            level+=1

    elif args.order=='tree':
        logging.error("Tree-based construction not yet supported.")
        sys.exit(1)
    else: #attempt simultaneous
        sys.stdout.write("reveal rem %s %s -m%d %s -o %s.gfa\n"% ("--64" if args.sa64 else ""," ".join(graphs),args.m, "-n "+str(args.n) if args.n!=None else "",args.output) )

    if len(tmpfiles)>0:
        step+=1
        sys.stdout.write("#(%d) Cleanup tempfiles\n"%step)
        for tmp in tmpfiles:
            sys.stdout.write("rm %s\n"%tmp)

    if args.unzip:
        step+=1
        sys.stdout.write("#(%d) Unzip all bubbles in the graph\n"%step)
        sys.stdout.write("reveal unzip %s.gfa -u10\n"%(args.output))

    if args.refine:
        step+=1
        sys.stdout.write("#(%d) Refine all bubbles in the graph using MSA\n"%step)
        sys.stdout.write("reveal refine %s.unzipped.gfa --nproc=%d --all --maxsize=10000 --minsize=2 --mindiff=0 --minconf=%d\n"%(args.output,args.nproc,args.minconf))

    if args.variants:
        step+=1
        sys.stdout.write("#(%d) Output variants\n"%step)
        sys.stdout.write("reveal variants %s.gfa --vcf > %s.anchored.vcf\n" %(args.output,args.output))
        if args.unzip:
            sys.stdout.write("reveal variants %s.unzipped.gfa --vcf > %s.unzipped.vcf\n" %(args.output,args.output))
        if args.refine:
            sys.stdout.write("reveal variants %s.unzipped.realigned.gfa --vcf > %s.refined.vcf\n" %(args.output,args.output))
        