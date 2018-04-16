
from utils import *
import reveallib
import reveallib64
import os


def matches(args): 
    
    if args.sa64:
        idx=reveallib64.index(sa=args.sa1, lcp=args.lcp1, cache=args.cache) #enable preconstruction of first SA and LCP array
    else:
        idx=reveallib.index(sa=args.sa1, lcp=args.lcp1, cache=args.cache) #enable preconstruction of first SA and LCP array
    
    G=nx.DiGraph()
    G.graph['paths']=[]
    t=IntervalTree()
    
    reffile=os.path.basename(args.reference)
    ctgfile=os.path.basename(args.contigs)
    
    ref2length=dict()
    idx.addsample(reffile)
    if args.reference.endswith(".gfa"):
        read_gfa(args.reference,idx,t,G)
    else:
        G.graph['paths'].append(reffile)
        for name,seq in fasta_reader(args.reference):
            ref2length[name]=len(seq)
            intv=idx.addsequence(seq)
            intv=Interval(intv[0],intv[1],name)
            t.add(intv)
            G.add_node(intv,offsets={reffile:0})
    
    contig2length=dict()
    idx.addsample(ctgfile)
    if args.contigs.endswith(".gfa"):
        read_gfa(args.contigs,idx,t,G)
    else:
        G.graph['paths'].append(ctgfile)
        for name,seq in fasta_reader(args.contigs):
            contig2length[name]=len(seq)
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
    
    #for each contig, print the length
    for name in contig2length:
        print "#%s\t%d"%(name,contig2length[name])
    
    idx.construct()
    
    if args.uniq:
        print "##refname\trefstart\tctgname\tctgstart\tlength\tn\torient"
        for mem in idx.getmums(args.minlength):
            refstart=mem[2][0]
            ctgstart=mem[2][1]
            rnode=t[refstart].pop() #start position on match to node in graph
            cnode=t[ctgstart].pop()
            print "%s\t%s\t%s\t%s\t%s\t%s\t%s" % (rnode[2], refstart-rnode[0], cnode[2], ctgstart-cnode[0], mem[0], mem[1], 0)
    else:
        print "##refname\trefstart\tctgname\tctgstart\tlength\tn\tunique\torient"
        for mem in idx.getmems(args.minlength):
            refstart=mem[2][0]
            ctgstart=mem[2][1]
            rnode=t[refstart].pop() #start position on match to node in graph
            cnode=t[ctgstart].pop()
            print "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (rnode[2], refstart-rnode[0], cnode[2], ctgstart-cnode[0], mem[0], mem[1], mem[3], 0)
    
    if args.rc:

        logging.debug("Indexing reverse complement...\n")
        
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
            rcG.graph['paths']=set([reffile])
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
            rcG.graph['paths']=set([ctgfile])
            for name,seq in fasta_reader(args.contigs):
                intv=idx.addsequence(rc(seq))
                intv=Interval(intv[0],intv[1],name)
                t.add(intv)
                rcG.add_node(intv,offsets={ctgfile:0},aligned=0)
        
        idx.construct()
        
        if args.uniq:
            for mem in idx.getmums(args.minlength):
                refstart=mem[2][0]
                ctgstart=mem[2][1]
                rnode=t[refstart].pop() #start position on match to node in graph
                cnode=t[ctgstart].pop()
                l=cnode[1]-cnode[0]
                print "%s\t%s\t%s\t%s\t%s\t%s\t%s" % (rnode[2], refstart-rnode[0], cnode[2], l-((ctgstart-cnode[0])+mem[0]), mem[0], mem[1], 1)
        else:
            for mem in idx.getmems(args.minlength):
                refstart=mem[2][0]
                ctgstart=mem[2][1]
                rnode=t[refstart].pop() #start position on match to node in graph
                cnode=t[ctgstart].pop()
                l=cnode[1]-cnode[0]
                print "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (rnode[2], refstart-rnode[0], cnode[2], l-((ctgstart-cnode[0])+mem[0]), mem[0], mem[1], mem[3], 1)
