#!/usr/bin/env python

from intervaltree import Interval, IntervalTree
import networkx as nx
from collections import defaultdict, deque
import threading
import reveallib
import argparse
import logging
import os
import schemes
import sys

def fasta_reader(fn):
    seq=""
    with open(fn,'r') as ff:
        for line in ff:
            if line.startswith(">"):
                if seq!="":
                    yield name,seq
                name=line.rstrip().replace(">","")
                seq=""
            else:
                seq+=line.rstrip()
        if seq!="":
            yield name,seq

def rc(seq):
    d = d = {'A':'T','C':'G','G':'C','T':'A','N':'N','a':'t','c':'g',\
        'g':'c','t':'a','n':'n','Y':'R','R':'Y','K':'M','M':'K',\
        'S':'S','W':'W','B':'V','V':'B','D':'H','H':'D','N':'N',\
        'X':'X','-':'-'}
    return "".join([d[b] for b in reversed(seq)])

def breaknode(node,pos,l):
    att=G.node[node]
    in_edges=G.in_edges(node)
    out_edges=G.out_edges(node)
    mn=Interval(pos,pos+l)
    other=set()
    if mn==node: #no breaking needed
        t.remove(node)
        return node,other
    G.add_node(mn,sample=att['sample'],contig=att['contig'],coordsample=att['coordsample'],coordcontig=att['coordcontig'],start=att['start']+(pos-node.begin),aligned=0)#create merge node
    if (node[0]!=pos):
        pn=Interval(node[0],pos)
        G.add_node(pn,sample=att['sample'],contig=att['contig'],coordsample=att['coordsample'],coordcontig=att['coordcontig'],start=att['start'],aligned=0)#create prefix node
        assert(not G.has_edge(pn,mn))
        assert(not G.has_edge(mn,pn))
        G.add_edge(pn,mn)
        t.add(pn)
        other.add(pn)
    else:
        pn=mn
    if (node[1]!=pos+l):
        sn=Interval(pos+l,node[1])
        G.add_node(sn,sample=att['sample'],contig=att['contig'],coordsample=att['coordsample'],coordcontig=att['coordcontig'],start=att['start']+((pos+l)-node.begin),aligned=0)#create suffix node
        assert(not G.has_edge(mn,sn))
        assert(not G.has_edge(sn,mn))
        G.add_edge(mn,sn)
        t.add(sn)
        other.add(sn)
    else:
        sn=mn
    G.remove_node(node)                     #update Graph
    t.remove(node)                          #update intervaltree
    for e in in_edges:
        assert(not G.has_edge(e[0],pn))
        assert(not G.has_edge(pn,e[0]))
        G.add_edge(e[0],pn)
    for e in out_edges:
        assert(not G.has_edge(sn,e[1]))
        assert(not G.has_edge(e[1],sn))
        G.add_edge(sn,e[1])
    return mn,other #return merge node

def mergenodes(mns):
    if reference!=None:
        for i,node in enumerate(mns):
            if reference in G.node[node]['sample']:
                refnode=node
                break
        else:
            refnode=mns[0]
    else:
        refnode=mns[0]
    
    G.node[refnode]['sample']=set.union(*[G.node[node]['sample'] for node in mns])
    G.node[refnode]['contig']=set.union(*[G.node[node]['contig'] for node in mns])
    G.node[refnode]['aligned']=1
    mns.remove(refnode)
    for mn in mns: #leave the first node, remove the rest
        for e in G.in_edges(mn):
            assert(e[0]!=refnode)
            #assert(not G.has_edge(e[0],refnode))
            #if G.has_edge(refnode,e[0]):
            #    print mns, refnode
            #    write_gml(G,T,outputfile="test.gml")
            assert(not G.has_edge(refnode,e[0]))
            G.add_edge(e[0],refnode)
        for e in G.out_edges(mn):
            assert(e[1]!=refnode)
            assert(not G.has_edge(e[1],refnode))
            #assert(not G.has_edge(refnode,e[1]))
            G.add_edge(refnode,e[1])
        G.remove_node(mn)
    return refnode

def bfs(G, source, reverse=False):
    if reverse and isinstance(G, nx.DiGraph):
        neighbors = G.predecessors_iter
    else:
        neighbors = G.neighbors_iter
    visited = set([source])
    queue = deque([(source, neighbors(source))])
    while queue:
        parent, children = queue[0]
        try:
            child = next(children)
            if child not in visited:
                visited.add(child)
                if not(G.node[child].has_key('aligned')):
                    queue.append((child, neighbors(child)))
                    yield child,0
                elif (G.node[child]['aligned']==0):
                    queue.append((child, neighbors(child)))
                    yield child,0
                else:
                    yield child,1
        except StopIteration:
            queue.popleft()

def segmentgraph(node,nodes):
    trailing=set()
    leading=set()
    reverse_trailing=set()
    reverse_leading=set()
    nodes=set(nodes)
    endpoints=set()
    for c,t in bfs(G,node):
        if t==0:
            trailing.add(c)
        else:
            endpoints.add(c)
        
    #reverse search for each endpoint
    if len(endpoints)>1:
        for endpoint in endpoints:
            for c,t in bfs(G,endpoint,reverse=True):
                if t==0:
                    reverse_trailing.add(c)
        trailing=trailing.intersection(reverse_trailing)
    
    endpoints=set()
    for c,t in bfs(G,node,reverse=True):
        if t==0:
            leading.add(c)
        else:
            endpoints.add(c)
    
    #reverse search for each endpoint
    if len(endpoints)>1:
        for endpoint in endpoints:
            for c,t in bfs(G,endpoint):
                if t==0:
                    reverse_leading.add(c)
        leading=leading.intersection(reverse_leading)
    
    leading = set([(i.begin,i.end) for i in leading if isinstance(i,Interval)]).intersection(nodes)
    trailing = set([(i.begin,i.end) for i in trailing if isinstance(i,Interval)]).intersection(nodes)
    rest = nodes - (leading | trailing)
    return list(leading), list(trailing), list(rest)

def graphalign(l,index,n,score,sp):
    nodes=index.nodes

    if l<schemes.minlength or len(nodes)==0 or score<schemes.minscore:
        return
    
    mns=[]
    topop=[]
    for i,pos in enumerate(sp):
        #assert(len(t[pos])==1) #be sure that there are no overlapping intervals in the tree!
        #assert(len(t[pos+l-1])==1) #be sure that there are no overlapping intervals in the tree!
        #assert(t[pos]==t[pos+l-1])
        #if i>0:
        #    assert(T[sp[i]:sp[i]+l]==T[sp[i-1]:sp[i-1]+l])
        old=t[pos].pop()
        mn,other=breaknode(old,pos,l)
        mns.append(mn)
        nodes.remove((old.begin,old.end))
        for node in other:
            nodes.append((node.begin,node.end))
    mn=mergenodes(mns)
    intervals=segmentgraph(mn,nodes)
    return intervals

def prune_nodes(G,T):
    converged=False
    reverse=False
    while not(converged):
        converged=True
        for node,data in G.nodes_iter(data=True):
            for run in [0,1]:
                if data['aligned']==1:
                    if run==0:
                        reverse=False
                        neis=G.successors(node)
                    else:
                        reverse=True
                        neis=G.predecessors(node)
                    seqs={}
                    for nei in neis:
                        if G.node[nei]['aligned']!=1:
                            if 'seq' not in G.node[nei]:
                                seq=T[nei.begin:nei.end]
                            else:
                                seq=G.node[nei]['seq']
                            if seq in seqs:
                                seqs[seq].append(nei)
                            else:
                                seqs[seq]=[nei]
                    for key in seqs.keys():
                        group=seqs[key]
                        if len(group)>1:
                            sink=[]
                            for v in group:
                                sink+=G.neighbors(v)
                            if len(set(sink))>1:
                                continue
                            if len(G.subgraph(group).edges())>0: #if interconnected don't merge!
                                continue
                            mergenodes(group)
                            converged=False

def read_gfa(gfafile, index, tree, graph, minsamples=1, maxsamples=None, targetsample=None):
    f=open(gfafile,'r')
    sep=";"
    mapping={} #temp mapping object
    edges=[] #tmp list for edges
    for line in f:
        if line.startswith('H'):
            h=line.split()
            tag=h[1].split(':')
            if tag[0]=="ORI":
                for sample in tag[2].rstrip(sep).split(sep):
                    if 'samples' in graph.graph:
                        graph.graph['samples'].append(sample)
                    else:
                        graph.graph['samples']=[sample]
        if line.startswith('S'):
            s=line.strip().split('\t')
            nodeid=s[1]
            ann={}
            for v in s[4:]:
                v=v.split(':')
                if v[2].find(sep)!=-1: #multi-valued
                    ann[v[0]]=v[2].rstrip(sep).split(sep)
                else:
                    ann[v[0]]=v[2]
            
            if "ORI" not in ann: #not a reveal graph, so no metadata on nodes, just index all
                if index!=None:
                    intv=index.addsequence(s[2].upper())
                    intv=Interval(intv[0],intv[1])
                    tree.add(intv)
                    graph.add_node(intv,sample={gfafile},contig={nodeid},aligned=0)
                    mapping[nodeid]=intv
                else:
                    graph.add_node(gfafile+'_'+nodeid,sample={gfafile},contig={nodeid},seq=s[2].upper(),aligned=0)
                    mapping[nodeid]=gfafile+'_'+nodeid
            else: #there are annotations on the nodes, so use them
                if not(isinstance(ann['ORI'], list)):
                    ann['ORI']=[ann['ORI']]
                if not(isinstance(ann['CTG'], list)):
                    ann['CTG']=[ann['CTG']]
                
                ann['ORI']=set(ann['ORI'])
                ann['CTG']=set(ann['CTG'])
                
                if len(ann['ORI'])<minsamples: #dont index these nodes, just add them to the graph
                    graph.add_node(gfafile+'_'+nodeid,sample=ann['ORI'],contig=ann['CTG'],coordsample=ann['CRD'],coordcontig=ann['CRDCTG'],start=int(ann['START']),seq=s[2].upper(),aligned=0)
                    mapping[nodeid]=gfafile+'_'+nodeid
                elif maxsamples!=None and len(ann['ORI'])>maxsamples: #dont index these nodes, just add them to the graph
                    graph.add_node(gfafile+'_'+nodeid,sample=ann['ORI'],contig=ann['CTG'],coordsample=ann['CRD'],coordcontig=ann['CRDCTG'],start=int(ann['START']),seq=s[2].upper(),aligned=0)
                    mapping[nodeid]=gfafile+'_'+nodeid
                elif (targetsample!=None and targetsample not in ann['ORI']): #dont index these nodes, just add them to the graph
                    graph.add_node(gfafile+'_'+nodeid,sample=ann['ORI'],contig=ann['CTG'],coordsample=ann['CRD'],coordcontig=ann['CRDCTG'],start=int(ann['START']),seq=s[2].upper(),aligned=0)
                    mapping[nodeid]=gfafile+'_'+nodeid
                else:
                    if index!=None:
                        intv=index.addsequence(s[2].upper())
                        intv=Interval(intv[0],intv[1])
                        tree.add(intv)
                        graph.add_node(intv,sample=ann['ORI'],contig=ann['CTG'],coordsample=ann['CRD'],coordcontig=ann['CRDCTG'],start=int(ann['START']),aligned=0)
                        mapping[nodeid]=intv
                    else:
                        graph.add_node(gfafile+'_'+nodeid,sample=ann['ORI'],contig=ann['CTG'],coordsample=ann['CRD'],coordcontig=ann['CRDCTG'],start=int(ann['START']),seq=s[2].upper(),aligned=0)
                        mapping[nodeid]=gfafile+'_'+nodeid
        
        #L	206	+	155	+	0M
        if line.startswith('L'):
	    edges.append(line)

    for line in edges:
	e=line.strip().split()
        assert(not graph.has_edge(mapping[e[1]],mapping[e[3]]))
        assert(not graph.has_edge(mapping[e[3]],mapping[e[1]]))
        graph.add_edge(mapping[e[1]],mapping[e[3]],ofrom=e[2],oto=e[4],cigar=e[5])



def write_gfa(G,T,outputfile="reference.gfa",vg=False):
    f=open(outputfile,'wb')
    sep=';'
    f.write('H\tVN:Z:1.0\n')
    f.write('H\tORI:Z:')
    for sample in G.graph['samples']:
        for subsample in sample.rstrip().rstrip(sep).split(sep):
            f.write(subsample+sep)
    f.write('\n')
    
    mapping={}

    for i,node in enumerate(nx.topological_sort(G)): #iterate once to get a mapping of ids to intervals
	mapping[node]=i+1 #one-based for vg
    
    for i,node in enumerate(nx.topological_sort(G)):
	i+=1
	data=G.node[node]
        if 'seq' in data:
            f.write('S\t'+str(i)+'\t'+data['seq'].upper())
        else:
            f.write('S\t'+str(i)+'\t'+T[node.begin:node.end].upper())
        
	if not vg:
            f.write("\t*\tORI:Z:%s\tCRD:Z:%s\tCRDCTG:Z:%s\tCTG:Z:%s\tSTART:Z:%s\n" % 
                    (
                    sep.join(data['sample']),
                    data['coordsample'].replace(sep,"").replace("\t",""),
                    data['coordcontig'].replace(sep,"").replace("\t",""),
                    sep.join(data['contig']),
                    data['start']
                    )
                )
        else:
            f.write("\n")
            for sample in data['sample']:
                f.write("P\t"+str(i)+"\t"+sample+"\t+\t"+str(node.end-node.begin)+"M\n")
	
	for to in G[node]:
	    f.write('L\t'+str(mapping[node])+'\t+\t'+str(mapping[to])+"\t+\t0M\n")
	
    f.close()

def write_gml(G,T,outputfile="reference",hwm=1000):
    mapping={}
    for n,d in G.nodes(data=True):
        mapping[n]=str(n)
        d=G.node[n]
        G.node[n]['n']=len(d['sample'])
        G.node[n]['sample']=str(d['sample'])
        G.node[n]['contig']=str(d['contig'])
        G.node[n]['coordsample']=str(d['coordsample'])
        G.node[n]['coordcontig']=str(d['coordcontig'])
        G.node[n]['start']=str(d['start'])
        if d.has_key('aligned'):
            G.node[n]['aligned']=str(d['aligned'])
        if 'seq' not in G.node[n]:
            G.node[n]['seq']=T[n.begin:n.end].upper()
        G.node[n]['l']=len(G.node[n]['seq'])
    G=nx.relabel_nodes(G,mapping,copy=True)
    outputfiles=[]
    
    i=0
    for subset in nx.connected_components(G.to_undirected()):
	sgn=[]
	g=G.subgraph(subset)
	for n in nx.topological_sort(g):
	    sgn.append(n)
	    if G.node[n]['n']==len(G.graph['samples']): #join node
		if len(sgn)>=hwm:
		    sg=G.subgraph(sgn)
		    fn=outputfile+'.'+str(i)+'.gml'
		    nx.write_gml(sg,fn)
		    outputfiles.append(fn)
		    sgn=[]
		    i+=1
	if len(sgn)>0:
	    sg=G.subgraph(sgn)
	    fn=outputfile+'.'+str(i)+'.gml'
	    nx.write_gml(sg,fn)
	    i+=1
	    outputfiles.append(fn)
    
    return outputfiles

def main():
    desc="""
    Type 'reveal <positional_argument> --help' for help on a specific subcommand.\n
    Reveal constructs population reference graphs by aligning multiple whole 
    genomes using recursive exact matching.
    http://www.biorxiv.org/content/early/2015/07/17/022715.
    """
    
    parser = argparse.ArgumentParser(prog="reveal", usage="reveal -h for usage", description=desc)
    parser.add_argument("-l", "--log-level", dest="loglevel", default=20, help="Log level: 10=debug 20=info (default) 30=warn 40=error 50=fatal.")
    
    subparsers = parser.add_subparsers()
    parser_aln = subparsers.add_parser('align',prog="reveal align", description="Construct a population graph from input genomes or other graphs.")
    parser_call = subparsers.add_parser('call',prog="reveal call", description="Extract variants from a graph.")
    parser_plot = subparsers.add_parser('plot', prog="reveal plot", description="Generate mumplot for two fasta files.")
    #parser_convert = subparsers.add_parser('convert')
    parser_extract = subparsers.add_parser('extract', prog="reveal extract", description="Extract the input sequence from a graph.")
    
    parser_aln.add_argument('inputfiles', nargs='*', help='Fasta or gfa files specifying either assembly/alignment graphs (.gfa) or sequences (.fasta). When only one gfa file is supplied, variants are called within the graph file.')
    parser_aln.add_argument("-o", "--output", dest="output", help="Prefix of the variant and alignment graph files to produce, default is \"sequence1_sequence2\"")
    #parser_aln.add_argument("-p", dest="pcutoff", type=float, default=1e-3, help="If, the probability of observing a MUM of the observed length by random change becomes larger than this cutoff the alignment is stopped (default 1e-3).")
    parser_aln.add_argument("-t", "--threads", dest="threads", type=int, default=0, help = "The number of threads to use for the alignment.")
    parser_aln.add_argument("-m", dest="minlength", type=int, default=20, help="Min length of an exact match (default 20).")
    parser_aln.add_argument("-c", dest="minscore", type=int, default=0, help="Min score of an exact match (default 0), exact maches are scored by their length and penalized by the indel they create with respect to previously accepted exact matches.")
    parser_aln.add_argument("-n", dest="minsamples", type=int, default=1, help="Only align nodes that occcur in this many samples (default 1).")
    parser_aln.add_argument("-x", dest="maxsamples", type=int, default=None, help="Only align nodes that have maximally this many samples (default None).")
    parser_aln.add_argument("-r", dest="reference", type=str, default=None, help="Name of the sequence that should be used as a coordinate system or reference.")
    parser_aln.add_argument("-s", dest="targetsample", type=str, default=None, help="Only align nodes in which this sample occurs.")
    parser_aln.add_argument("--gml", dest="gml", action="store_true", default=False, help="Produce a gml graph instead gfa.")
    parser_aln.add_argument("--gml-max", dest="hwm", default=4000, help="Max number of nodes per graph in gml output.")
    parser_aln.add_argument("--vg", dest="vg", action="store_true", default=False, help="Produce a gfa graph without node annotations, to ensure it's parseable by vg.")
    parser_aln.add_argument("--align-contigs", dest="contigs", action="store_true", default=False, help="Use when pairwise aligning a set of contigs to a single genome or graph. Contigs are aligned one by one (slow).")
    parser_aln.set_defaults(func=align)

    parser_call.add_argument('graph', nargs=1, help='A graph in gfa format from which variants should be extracted.')    
    parser_call.add_argument("--minlen", dest="minlen", type=int, default=1, help="Output variants in a gfa graph that are larger or equal to this size (default=1).")
    parser_call.add_argument("--maxlen", dest="maxlen", type=int, default=1000, help="Output variants in a gfa graph that are smaller or equal to this size (to limit memory consumption for nw-alignments for inversions) (default=1000).")
    parser_call.add_argument("--snps", dest="snps", action="store_true", default=False, help="Output snps in a gfa graph.")
    parser_call.add_argument("--indels", dest="indels", action="store_true", default=False, help="Output indels in a gfa graph.")
    parser_call.add_argument("--inv", dest="inv", action="store_true", default=False, help="Output inversions in a gfa graph.")
    parser_call.add_argument("--multi", dest="multi", action="store_true", default=False, help="Output variants with multiple alleles in a gfa graph.")
    parser_call.add_argument("--all", dest="allvar", action="store_true", default=False, help="Output all variants in a gfa graph.")
    parser_call.set_defaults(func=call)
    
    parser_extract.add_argument('graph', nargs=1, help='gfa file specifying the graph from which the genome should be extracted.')
    parser_extract.add_argument('samples', nargs='*', help='Name of the sample to be extracted from the graph.')
    parser_extract.add_argument("--width", dest="width", type=int, default=100 , help='Line width for fasta output.')
    parser_extract.set_defaults(func=extract)
    
    parser_plot.add_argument('fastas', nargs=2, help='Two fasta files for which a dotplot should be generated.')
    parser_plot.add_argument("-m", dest="minlength", type=int, default=100, help="Minimum length of exact matches to vizualize (default=100).")
    #parser_plot.add_argument("", dest="pos", default=None, help="Position on the reference genome to vizualize.")
    #parser_plot.add_argument("", dest="env", default=100, help="Size of the region aroung the targeted position to vizualize.")
    parser_plot.set_defaults(func=plot)

    #parser_plot.set_defaults(func=convert)
    
    args = parser.parse_args()
    logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p', level=args.loglevel)
    args.func(args)


def align(args):

    if len(args.inputfiles)<=1:
        logging.fatal("Specify at least 2 (g)fa files for creating a reference graph.")
        return
    
    if args.contigs:
	if len(args.inputfiles)!=2:
	    logging.fatal("Only pairwise alignment of contigs is possible.")
            return
	else:
	    G,idx=align_contigs(args)
    else:
	G,idx=align_genomes(args)
    
    if args.output==None:
        args.output="_".join([os.path.basename(f)[:os.path.basename(f).rfind('.')] for f in args.inputfiles])
    
    print "Merging nodes..."
    T=idx.T
    prune_nodes(G,T) #TODO: do this after every alignment step to reduce graph size during alignment
    print "Done."
    
    print "Writing graph..."
    if args.gml:
        graph=write_gml(G,T, hwm=args.hwm, outputfile=args.output)
    else:
        write_gfa(G,T,vg=args.vg, outputfile=args.output+'.gfa')
        graph=args.output+'.gfa'
    print "Done."
    
    print "Alignment graph written to:",graph


def align_genomes(args):
    #globale variables to simplify callbacks from c extension
    global t,G,reference
     
    reference=args.reference
    
    t=IntervalTree()
    idx=reveallib.index()
    G=nx.DiGraph()
    G.graph['samples']=[]
    #schemes.pcutoff=args.pcutoff
    schemes.minlength=args.minlength
    schemes.minscore=args.minscore
    
    for i,sample in enumerate(args.inputfiles):
        idx.addsample(os.path.basename(sample))
        if sample.endswith(".gfa"):
            #TODO: now applies to all graphs! probably want to have this graph specific if at all...
            read_gfa(sample,idx,t,G,minsamples=args.minsamples,
                                    maxsamples=args.maxsamples,
                                    targetsample=args.targetsample)
        else: #consider it to be a fasta file
            G.graph['samples'].append(sample)
            for name,seq in fasta_reader(sample):
                intv=idx.addsequence(seq.upper())
                Intv=Interval(intv[0],intv[1])
                t.add(Intv)
                schemes.ts.add(Intv)
                schemes.interval2sampleid[Intv]=i
                G.add_node(Intv,sample={sample},contig={name.replace(";","")},coordsample=sample,coordcontig=name.replace(";",""),start=0,aligned=0)
     
    print "Constructing index..."
    idx.construct()
    print "Done."
    
    print "Aligning genomes..."
    if len(args.inputfiles)>2:
        idx.align(schemes.mumpicker2,graphalign,threads=args.threads)
    else:
        idx.align(None,graphalign,threads=args.threads)
    print "Done."
    return G,idx

def align_contigs(args):
    global t,G,reference
    
    reference=args.reference
    schemes.minlength=args.minlength
    schemes.minscore=args.minscore
    G=nx.DiGraph()
    G.graph['samples']=[]
    
    ref=args.inputfiles[0]
    contigs=args.inputfiles[1]
    
    t=IntervalTree()
    idx=reveallib.index()
    
    idx.addsample(os.path.basename(ref))
    if ref.endswith(".gfa"):
        read_gfa(ref,idx,t,G,minsamples=args.minsamples,
                                maxsamples=args.maxsamples,
                                targetsample=args.targetsample)
    else:
        G.graph['samples'].append(ref)
        for chromname,chrom in fasta_reader(ref):
            intv=idx.addsequence(chrom.upper())
            Intv=Interval(intv[0],intv[1])
            t.add(Intv)
            G.add_node(Intv,sample={ref},contig={chromname.replace(";","")},coordsample=ref,coordcontig=chromname.replace(";",""),start=0,aligned=0)
    
    i=0
    G.graph['samples'].append(contigs)
    for contigname,contig in fasta_reader(contigs):
	idx.addsample(os.path.basename(contigs))
	intv=idx.addsequence(contig.upper())
	Intv=Interval(intv[0],intv[1])
	t.add(Intv)
	G.add_node(Intv,sample={contigs},contig={contigname.replace(";","")},coordsample=contigs,coordcontig=contigname.replace(";",""),start=0,aligned=0)
	
	print "Constructing index...",idx.n
        idx.construct()
	print "Done."
	
	print "Aligning contig",contigname,len(contig),"..."
        idx.align(None,graphalign,threads=args.threads)
        identity=0
        for node,data in G.nodes(data=True):
            if data['aligned']==1 and isinstance(node,Interval):
                identity+=node.end-node.begin
	print "Done, %d bp out of %d bp aligned (%d%%)." % (identity,len(contig),round(100*(identity/float(len(contig))),2))
        
        t=IntervalTree()
        T=idx.T
        del(idx)
	idx=reveallib.index()
	idx.addsample(os.path.basename(ref))
        #add sequence of all nodes in G that are not spanned by the alignment
        mapping={}
        for node,data in G.nodes(data=True):
            if isinstance(node,Interval):
                neis=G.neighbors(node)
                if len(neis)<=1 and data['aligned']==0 and contigs not in data['sample']: #use unconnected and unaligned nodes for next run
                    intv=idx.addsequence(T[node.begin:node.end])
                    Intv=Interval(intv[0],intv[1])
                    t.add(Intv)
                    mapping[node]=Intv
                    continue
                
                if len(neis)==2 and data['aligned']==0 and contigs not in data['sample']: #also index nodes that connect two alignments
                    if len(nx.all_simple_paths(G, source=neis[0], target=neis[1]))==1 and G.node[neis[0]]['aligned']==1 and G.node[neis[1]]['aligned']==1:
                        intv=idx.addsequence(T[node.begin:node.end])
                        Intv=Interval(intv[0],intv[1])
                        t.add(Intv)
                        mapping[node]=Intv
                        continue
                
                #all other cases, simply add the sequence to the graph, but dont index it
                data['seq']=T[node.begin:node.end].upper()
                mapping[node]=i
                i+=1
        
        G=nx.relabel_nodes(G,mapping, copy=True)
    return G,idx

#extract variants from gfa graphs
def call(args):
    if args.graph[0].endswith(".gfa"):
        import caller
        logging.info("Parsing GFA file.")
        c=caller.Caller(args.graph[0])
        logging.info("Writing variants.")
        c.call(minlength=args.minlen,maxlength=args.maxlen,allvar=args.allvar,inversions=args.inv,indels=args.indels,snps=args.snps,multi=args.multi)

def plot(args):
    from matplotlib import pyplot as plt
    
    if len(args.fastas)!=2:
        logging.fatal("Can only create mumplot for 2 sequences.")
        return
    
    idx=reveallib.index()
    for sample in args.fastas:
        idx.addsample(sample)
        for name,seq in fasta_reader(sample):
            intv=idx.addsequence(seq.upper())
    idx.construct()
    mmems=[(mem[0],mem[1],mem[2],0) for mem in idx.getmums() if mem[0]>args.minlength]
    sep=idx.nsep[0]
    
    idx=reveallib.index()
    sample=args.fastas[0]
    idx.addsample(sample)
    for name,seq in fasta_reader(sample):
        intv=idx.addsequence(seq.upper())
    
    sample=args.fastas[1]
    idx.addsample(sample)
    ls=0
    for name,seq in fasta_reader(sample):
        ls+=len(seq)
        intv=idx.addsequence(rc(seq.upper()))
    idx.construct()
    mmems+=[(mem[0],mem[1],mem[2],1) for mem in idx.getmums() if mem[0]>args.minlength]
    
    print len(mmems),"max exact matches."

    pos=None
    
    for mem in mmems:
        sps=sorted(mem[1:3])
        l=mem[0]
        sp1=sps[0]
        sp2=sps[1]-sep
        ep1=sp1+l
        ep2=sp2+l
        
        if pos!=None:
            if abs(sp1-pos)>dist:
                continue
            
            if mem[3]==0 and pos_==None:
                pos_=sp2
            
            if mem[3]==1 and pos_==None:
                pos_=ls-sp2
            
            if mem[3]==0 and abs(sp2-pos_)>dist:
                #print sp2, pos_, l
                continue
            
            if mem[3]==1 and abs((ls-sp2)-pos_)>dist:
                #print ls-sp2, pos_, l
                continue
        
        if mem[3]==0:
            plt.plot([sp1,ep1],[sp2,ep2],'r-')
        else:
            assert(mem[3]==1)
            plt.plot([sp1,ep1],[ls-sp2,ls-ep2],'g-')

    #plt.savefig(str(pos)+'.png')
    plt.title(args.fastas[0]+" vs. "+args.fastas[1])
    plt.xlabel(args.fastas[0])
    plt.ylabel(args.fastas[1])
    plt.autoscale(enable=False)
    plt.show()


#TODO: convert gfa graph to gml (sub)graphs
#def convert(args):
#    return

def extract(args):
    if not args.graph[0].endswith(".gfa"):
        logging.fatal("Invalid gfa file.")
        return
    G=nx.DiGraph()
    read_gfa(args.graph[0], None, None, G)
    for sample in args.samples:
        sg=[]
        for node,data in G.nodes(data=True):
            if sample in data['sample']:
                sg.append(node)
        width=args.width
        for i,contig in enumerate(nx.connected_components(G.subgraph(sg).to_undirected())):
            seq=""
            try:
                sys.stdout.write(">"+sample+" "+str(i)+"\n")
                for node in nx.topological_sort(G.subgraph(contig)):
                    seq+=G.node[node]['seq']
                f=0
                for i in xrange(width,len(seq),width):
                    sys.stdout.write(seq[f:i]+'\n')
                    f=i
                sys.stdout.write(seq[f:]+'\n')
            except IOError:
                try:
                    sys.stdout.close()
                except IOError:
                    pass
                try:
                    sys.stderr.close()
                except IOError:
                    pass
