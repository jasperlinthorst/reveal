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
    if l<schemes.minlength or len(nodes)==0 or score<0:
        #print l, schemes.minlength, len(nodes), n, sp, score
        return
    #assert(n==len(sp))
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
                    graph.graph['samples'].append(sample)
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
                intv=index.addsequence(s[2].upper())
                intv=Interval(intv[0],intv[1])
                tree.add(intv)
                graph.add_node(intv,sample={gfafile},contig={nodeid},aligned=0)
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
                    intv=index.addsequence(s[2].upper())
                    intv=Interval(intv[0],intv[1])
                    tree.add(intv)
                    graph.add_node(intv,sample=ann['ORI'],contig=ann['CTG'],coordsample=ann['CRD'],coordcontig=ann['CRDCTG'],start=int(ann['START']),aligned=0)
                    mapping[nodeid]=intv
                    
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
                    data['coordsample'],
                    data['coordcontig'],
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
	
    #for node1,node2 in G.edges_iter():
    #    f.write('L\t'+str(mapping[node1])+'\t+\t'+str(mapping[node2])+"\t+\t0M\n")
    f.close()

def write_gml(G,T,outputfile="reference",hwm=1000):
    mapping={}
    for n in nx.nodes(G):
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

#return a subgraph or sequence that contains (if all is well) the input sequence/graph from the specified sample
#def extract(sample,graph):
#    for node in graph.


def main():
    usage = "reveal [options] <sequence1.(g)fa> <sequence2.(g)fa> ... <sequenceN.(g)fa>"
    parser = argparse.ArgumentParser(usage)
    parser.add_argument('inputfiles', nargs='*', help='Fasta or gfa files specifying either assembly/alignment graphs (.gfa) or sequences (.fasta). When only one gfa file is supplied, variants are called within the graph file.')
    parser.add_argument("-o", "--output", dest="output", help="Prefix of the variant and alignment graph files to produce, default is \"sequence1_sequence2\"")
    parser.add_argument("-p", dest="pcutoff", type=float, default=1e-3, help="If, the probability of observing a MUM of the observed length by random change becomes larger than this cutoff the alignment is stopped (default 1e-3).")
    parser.add_argument("-t", "--threads", dest="threads", type=int, default=0, help = "The number of threads to use for the alignment.")
    parser.add_argument("-l", "--log-level", dest="loglevel", default=20, help="Log level: 10=debug 20=info (default) 30=warn 40=error 50=fatal.")
    parser.add_argument("-m", dest="minlength", type=int, default=20, help="Min length of an exact match (default 20).")
    parser.add_argument("-n", dest="minsamples", type=int, default=1, help="Only align nodes that occcur in this many samples (default 1).")
    parser.add_argument("-x", dest="maxsamples", type=int, default=None, help="Only align nodes that have maximally this many samples (default None).")
    parser.add_argument("-r", dest="reference", type=str, default=None, help="Name of the sequence that should be used as a coordinate system or reference.")
    parser.add_argument("-s", dest="targetsample", type=str, default=None, help="Only align nodes in which this sample occurs.")
    parser.add_argument("--gml", dest="gml", action="store_true", default=False, help="Produce a gml graph instead gfa.")
    parser.add_argument("--gml-max", dest="hwm", default=4000, help="Max number of nodes per graph in gml output.")
    parser.add_argument("--vg", dest="vg", action="store_true", default=False, help="Produce a gfa graph without node annotations, to ensure it's parseable by vg.")
    
    parser.add_argument("--minlen", dest="minlen", type=int, default=1, help="Output variants in a gfa graph that are larger or equal to this size (default=1).")
    parser.add_argument("--maxlen", dest="maxlen", type=int, default=1000, help="Output variants in a gfa graph that are smaller or equal to this size (to limit memory consumption for nw-alignments for inversions) (default=1000).")
    parser.add_argument("--snps", dest="snps", action="store_true", default=False, help="Output snps in a gfa graph.")
    parser.add_argument("--indels", dest="indels", action="store_true", default=False, help="Output indels in a gfa graph.")
    parser.add_argument("--inv", dest="inv", action="store_true", default=False, help="Output inversions in a gfa graph.")
    parser.add_argument("--multi", dest="multi", action="store_true", default=False, help="Output variants with multiple alleles in a gfa graph.")
    parser.add_argument("--all", dest="allvar", action="store_true", default=False, help="Output all variants in a gfa graph.")
    
    args = parser.parse_args()
    logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p', level=args.loglevel)
    
    if len(args.inputfiles)==1 and args.inputfiles[0].endswith(".gfa"):
        import caller
        logging.info("Parsing GFA file.")
        c=caller.Caller(args.inputfiles[0])
        logging.info("Writing variants.")
        c.call(minlength=args.minlen,maxlength=args.maxlen,allvar=args.allvar,inversions=args.inv,indels=args.indels,snps=args.snps,multi=args.multi)
        return
    
    if len(args.inputfiles)<=1:
        logging.fatal("Specify at least 2 (g)fa files for creating a reference graph.")
        return
    
    if args.output==None:
        args.output="_".join([os.path.basename(f)[:os.path.basename(f).rfind('.')] for f in args.inputfiles])
    
    #globale variables to simplify callbacks from c extension
    global t,G,reference
    
    reference=args.reference
    
    t=IntervalTree()
    idx=reveallib.index()
    G=nx.DiGraph()
    G.graph['samples']=[]
    schemes.pcutoff=args.pcutoff
    schemes.minlength=args.minlength
    
    for i,sample in enumerate(args.inputfiles):
        idx.addsample(sample)
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
                G.add_node(Intv,sample={sample},contig={name},coordsample=sample,coordcontig=name,start=0,aligned=0)

    print "Constructing index..."
    idx.construct()
    print "Done."
    
    print "Aligning genomes..."
    if len(args.inputfiles)>2:
        idx.align(schemes.mumpicker2,graphalign,threads=args.threads)
    else:
        idx.align(None,graphalign,threads=args.threads)
    print "Done."
    
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
