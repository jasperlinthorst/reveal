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
import time


def fasta_reader(fn,truncN=False):
    seq=""
    with open(fn,'r') as ff:
        for line in ff:
            if line.startswith(">"):
                if seq!="":
                    yield name,seq
                name=line.rstrip().replace(">","")
                seq=""
            else:
                if truncN:
                    for base in line.rstrip():
                        if base.upper()=='N':
                            if len(seq)==0:
                                continue
                            elif seq[-1]=='N':
                                continue
                            else:
                                seq+='N'
                        else:
                            seq+=base
                else:
                    seq+=line.rstrip()
        if seq!="":
            yield name,seq

def rc(seq):
    d = {'A':'T','C':'G','G':'C','T':'A','N':'N','a':'t','c':'g',\
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
    global o
    ri=0
    if reference!=None:
        for i,node in enumerate(mns):
            if reference in G.node[node]['sample']:
                refnode=node
                ri=i
                break
        else:
            refnode=mns[ri]
    else:
        refnode=mns[ri]
    
    G.node[refnode]['sample']=set.union(*[G.node[node]['sample'] for node in mns])
    G.node[refnode]['contig']=set.union(*[G.node[node]['contig'] for node in mns])
    o+=1 #increment counter, to be able to keep track of when nodes were aligned 
    G.node[refnode]['aligned']=o
    #mns.remove(refnode)
    tmp=mns.pop(ri)
    assert(tmp==refnode)
    for mn in mns: #leave the first node, remove the rest
        for e in G.in_edges(mn):
            #assert(e[0]!=refnode)
            #assert(not G.has_edge(e[0],refnode))
            #if G.has_edge(refnode,e[0]):
            #    print mns, refnode
            #    write_gml(G,T,outputfile="test.gml")
            assert(not G.has_edge(refnode,e[0]))
            G.add_edge(e[0],refnode)
        for e in G.out_edges(mn):
            #assert(e[1]!=refnode)
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
    try:
        nodes=index.nodes
        global o
        
        #print "PICKED MUM",l,n,score,sp
       
        if len(nodes)==0:
            return
        
        if l==0:
            return
        
        if l<schemes.minlength or score<schemes.minscore:
            #print "rejected",l,score
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
            if isinstance(old,Interval):
                nodes.remove((old.begin,old.end))
            for node in other:
                if isinstance(node,Interval):
                    nodes.append((node.begin,node.end))
        
        mn=mergenodes(mns)
        intervals=segmentgraph(mn,nodes)
        return intervals
    except Exception as e:
        print "Exception in graphalign",type(e),str(e)
        return None

def prune_nodes(G,T):
    converged=False
    reverse=False
    while not(converged):
        converged=True
        for node,data in G.nodes_iter(data=True):
            for run in [0,1]:
                if data['aligned']!=0:
                    if run==0:
                        reverse=False
                        neis=G.successors(node)
                    else:
                        reverse=True
                        neis=G.predecessors(node)
                    seqs={}
                    for nei in neis:
                        if G.node[nei]['aligned']==0:
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
            
            if len(s)>3:
                for v in s[4:]:
                    v=v.split(':')
                    if v[2].find(sep)!=-1: #multi-valued
                        ann[v[0]]=v[2].rstrip(sep).split(sep)
                    else:
                        ann[v[0]]=v[2]
            
            if "ORI" not in ann: #not a reveal graph, so no metadata on nodes, just index all
                if index!=None:
                    if len(s)<3:
                        intv=None
                    else:
                        intv=index.addsequence(s[2].upper())
                        intv=Interval(intv[0],intv[1])
                        tree.add(intv)
                    graph.add_node(intv,sample={gfafile},contig={nodeid},coordsample={gfafile},coordcontig={nodeid},start=0,aligned=0)
                    mapping[nodeid]=intv
                else:
                    if len(s)<3:
                        s.append("")
                    graph.add_node(nodeid,sample={gfafile},contig={nodeid},coordsample={gfafile},coordcontig={nodeid},start=0,seq=s[2].upper(),aligned=0)
                    mapping[nodeid]=nodeid
            else: #there are annotations on the nodes, so use them
                if not(isinstance(ann['ORI'], list)):
                    ann['ORI']=[ann['ORI']]
                if not(isinstance(ann['CTG'], list)):
                    ann['CTG']=[ann['CTG']]
                
                ann['ORI']=set(ann['ORI'])
                ann['CTG']=set(ann['CTG'])
                
                if len(ann['ORI'])<minsamples: #dont index these nodes, just add them to the graph
                    graph.add_node(nodeid,sample=ann['ORI'],contig=ann['CTG'],coordsample=ann['CRD'],coordcontig=ann['CRDCTG'],start=int(ann['START']),seq=s[2].upper(),aligned=0)
                    mapping[nodeid]=nodeid
                elif maxsamples!=None and len(ann['ORI'])>maxsamples: #dont index these nodes, just add them to the graph
                    graph.add_node(nodeid,sample=ann['ORI'],contig=ann['CTG'],coordsample=ann['CRD'],coordcontig=ann['CRDCTG'],start=int(ann['START']),seq=s[2].upper(),aligned=0)
                    mapping[nodeid]=nodeid
                elif (targetsample!=None and targetsample not in ann['ORI']): #dont index these nodes, just add them to the graph
                    graph.add_node(nodeid,sample=ann['ORI'],contig=ann['CTG'],coordsample=ann['CRD'],coordcontig=ann['CRDCTG'],start=int(ann['START']),seq=s[2].upper(),aligned=0)
                    mapping[nodeid]=nodeid
                else:
                    if index!=None:
                        intv=index.addsequence(s[2].upper())
                        intv=Interval(intv[0],intv[1])
                        tree.add(intv)
                        graph.add_node(intv,sample=ann['ORI'],contig=ann['CTG'],coordsample=ann['CRD'],coordcontig=ann['CRDCTG'],start=int(ann['START']),aligned=0)
                        mapping[nodeid]=intv
                    else:
                        graph.add_node(nodeid,sample=ann['ORI'],contig=ann['CTG'],coordsample=ann['CRD'],coordcontig=ann['CRDCTG'],start=int(ann['START']),seq=s[2].upper(),aligned=0)
                        mapping[nodeid]=nodeid
        
        #L      206     +       155     +       0M
        if line.startswith('L'):
            edges.append(line)

    for line in edges:
        e=line.strip().split()
        assert(not graph.has_edge(mapping[e[1]],mapping[e[3]]))
        assert(not graph.has_edge(mapping[e[3]],mapping[e[1]]))
        graph.add_edge(mapping[e[1]],mapping[e[3]],ofrom=e[2],oto=e[4],cigar=e[5])



def write_gfa(G,T,outputfile="reference.gfa", nometa=False, path=False):
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
            if isinstance(node,Interval):
                f.write('S\t'+str(i)+'\t'+T[node.begin:node.end].upper())
            else:
                f.write('S\t'+str(i)+'\t')
        
        if nometa:
            f.write("\n")
        else:
            f.write("\t*\tORI:Z:%s\tCRD:Z:%s\tCRDCTG:Z:%s\tCTG:Z:%s\tSTART:Z:%s\tALIGNED:Z:%s\n" % 
                    (
                    sep.join(data['sample']),
                    data['coordsample'].replace(sep,"").replace("\t","") if data['coordsample']!=None else "",
                    data['coordcontig'].replace(sep,"").replace("\t","") if data['coordcontig']!=None else "",
                    sep.join(data['contig']),
                    data['start'],
                    data['aligned']
                    )
                )
        
        if path:
            f.write("\n")
            for sample in data['sample']:
                f.write("P\t"+str(i)+"\t"+sample+"\t+\t"+str(node.end-node.begin)+"M\n")
        
        for to in G[node]:
            f.write('L\t'+str(mapping[node])+'\t+\t'+str(mapping[to])+"\t+\t0M\n")
    
    f.close()

def write_gml(G,T,outputfile="reference",hwm=4000):
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
            if isinstance(n,Interval):
                G.node[n]['seq']=T[n.begin:n.end].upper()
            else:
                G.node[n]['seq']=""

        G.node[n]['l']=len(G.node[n]['seq'])
    G=nx.relabel_nodes(G,mapping,copy=True)
    outputfiles=[]
    
    i=0
    for subset in nx.connected_components(G.to_undirected()):
        sgn=[]
        g=G.subgraph(subset)
        gn=len(G.graph['samples'])
        for n in nx.topological_sort(g):
            if sgn==[]:
                fr=n
            sgn.append(n)
            if G.node[n]['n']==gn: #join/split node
                if len(sgn)>=hwm:
                    sg=G.subgraph(sgn)
                    fn=outputfile+'.'+str(i)+'_'+G.node[fr]['start']+'_'+G.node[n]['start']+'.gml'
                    nx.write_gml(sg,fn)
                    outputfiles.append(fn)
                    sgn=[]
                    i+=1
        if len(sgn)>0:
            sg=G.subgraph(sgn)
            fn=outputfile+'.'+str(i)+'_'+G.node[fr]['start']+'_'+G.node[n]['start']+'.gml'
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
    parser_convert = subparsers.add_parser('convert', prog="reveal convert", description="Convert gfa graph to gml.")
    parser_extract = subparsers.add_parser('extract', prog="reveal extract", description="Extract the input sequence from a graph.")
    parser_comp = subparsers.add_parser('comp', prog="reveal comp", description="Reverse complement the graph.")
    parser_subgraph = subparsers.add_parser('subgraph', prog="reveal subgraph", description="Extract subgraph from gfa by specified node ids.")
    parser_bubbles = subparsers.add_parser('bubbles', prog="reveal bubbles", description="Extract all bubbles from the graph.")
    
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
    parser_aln.add_argument("--nometa", dest="nometa", action="store_true", default=False, help="Produce a gfa graph without node annotations, to ensure it's parseable by other programs.")
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
    parser_plot.add_argument("-i", dest="interactive", action="store_true", default=False, help="Wheter to produce interactive plots which allow zooming on the dotplot (default=False).")
    #parser_plot.add_argument("", dest="pos", default=None, help="Position on the reference genome to vizualize.")
    #parser_plot.add_argument("", dest="env", default=100, help="Size of the region aroung the targeted position to vizualize.")
    parser_plot.set_defaults(func=plot)
    
    parser_comp.add_argument('graph', nargs=1, help='The graph to be reverse complemented.')
    parser_comp.set_defaults(func=comp)

    parser_convert.add_argument('graph', nargs=1, help='The gfa graph to convert to gml.')
    parser_convert.add_argument("-n", dest="minsamples", type=int, default=1, help="Only align nodes that occcur in this many samples (default 1).")
    parser_convert.add_argument("-x", dest="maxsamples", type=int, default=None, help="Only align nodes that have maximally this many samples (default None).")
    parser_convert.add_argument("-s", dest="targetsample", type=str, default=None, help="Only align nodes in which this sample occurs.")
    parser_convert.add_argument("--gml-max", dest="hwm", default=4000, help="Max number of nodes per graph in gml output.")
    parser_convert.set_defaults(func=convert)

    parser_subgraph.add_argument('inputfiles', nargs='*', help='The gfa graph followed by node ids that make up the subgraph.')
    parser_subgraph.add_argument("-o", dest="outfile", type=str, default="~tmp", help="Prefix of the file to which subgraph will be written.")
    parser_subgraph.add_argument("--gml", dest="gml", action="store_true", default=False, help="Produce a gml graph instead of gfa.")
    parser_subgraph.set_defaults(func=subgraph)

    parser_bubbles.add_argument("graph", nargs=1, help='Graph in gfa format from which bubbles are to be extracted.')
    parser_bubbles.set_defaults(func=bubbles)

    args = parser.parse_args()
    logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p', level=args.loglevel)
    args.func(args)

def bubbles(args):

    if len(args.graph)<1:
        logging.fatal("Specify a gfa file to extract bubbles.")
        return

    G=nx.DiGraph()
    read_gfa(args.graph[0],None,"",G)
    
    #G.add_edge(1,2)
    #G.add_edge(1,3)
    #G.add_edge(2,3)
    #G.add_edge(3,4)
    #G.add_edge(4,8)
    #G.add_edge(3,5)
    #G.add_edge(5,6)
    #G.add_edge(5,9)
    #G.add_edge(9,10)
    #G.add_edge(6,10)
    #G.add_edge(10,7)
    #G.add_edge(6,7)
    #G.add_edge(7,8)
    #G.add_edge(3,11)
    #G.add_edge(11,12)
    #G.add_edge(12,8)
    #G.add_edge(8,13)
    #G.add_edge(8,14)
    #G.add_edge(13,14)
    #G.add_edge(13,15)
    #G.add_edge(15,14)
    
    def entrance(G,v):
        for c in G.successors(v):
            if len(G.predecessors(c))==1:
                return True
        return False 

    def exit(G,v):
        for p in G.predecessors(v):
            if len(G.successors(p))==1:
                return True
        return False
    
    def nextentrance(candidates,v):
        #TODO: rewrite this
        for candidate in candidates[candidates.index((v,0))+1:]:
            if candidate[1]==0:
                return candidate
    
    def superbubble(G):
        candidates=[]
        sspairs=[]
        #prevEnt=None
        prevEnti=None
        alternativeEntrance={}
        previousEntrance={}
        ordD={}
        ordD_=nx.topological_sort(G)
        
        #construct candidates array
        for i,v in enumerate(ordD_):
            ordD[v]=i
            alternativeEntrance[v]=None
            previousEntrance[v]=prevEnti
            if exit(G,v):
                candidates.append((v,1))
            if entrance(G,v):
                candidates.append((v,0))
                #prevEnt=v
                prevEnti=i
        
        #construct outparent
        outparent=[None]*(len(ordD))
        for i,c in enumerate(ordD):
            tmp=[]
            for p in G.predecessors(c):
                tmp.append(ordD[p])
            if len(tmp)>0:
                outparent[ordD[c]]=min(tmp)
        
        #construct outchild
        outchild=[None]*(len(ordD))
        for i,c in enumerate(ordD):
            tmp=[]
            for p in G.successors(c):
                tmp.append(ordD[p])
            if len(tmp)>0:
                outchild[ordD[c]]=max(tmp)
        
        #loop
        while len(candidates)!=0:
            if candidates[-1][1]==0:
                del candidates[-1]
            else:
                reportsuperbubble(candidates[0],candidates[-1],candidates,previousEntrance,alternativeEntrance,G,ordD,ordD_,outchild,outparent,sspairs)
        
        return ordD,ordD_,sspairs
    
    def reportsuperbubble(vstart,vexit,candidates,previousEntrance,alternativeEntrance,G,ordD,ordD_,outchild,outparent,sspairs):
        if (vstart[0] == None) or (vexit[0] == None) or (ordD[vstart[0]] >= ordD[vexit[0]]):
            del candidates[-1]
            return
        si=previousEntrance[vexit[0]]
        s=ordD_[si]
        while ordD[s] >= ordD[vstart[0]]:
            valid = validatesuperbubble(s, vexit[0], ordD, ordD_, outchild, outparent, previousEntrance, G)
            if valid==s or valid==alternativeEntrance[s] or valid==-1:
                break
            alternativeEntrance[s] = valid
            s = valid

        del candidates[-1]
        if (valid == s):
            sspairs.append((s, vexit[0]))
            while (candidates[-1][0] is not s):
                if candidates[-1][1]==1:
                    ne=nextentrance(candidates,s)
                    if ne!=None:
                        reportsuperbubble(ne, candidates[-1], candidates, previousEntrance, alternativeEntrance, G, ordD, ordD_, outchild, outparent,sspairs)
                    else:
                        del candidates[-1]
                else:
                    del candidates[-1]
    
    def validatesuperbubble(startVertex, endVertex, ordD, ordD_, outchild, outparent, previousEntrance, G):
        start=ordD[startVertex]
        end=ordD[endVertex]
        if start+1!=end:
            oc=max(outchild[start:end])
            op=min(outparent[start+1:end+1])
        else:
            oc=outchild[start]
            op=outparent[end]
        if oc!=end:
            return -1
        if op==start:
            return startVertex
        elif entrance(G, ordD_[op]):
            return ordD_[op]
        elif previousEntrance[ordD_[op]]==None: #
            return -1
        else:
            #print ordD_[op]
            #print previousEntrance[ordD_[op]]
            return ordD_[previousEntrance[ordD_[op]]]
        return startVertex
    
    ordD,ordD_,sspairs=superbubble(G)
    
    if 'samples' in G.graph:
        gori=sorted(G.graph['samples'])
    else:
        gori=[]
    
    sg=set()
    for pair in sspairs:
        size=(ordD[pair[1]]-ordD[pair[0]])-1
        bubblenodes=ordD_[ordD[pair[0]]:ordD[pair[1]]+1]

        source=G.node[pair[0]]['sample']
        sink=G.node[pair[1]]['sample']
        
        #TODO: supbub algorithm detects invalid bubbles; filter them out here
        if source!=sink:
            #print "ERROR, skipping invalid bubble"
            continue
        
        #determine genotypes...
        gt=G.successors(pair[0])
        gt.sort(key=lambda l: ordD[l]) #topological sort neighbors
        
        if len(gt)<size:
            #complex bubble; not possible to output the actual alleles, instead we can make a call based on branching at source node and output the nodes that make up the complex bubble so we can export and vizualize the subgraph
            genotypes=['N']*len(gt)
            calls={}
            for i,node in enumerate(gt):
                n=G.node[node]
                for sample in n['sample']:
                    calls[sample]=i
        else:
            genotypes=[None]*len(gt)
            calls={}
            for i,node in enumerate(gt):
                n=G.node[node]
                if node==pair[1]: #should always be last of the iteration
                    genotypes[i]="-"
                    for sample in source: #all thats left..
                        calls[sample]=i #TODO: need sample annotations on edges to properly make this call! For now, use this... works in case of monoploids
                else:
                    genotypes[i]=n['seq']
                    for sample in n['sample']:
                        if sample in calls:
                            print "WARNING, multiple paths (diploid variant?), report only one allele..."
                        else:
                            calls[sample]=i
        
        d=G.node[pair[0]]
        if 'coordcontig' in d:
            crdctg=d['coordcontig']
        else:
            crdctg='N/A'
        
        if 'coordsample' in d:
            crdsmpl=d['coordsample']
        else:
            crdsmpl='N/A'
        
        if 'start' in d:
            pos=d['start']+len(d['seq'])+1
        else:
            pos='N/A'
        
        sys.stdout.write("%s\t%s\t%s\t%s\t%s\t%s\t%s"%(pair[0],pair[1],",".join(bubblenodes),crdsmpl,crdctg,pos,",".join(genotypes)))
         
        for sample in gori:
            if sample in calls:
                sys.stdout.write("\t%s"%calls[sample])
            else:
                sys.stdout.write("\t-")
        sys.stdout.write("\n")

def subgraph(args):
    if len(args.inputfiles)<=1:
        logging.fatal("Specify 1 gfa file followed by node ids for which subgraph is to be extracted.")
        return
    if not args.inputfiles[0].endswith('.gfa'):
        logging.fatal("Specify gfa file as first argument of subgraph subcommand.")
        return
    G=nx.DiGraph()
    read_gfa(args.inputfiles[0],None,"",G)
    nodes=set()
    for node in args.inputfiles[1:]:
        nodes.add(node)
    sg=G.subgraph(nodes)
    if args.gml:
        write_gml(sg,"",outputfile=args.outfile)
    else:
        write_gfa(sg,"",outputfile=args.outfile+".gfa")

def align_seq(s1,s2,minlength=1,minscore=0):
    global t,G,reference,o

    t=IntervalTree()
    idx=reveallib.index()
    G=nx.DiGraph()
    G.graph['samples']=[]
    reference=None
    o=0
    
    idx.addsample("s1")
    intv=idx.addsequence(s1.upper())
    Intv=Interval(intv[0],intv[1])
    t.add(Intv)
    G.add_node(Intv,sample={"s1"},contig={"s1"},coordsample="s1",coordcontig="s1",start=0,aligned=0)

    idx.addsample("s2")
    intv=idx.addsequence(s2.upper())
    Intv=Interval(intv[0],intv[1])
    t.add(Intv)
    G.add_node(Intv,sample={"s2"},contig={"s2"},coordsample="s2",coordcontig="s2",start=0,aligned=0)
    
    schemes.ts=t
    schemes.minlength=minlength
    schemes.minscore=minscore

    idx.construct()
    
    idx.align(None,graphalign)
    
    alignedbases=0
    totbases=min([idx.n-idx.nsep[0],idx.nsep[0]])
    for node,data in G.nodes(data=True):
        if data['aligned']!=0:
            alignedbases+=(node.end-node.begin)
    
    return alignedbases/float(totbases)

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
    
    logging.info("Merging nodes...")
    T=idx.T
    if len(G.graph['samples'])>2:
        prune_nodes(G,T) #TODO: do this after every alignment step to reduce graph size during alignment
    logging.info("Done.")
    
    logging.info("Writing graph...")
    if args.gml:
        graph=write_gml(G,T, hwm=args.hwm, outputfile=args.output)
    else:
        write_gfa(G,T,nometa=args.nometa, outputfile=args.output+'.gfa')
        graph=args.output+'.gfa'
    logging.info("Done.")
    
    logging.info("Alignment graph written to: %s"%graph)

def align_genomes(args):
    #globale variables to simplify callbacks from c extension
    global t,G,reference,o
    
    reference=args.reference
    
    t=IntervalTree()
    idx=reveallib.index()
    G=nx.DiGraph()
    G.graph['samples']=[]
    o=0
    #schemes.pcutoff=args.pcutoff
    schemes.minlength=args.minlength
    schemes.minscore=args.minscore
    
    for i,sample in enumerate(args.inputfiles):
        idx.addsample(os.path.basename(sample))
        if sample.endswith(".gfa"):
            #TODO: now applies to all graphs! probably want to have this graph specific if at all...
            logging.info("Reading graph: %s ..." % sample)
            read_gfa(sample,idx,t,G,minsamples=args.minsamples,
                                    maxsamples=args.maxsamples,
                                    targetsample=args.targetsample)
            logging.info("Done.")
        else: #consider it to be a fasta file
            G.graph['samples'].append(os.path.basename(sample))
            for name,seq in fasta_reader(sample):
                intv=idx.addsequence(seq.upper())
                Intv=Interval(intv[0],intv[1])
                t.add(Intv)
                #schemes.ts.add(Intv)
                schemes.interval2sampleid[Intv]=i
                G.add_node(Intv,sample={os.path.basename(sample)},contig={name.replace(";","")},coordsample=os.path.basename(sample),coordcontig=name.replace(";",""),start=0,aligned=0)
    
    if not nx.is_directed_acyclic_graph(G):
        logging.error("*** Input is not a DAG! Experimental output...")
    
    schemes.ts=t
    
    logging.info("Constructing index...")
    idx.construct()
    logging.info("Done.")

    if len(args.inputfiles)>2:
        logging.info("Constructing multi-alignment...")
        idx.align(schemes.multimumpicker,graphalign,threads=args.threads)
    else:
        logging.info("Constructing pairwise-alignment...")
        idx.align(None,graphalign,threads=args.threads)
    
    alignedbases=0
    if idx.nsamples>2:
        totbases=idx.n-nx.number_of_nodes(G)
        for node,data in G.nodes(data=True):
            if data['aligned']!=0:
                alignedbases+=(node.end-node.begin)*len(data['sample'])
    else:
        totbases=min([idx.n-idx.nsep[0],idx.nsep[0]])
        for node,data in G.nodes(data=True):
            if data['aligned']!=0:
                alignedbases+=(node.end-node.begin)
    
    logging.info("Done (%.2f%% identity, %d bases out of %d aligned)."%((alignedbases/float(totbases))*100,alignedbases,totbases))
    
    return G,idx

def align_contigs(args):
    global t,G,reference,o
    
    reference=args.reference
    schemes.minlength=args.minlength
    schemes.minscore=args.minscore
    G=nx.DiGraph()
    G.graph['samples']=[]
    
    ref=args.inputfiles[0]
    contigs=args.inputfiles[1]
    
    t=IntervalTree()
    idx=reveallib.index()
    
    totbases=0
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
    
    totbases=idx.n-nx.number_of_nodes(G)
    alignedbases=0
    i=0
    G.graph['samples'].append(contigs)
    for contigname,contig in fasta_reader(contigs):
        totbases+=len(contig)
        idx.addsample(os.path.basename(contigs))
        intv=idx.addsequence(contig.upper())
        Intv=Interval(intv[0],intv[1])
        t.add(Intv)
        G.add_node(Intv,sample={contigs},contig={contigname.replace(";","")},coordsample=contigs,coordcontig=contigname.replace(";",""),start=0,aligned=0)

    logging.info("Constructing index... (size=%d)"%idx.n)
    idx.construct()
    logging.info("Done.")

    logging.info("Aligning %s (%d bp)" % (contigname,len(contig)))
    idx.align(None,graphalign,threads=args.threads)
    identity=0
    for node,data in G.nodes(data=True):
        if data['aligned']!=0 and isinstance(node,Interval):
            identity+=node.end-node.begin

    alignedbases+=identity
    logging.info("%s aligned, %d bp out of %d bp aligned (%.2f%%)." % (contigname,identity,len(contig),100*(identity/float(len(contig)))))
    
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
                if len(nx.all_simple_paths(G, source=neis[0], target=neis[1]))==1 and G.node[neis[0]]['aligned']!=0 and G.node[neis[1]]['aligned']!=0:
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
    logging.info("Done (%.2f%% identity, %d bases out of %d aligned)."%(((alignedbases*2)/float(totbases))*100,alignedbases,totbases))
    return G,idx

#extract variants from gfa graphs
def call(args):
    if args.graph[0].endswith(".gfa"):
        import caller
        logging.info("Parsing GFA file.")
        c=caller.Caller(args.graph[0],logger=logging)
        logging.info("Writing variants.")
        c.call(minlength=args.minlen,maxlength=args.maxlen,allvar=args.allvar,inversions=args.inv,indels=args.indels,snps=args.snps,multi=args.multi)

def plot(args):
    from matplotlib import pyplot as plt
    
    if len(args.fastas)==2:    
        #get mmems for forward orientation
        idx=reveallib.index()
        for sample in args.fastas:
            idx.addsample(sample)
            for name,seq in fasta_reader(sample,truncN=False):
                intv=idx.addsequence(seq.upper())
                break #expect only one sequence per fasta for now
        idx.construct()
        mmems=[(mem[0],mem[1],mem[2],0) for mem in idx.getmums() if mem[0]>args.minlength]
        sep=idx.nsep[0]
        
        #get mmems for reverse orientation
        idx=reveallib.index()
        sample=args.fastas[0]
        idx.addsample(sample)
        for name,seq in fasta_reader(sample,truncN=False):
            horzgaps=[i for i,c in enumerate(seq) if c=='N']
            intv=idx.addsequence(seq.upper())
            break #expect only one sequence per fasta for now
        
        sample=args.fastas[1]
        idx.addsample(sample)
        ls=0
        for name,seq in fasta_reader(sample,truncN=False):
            vertgaps=[i for i,c in enumerate(seq) if c=='N']
            ls+=len(seq)
            intv=idx.addsequence(rc(seq.upper()))
            break #expect only one sequence per fasta for now
        idx.construct()
        mmems+=[(mem[0],mem[1],mem[2],1) for mem in idx.getmums() if mem[0]>args.minlength]
        
    elif len(args.fastas)==1:

        idx=reveallib.index()
        sample=args.fastas[0]
        idx.addsample(sample)
        for name,seq in fasta_reader(sample, truncN=False):
            horzgaps=[i for i,c in enumerate(seq) if c=='N']
            intv=idx.addsequence(seq.upper())
            break #expect only one sequence per fasta for now
        
        sample=args.fastas[0]+"_rc"
        idx.addsample(sample)
        ls=0
        for name,seq in fasta_reader(sample, truncN=False):
            vertgaps=[i for i,c in enumerate(seq) if c=='N']
            ls+=len(seq)
            intv=idx.addsequence(rc(seq.upper()))
            break #expect only one sequence per fasta for now
        idx.construct()
        mmems=[(mem[0],mem[1],mem[2],0) for mem in idx.getmums() if mem[0]>args.minlength]
    else:
        logging.fatal("Can only create mumplot for 2 sequences or self plot for 1 sequence.")
        return
    
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
    
    lgap=0
    g=0
    for gap in horzgaps:
        if lgap==gap-1:
            lgap=gap
            continue
        plt.axvline(gap)
        g+=1
        lgap=gap
    
    g=0
    lgap=0
    for gap in vertgaps:
        if lgap==gap-1:
            lgap=gap
            continue
        plt.axhline(gap)
        g+=1
        lgap=gap

    plt.title(args.fastas[0]+" vs. "+args.fastas[1])
    plt.xlabel(args.fastas[0])
    plt.ylabel(args.fastas[1])
    plt.autoscale(enable=False)
    
    if args.interactive:
        plt.show()
    else:
        fn1=args.fastas[0][0:args.fastas[0].rfind('.')] if args.fastas[0].find('.')!=-1 else args.fastas[0]
        fn2=args.fastas[1][0:args.fastas[1].rfind('.')] if args.fastas[1].find('.')!=-1 else args.fastas[1]
        plt.savefig(fn1+"_"+fn2+".png")

def comp(args):
    G=nx.DiGraph()
    G.graph['samples']=[]
    t=IntervalTree()
    read_gfa(args.graph[0],None,t,G,targetsample=None)
    for node in G.node:
        G.node[node]['seq']=rc(G.node[node]['seq'])
    G.reverse(copy=False)
    write_gfa(G,"",outputfile=args.graph[0].replace('.gfa','.rc.gfa'), nometa=True)

def convert(args):
    if not args.graph[0].endswith(".gfa") or len(args.graph)!=1:
        logging.fatal("Specify a gfa file.")
        return
    
    G=nx.DiGraph()
    G.graph['samples']=[]
    t=IntervalTree()
    
    read_gfa(args.graph[0],None,t,G,minsamples=args.minsamples,
                         maxsamples=args.maxsamples,
                         targetsample=args.targetsample)
    
    graph=write_gml(G,"", hwm=args.hwm, outputfile=args.graph[0].rstrip(".gfa"))
    logging.info("GML graph written to: %s"%graph)

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
