#!/usr/bin/env python

from intervaltree import Interval, IntervalTree
import networkx as nx
from collections import defaultdict, deque
import threading
import reveallib
import reveallib64
import argparse
import logging
import os
import schemes
import sys
import time

try:
    from matplotlib import pyplot as plt
except:
    pass

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

def fasta_writer(fn,name_seq,lw=100):
    seq=""
    with open(fn,'w') as ff:
        for name,seq in name_seq:
            if not name.startswith(">"):
                name=">"+name+"\n"
            ff.write(name)
            for i in range( (len(seq)/lw)+(len(seq) % lw > 0)):
                ff.write(seq[i*lw:(i+1)*lw]+"\n")

def rc(seq):
    d = {'A':'T','C':'G','G':'C','T':'A','N':'N','a':'t','c':'g',\
        'g':'c','t':'a','n':'n','Y':'R','R':'Y','K':'M','M':'K',\
        'S':'S','W':'W','B':'V','V':'B','D':'H','H':'D','N':'N',\
        'X':'X','-':'-'}
    return "".join([d[b] for b in reversed(seq)])

def breaknode(node,pos,l):
    logging.debug("Breaking node: %s"%str(node))
    att=G.node[node]
    in_edges=G.in_edges(node)
    out_edges=G.out_edges(node)
    mn=Interval(pos,pos+l)
    other=set()
    
    if mn==node: #no breaking needed
        t.remove(node)
        return node,other
    
    moffsets=dict()
    for s in att['offsets']:
        moffsets[s]=att['offsets'][s]+(pos-node.begin)
    
    soffsets=dict()
    for s in att['offsets']:
        soffsets[s]=att['offsets'][s]+((pos+l)-node.begin)
    
    G.add_node(mn,sample=att['sample'],offsets=moffsets,aligned=0)#create merge node
    
    if (node[0]!=pos):
        pn=Interval(node[0],pos)
        G.add_node(pn,sample=att['sample'],offsets=att['offsets'],aligned=0)#create prefix node
        assert(not G.has_edge(pn,mn))
        assert(not G.has_edge(mn,pn))
        G.add_edge(pn,mn)
        t.add(pn)
        other.add(pn)
    else:
        pn=mn
    if (node[1]!=pos+l):
        sn=Interval(pos+l,node[1])
        G.add_node(sn,sample=att['sample'],offsets=soffsets,aligned=0)#create suffix node
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
    
    logging.debug("Leading/Trailing node(s): %s"%str(other))
    logging.debug("Matching node: %s"%str(mn))

    return mn,other #return merge node

def mergenodes(mns,mark=True):
    logging.debug("Merging nodes %s"%str(mns))
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
    
    #create offset dict for gap penalty calculation from alignment graph
    newoffsets=dict()
    for node in mns:
        d=G.node[node]
        for sample in d['offsets']:
            if sample in newoffsets:
                print "Error, merging nodes that originate from the same sample."
            assert(sample not in newoffsets)
            newoffsets[sample]=d['offsets'][sample]
    
    G.node[refnode]['offsets']=newoffsets
    G.node[refnode]['sample']=set.union(*[G.node[node]['sample'] for node in mns])
    assert(len(G.node[refnode]['sample'])==len(newoffsets))

    if mark:
        o+=1 #increment counter, to be able to keep track of when nodes were aligned 
        G.node[refnode]['aligned']=o

    #mns.remove(refnode)
    tmp=mns.pop(ri)
    assert(tmp==refnode)
    
    for mn in mns: #leave the first node, remove the rest
        for e in G.in_edges(mn):
            assert(not G.has_edge(refnode,e[0]))
            G.add_edge(e[0],refnode)
        for e in G.out_edges(mn):
            assert(not G.has_edge(e[1],refnode))
            G.add_edge(refnode,e[1])
        G.remove_node(mn)

    return refnode

def bfs(G, source, reverse=False, ignore=set()):
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
                elif (G.node[child]['aligned']!=0 and child in ignore): #keep searching
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
    trace=False
    
    #forward search
    endpoints=set()
    for c,t in bfs(G,node):
        if t==0:
            trailing.add(c)
        else:
            endpoints.add(c)
    
    #reverse search for each endpoint
    if len(endpoints)>1:
        for endpoint in endpoints:
            for c,t in bfs(G,endpoint,reverse=True,ignore=endpoints):
                if t==0:
                    reverse_trailing.add(c)
        trailing=trailing.intersection(reverse_trailing)
    
    #backward search
    endpoints=set()
    for c,t in bfs(G,node,reverse=True):
        if t==0:
            leading.add(c)
        else:
            endpoints.add(c)
    
    #reverse search for each endpoint
    if len(endpoints)>1:
        for endpoint in endpoints:
            for c,t in bfs(G,endpoint,ignore=endpoints):
                if t==0:
                    reverse_leading.add(c)
        leading=leading.intersection(reverse_leading)
    
    leading = set([(i.begin,i.end) for i in leading if isinstance(i,Interval)]).intersection(nodes) #TODO: remove "if isinstance(i,Interval)]"
    trailing = set([(i.begin,i.end) for i in trailing if isinstance(i,Interval)]).intersection(nodes)
    
    rest = nodes - (leading | trailing)
    
    return list(leading), list(trailing), list(rest), (node.begin,node.end)

def graphalign(l,index,n,score,sp,penalty):
    
    nodes=index.nodes
    isize=index.n
    
    if len(nodes)==0:
        logging.debug("Invalid set of nodes (length=%d, samples=%d, score=%d, penalty=%d, sp=%s, indexsize=%d)",l,n,score,penalty,sp,isize)
        return
    
    if l==0:
        logging.debug("Invalid length (length=%d, samples=%d, score=%d, penalty=%d, sp=%s, indexsize=%d)",l,n,score,penalty,sp,isize)
        return
    
    if score<schemes.minscore:
        logging.debug("Reject MUM, score too low (length=%d, samples=%d, score=%d, penalty=%d, sp=%s, indexsize=%d)",l,n,score,penalty,sp,isize)
        return
    
    if l<schemes.minlength:
        logging.debug("Reject MUM, too short (length=%d, samples=%d, score=%d, penalty=%d, sp=%s, indexsize=%d)",l,n,score,penalty,sp,isize)
        return
    
    logging.debug("Align graph to MUM of length %d (samples=%d, score=%d, penalty=%d, sp=%s, indexsize=%d)",l,n,score,penalty,sp,isize)
    
    mns=[]
    topop=[]
    
    for i,pos in enumerate(sp):
        old=t[pos].pop()
        assert(old.end-old.begin>=l)
        mn,other=breaknode(old,pos,l)
        mns.append(mn)
        if isinstance(old,Interval):
            nodes.remove((old.begin,old.end))
        for node in other:
            if isinstance(node,Interval):
                nodes.append((node.begin,node.end))
    mn=mergenodes(mns)
    
    msamples=set(G.node[Interval(mn[0],mn[1])]['offsets'].keys())
    
    intervals=segmentgraph(mn,nodes)
    
    leading,trailing,rest,merged=intervals
    
    logging.debug("leading intervals: %s"%leading)
    logging.debug("trailing intervals: %s"%trailing)
    #logging.debug("parallel intervals: %s"%rest)
    
    newleft=mn
    newright=mn
    
    for intv in leading:
        if not G.node[Interval(intv[0],intv[1])]['sample'].issubset(msamples): #no clean dissection of all paths on the left
            newright=index.right
            break
    
    for intv in trailing:
        if not G.node[Interval(intv[0],intv[1])]['sample'].issubset(msamples): #no clean dissection of all paths on the right
            newleft=index.left
            break
    
    return leading,trailing,rest,merged,newleft,newright

def prune(node,T,reverse=False):
    seqs={}
    pruned=[]
    if reverse:
        neis=G.predecessors(node)
    else:
        neis=G.successors(node)
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
        tmpgroup=list(group)
        if len(group)>1:
            sink=[]
            merge=True
            for v in group:
                if reverse:
                    if len(G.predecessors(v))>1:
                        merge=False
                        break
                else:
                    if len(G.successors(v))>1:
                        merge=False
                        break
            if merge:
                merged=mergenodes(group)
                converged=False
                pruned+=group
    return pruned

def prune_nodes(G,T):
    trace=False
    converged=False
    while not(converged):
        converged=True
        for node,data in G.nodes_iter(data=True):
            if node not in G:
                continue
            
            for run in [0,1]:
                if data['aligned']!=0:
                    if run==0:
                        neis=G.successors(node)
                    else:
                        neis=G.predecessors(node)
                        if trace:
                            print "predecessors",neis
                    seqs={}
                    for nei in neis:
                        if 'seq' not in G.node[nei]:
                            if not isinstance(nei,Interval):
                                continue
                            seq=T[nei.begin:nei.end]
                        else:
                            seq=G.node[nei]['seq']
                        if seq in seqs:
                            seqs[seq].append(nei)
                        else:
                            seqs[seq]=[nei]
                    
                    if trace:
                        print "seqs",seqs
                    
                    for key in seqs.keys():
                        group=seqs[key]
                        if len(group)>1:
                            merge=True
                            for v in group:
                                if run==0:
                                    if len(G.predecessors(v))>1:
                                        merge=False
                                        break
                                else:
                                    if len(G.successors(v))>1:
                                        merge=False
                                        break
                            if merge:
                                if trace:
                                    print "mergenode",group
                                mergenodes(group,mark=True)
                                converged=False

def read_gfa(gfafile, index, tree, graph, minsamples=1, maxsamples=None, targetsample=None, revcomp=False):
    f=open(gfafile,'r')
    sep=";"
    mapping={} #temp mapping object
    edges=[] #tmp list for edges
    id2sample={}
    i=0
    for line in f:
        if line.startswith('H'):
            h=line.split()
            tag=h[1].split(':')
            if tag[0]=="ORI":
                for sample in tag[2].rstrip(sep).split(sep):
                    id2sample[i]=sample
                    i+=1
                    if 'samples' in graph.graph:
                        graph.graph['samples'].append(sample)
                    else:
                        graph.graph['samples']=[sample]
        
        if line.startswith('S'):
            s=line.strip().split('\t')
            nodeid=int(s[1])
            gnodeid=graph.number_of_nodes()+1
            
            if graph.has_node(gnodeid):
                logging.fatal("Id space for nodes is larger than total number of nodes in the graph.")
                return
            
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
                        if revcomp:
                            intv=index.addsequence(rc(s[2]).upper())
                        else:
                            intv=index.addsequence(s[2].upper())
                        intv=Interval(intv[0],intv[1])
                        tree.add(intv)
                    graph.add_node(intv,sample={gfafile},offsets={gfafile:0},aligned=0)
                    mapping[nodeid]=intv
                else:
                    if len(s)<3:
                        s.append("")
                    graph.add_node(gnodeid,sample={gfafile},offsets={gfafile:0},seq=s[2].upper(),aligned=0)
                    mapping[nodeid]=gnodeid

            else: #there are annotations on the nodes, so use them
                if not(isinstance(ann['ORI'], list)):
                    ann['ORI']=[ann['ORI']]
                if not(isinstance(ann['OFFSETS'], list)):
                    ann['OFFSETS']=[ann['OFFSETS']]
                offsets=dict()
                if 'OFFSETS' in ann: #create dictionary by sample name
                    for sampleid,offset in zip(ann['ORI'],ann['OFFSETS']):
                        offsets[id2sample[int(sampleid)]]=int(offset)
                tmp=set()
                for ori in ann['ORI']: #map sample ids back to sample names
                    tmp.add(id2sample[int(ori)])
                ann['ORI']=tmp
                
                assert(isinstance(ann['ORI'],set))
                
                if len(ann['ORI'])<minsamples: #dont index these nodes, just add them to the graph
                    graph.add_node(gnodeid,sample=ann['ORI'],seq=s[2].upper(),aligned=0,offsets=offsets)
                    mapping[nodeid]=gnodeid
                elif maxsamples!=None and len(ann['ORI'])>maxsamples: #dont index these nodes, just add them to the graph
                    graph.add_node(gnodeid,sample=ann['ORI'],seq=s[2].upper(),aligned=0,offsets=offsets)
                    mapping[nodeid]=gnodeid
                elif (targetsample!=None and targetsample not in ann['ORI']): #dont index these nodes, just add them to the graph
                    graph.add_node(gnodeid,sample=ann['ORI'],seq=s[2].upper(),aligned=0,offsets=offsets)
                    mapping[nodeid]=gnodeid
                else:
                    if index!=None:
                        if revcomp:
                            intv=index.addsequence(rc(s[2]).upper())
                        else:
                            intv=index.addsequence(s[2].upper())
                        intv=Interval(intv[0],intv[1])
                        tree.add(intv)
                        graph.add_node(intv,sample=ann['ORI'],aligned=0,offsets=offsets)
                        mapping[nodeid]=intv
                    else:
                        graph.add_node(gnodeid,sample=ann['ORI'],seq=s[2].upper(),aligned=0,offsets=offsets)
                        mapping[nodeid]=gnodeid
        
        #L      206     +       155     +       0M
        if line.startswith('L'):
            edges.append(line)

    for line in edges:
        e=line.strip().split()
        assert(not graph.has_edge(mapping[int(e[1])],mapping[int(e[3])]))
        assert(not graph.has_edge(mapping[int(e[3])],mapping[int(e[1])]))
        graph.add_edge(mapping[int(e[1])],mapping[int(e[3])],ofrom=e[2],oto=e[4],cigar=e[5])
    
    if revcomp:
        genome2length=dict()
        #relabel the offsets, determine the length of all genomes in the graph, then l-pos
        for sample in graph.graph['samples']:
            maxp=0
            for node,data in graph.nodes(data=True):
                if sample in data['offsets']:
                    if data['offsets'][sample]+ (node[1]-node[0]) >maxp:
                        maxp=data['offsets'][sample]+(node[1]-node[0])
            genome2length[sample]=maxp
        
        for sample in graph.graph['samples']:
            for node,data in graph.nodes(data=True):
                if sample in data['offsets']:
                    graph.node[node]['offsets'][sample]=genome2length[sample]-(graph.node[node]['offsets'][sample]+(node[1]-node[0]))
        
        graph.reverse(copy=False)

#simply write sequence without the graph topology
def write_fasta(G,T,outputfile="reference.fasta"):
    f=open(outputfile,'wb')
    for i,node in enumerate(nx.topological_sort(G)):
        if isinstance(node,str):
            if node=='start' or node=='end':
                continue
        i+=1
        data=G.node[node]
        seq=""
        if len(node)==3:
            nodename=node[2]
        else:
            nodename=str(node)
        if 'seq' in data:
            f.write(">%s\n"%nodename)
            f.write(data['seq'].upper()+"\n")
        else:
            if isinstance(node,Interval):
                f.write(">%s\n"%nodename)
                f.write(T[node.begin:node.end].upper()+"\n")
            else:
                f.write(">%s\n"%nodename)
                logging.warn("No sequence for node: %s"%nodename)
    f.close()

def write_gfa(G,T,outputfile="reference.gfa", nometa=False, paths=False):
    
    if not outputfile.endswith(".gfa"):
        outputfile+=".gfa"
    
    f=open(outputfile,'wb')
    sep=';'
    f.write('H\tVN:Z:1.0\n')
    f.write('H\tORI:Z:')
    
    sample2id=dict()
    i=0
    
    for sample in G.graph['samples']:
        for subsample in sample.rstrip().rstrip(sep).split(sep):
            f.write(subsample+sep)
            sample2id[subsample]=i
            i+=1
    
    f.write('\n')
    
    mapping={}
    
    for i,node in enumerate(nx.topological_sort(G)): #iterate once to get a mapping of ids to intervals
        mapping[node]=i+1 #one-based for vg
    
    for i,node in enumerate(nx.topological_sort(G)):
        if isinstance(node,str):
            if node=='start' or node=='end':
                continue

        i+=1
        data=G.node[node]
        seq=""
        if 'seq' in data:
            f.write('S\t'+str(i)+'\t'+data['seq'].upper())
            seq=data['seq']
        else:
            if isinstance(node,Interval):
                f.write('S\t'+str(i)+'\t'+T[node.begin:node.end].upper())
                seq=T[node.begin:node.end].upper()
            else:
                f.write('S\t'+str(i)+'\t')
        
        #make a list out of samples and sort by id
        if 'sample' not in data:
            print "ERROR",node,data
        
        tmp=list(data['sample'])
        tmp.sort(key=lambda s: sample2id[s])
        data['sample']=tmp
        
        if not isinstance(data['offsets'],dict):
            print "ERROR offset data:", data['offsets']
        
        if nometa:
            f.write("\n")
        else:
            f.write("\t*\tORI:Z:%s\tOFFSETS:Z:%s\tRC:i:%s\n" % 
                    (
                    sep.join([str(sample2id[s]) for s in data['sample']]),
                    sep.join([str(data['offsets'][s]) for s in data['sample']]),
                    len(data['sample'])*len(seq)
                    )
                )
        
        for to in G[node]:
            f.write('L\t'+str(mapping[node])+'\t+\t'+str(mapping[to])+"\t+\t0M\n")
    
    if paths:
        for sample in G.graph['samples']:
            path=[]
            for node in nx.topological_sort(G):
                if sample in G.node[node]['sample']:
                    i=mapping[node]
                    path.append(str(i)+"+")
            f.write("P\t"+sample+"\t"+",".join(path)+"\t"+",".join(["0M"]*len(path))+"\n")
    
    f.close()

def write_gml(G,T,outputfile="reference",partition=True,hwm=4000):
    G=G.copy()
    mapping={}
    totn=len(G.graph['samples'])

    for key in G.graph:
        G.graph[key]=str(G.graph[key])
    
    for n,d in G.nodes(data=True):
        mapping[n]=str(n)
        d=G.node[n]
         
        if 'sample' in d:
            G.node[n]['n']=len(d['sample'])
        
        for key in d:
            v=d[key]
            if not isinstance(v,str) and not isinstance(v,int):
                G.node[n][key]=str(v)
        
        if 'seq' not in G.node[n]:
            if isinstance(n,Interval):
                G.node[n]['seq']=T[n.begin:n.end].upper()
            else:
                G.node[n]['seq']=""
        G.node[n]['l']=len(G.node[n]['seq'])
    
    G=nx.relabel_nodes(G,mapping)
    outputfiles=[]
    
    if partition:
        i=0
        for subset in nx.connected_components(G.to_undirected()):
            sgn=[]
            g=G.subgraph(subset)
            gn=G.graph['samples']
            for n in nx.topological_sort(g):
                if sgn==[]:
                    fr=n
                sgn.append(n)
                if G.node[n]['n']==totn: #join/split node
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
    else:
        if not outputfile.endswith(".gml"):
            outputfile=outputfile+'.gml'
        nx.write_gml(G,outputfile)
        outputfiles.append(outputfile)
    
    return outputfiles

def main():
    desc="""
    Type 'reveal <positional_argument> --help' for help on a specific subcommand.\n
    Reveal constructs population reference graphs by aligning multiple whole 
    genomes using recursive exact matching.
    http://www.biorxiv.org/content/early/2015/07/17/022715.
    """
    
    parser = argparse.ArgumentParser(prog="reveal", usage="reveal -h for usage", description=desc)
    parser.add_argument("-l", "--log-level", type=int, dest="loglevel", default=20, help="Log level: 10=debug 20=info (default) 30=warn 40=error 50=fatal.")
    parser.add_argument("--64", dest="sa64", default=False, action="store_true", help="Use 64bit suffix array in the index.")
    
    subparsers = parser.add_subparsers()
    parser_aln = subparsers.add_parser('align',prog="reveal align", description="Construct a population graph from input genomes or other graphs.")
    parser_plot = subparsers.add_parser('plot', prog="reveal plot", description="Generate mumplot for two fasta files.")
    parser_convert = subparsers.add_parser('convert', prog="reveal convert", description="Convert gfa graph to gml.")
    parser_extract = subparsers.add_parser('extract', prog="reveal extract", description="Extract the input sequence from a graph.")
    parser_comp = subparsers.add_parser('comp', prog="reveal comp", description="Reverse complement the graph.")
    parser_finish = subparsers.add_parser('finish', prog="reveal finish", description="Finish a draft assembly by ordering and orienting contigs with respect to a finished reference assembly.")
    parser_matches = subparsers.add_parser('matches', prog="reveal matches", description="Outputs all (multi) m(u/e)ms.")
    parser_subgraph = subparsers.add_parser('subgraph', prog="reveal subgraph", description="Extract subgraph from gfa by specified node ids.")
    parser_bubbles = subparsers.add_parser('bubbles', prog="reveal bubbles", description="Extract all bubbles from the graph.")
    parser_realign = subparsers.add_parser('realign', prog="reveal realign", description="Realign between two nodes in the graph.")
    parser_merge = subparsers.add_parser('merge', prog="reveal merge", description="Combine multiple gfa graphs into a single gfa graph.")
    
    parser_aln.add_argument('inputfiles', nargs='*', help='Fasta or gfa files specifying either assembly/alignment graphs (.gfa) or sequences (.fasta). When only one gfa file is supplied, variants are called within the graph file.')
    parser_aln.add_argument("-o", "--output", dest="output", help="Prefix of the variant and alignment graph files to produce, default is \"sequence1_sequence2\"")
    #parser_aln.add_argument("-p", dest="pcutoff", type=float, default=1e-3, help="If, the probability of observing a MUM of the observed length by random change becomes larger than this cutoff the alignment is stopped (default 1e-3).")
    parser_aln.add_argument("-t", "--threads", dest="threads", type=int, default=0, help = "The number of threads to use for the alignment.")
    parser_aln.add_argument("-m", dest="minlength", type=int, default=15, help="Min length of an exact match (default 20).")
    parser_aln.add_argument("-c", dest="minscore", type=int, default=0, help="Min score of an exact match (default 0), exact maches are scored by their length and penalized by the indel they create with respect to previously accepted exact matches.")
    parser_aln.add_argument("-n", dest="minn", type=int, default=2, help="Only align graph on exact matches that occur in at least this many samples.")
    parser_aln.add_argument("--wp", dest="wpen", type=int, default=1, help="Multiply penalty for a MUM by this number in scoring scheme.")
    parser_aln.add_argument("--ws", dest="wscore", type=int, default=3, help="Multiply length of MUM by this number in scoring scheme.")
    parser_aln.add_argument("--mumplot", dest="mumplot", action="store_true", default=False, help="Save a mumplot for the actual aligned chain of anchors (depends on matplotlib).")
    parser_aln.add_argument("-i", dest="interactive", action="store_true", default=False, help="Show an interactive visualisation of the mumplot (depends on matplotlib).")
    
    parser_aln.add_argument("--sa", dest="sa", default="", help="Specify a preconstructed suffix array to decouple suffix array construction.")
    parser_aln.add_argument("--lcp", dest="lcp", default="", help="Specify a preconstructed lcp array to decouple lcp array construction.")
    parser_aln.add_argument("--cache", dest="cache", default=False, action="store_true", help="When specified, it caches the suffix and lcp array to disk after construction.")
    
    parser_aln.add_argument("-g", dest="minsamples", type=int, default=1, help="Only index nodes that occur in this many samples or more (default 1).")
    parser_aln.add_argument("-x", dest="maxsamples", type=int, default=None, help="Only align nodes that have maximally this many samples (default None).")
    parser_aln.add_argument("-r", dest="reference", type=str, default=None, help="Name of the sequence that should be used as a coordinate system or reference.")
    parser_aln.add_argument("-s", dest="targetsample", type=str, default=None, help="Only align nodes in which this sample occurs.")
    parser_aln.add_argument("--gml", dest="gml", action="store_true", default=False, help="Produce a gml graph instead gfa.")
    parser_aln.add_argument("--gml-max", dest="hwm", default=4000, help="Max number of nodes per graph in gml output.")
    parser_aln.add_argument("--nometa", dest="nometa", action="store_true", default=False, help="Produce a gfa graph without node annotations, to ensure it's parseable by other programs.")
    parser_aln.add_argument("--paths", dest="paths", action="store_true", default=False, help="Output paths in GFA.")
    #parser_aln.add_argument("--align-contigs", dest="contigs", action="store_true", default=False, help="Use when pairwise aligning a set of contigs to a single genome or graph. Contigs are aligned one by one (slow).")
    parser_aln.set_defaults(func=align_cmd)
    
    parser_extract.add_argument('graph', nargs=1, help='gfa file specifying the graph from which the genome should be extracted.')
    parser_extract.add_argument('samples', nargs='*', help='Name of the sample to be extracted from the graph.')
    parser_extract.add_argument("--width", dest="width", type=int, default=100 , help='Line width for fasta output.')
    parser_extract.set_defaults(func=extract_cmd)
    
    parser_plot.add_argument('fastas', nargs=2, help='Two fasta files for which a mumplot should be generated.')
    parser_plot.add_argument("-m", dest="minlength", type=int, default=100, help="Minimum length of exact matches to vizualize (default=100).")
    parser_plot.add_argument("-i", dest="interactive", action="store_true", default=False, help="Wheter to produce interactive plots which allow zooming on the dotplot (default=False).")
    parser_plot.add_argument("-u", dest="uniq", action="store_true", default=False, help="Plot only maximal unique matches.")
    parser_plot.add_argument("-r", dest="region", default=None, help="Highlight interval (as \"<start>:<end>\") with respect to x-axis (first sequence).")
    
    #parser_plot.add_argument("-p", dest="pos", default=None, type=int, help="Position on the reference genome to vizualize.")
    #parser_plot.add_argument("-e", dest="env", default=1000, type=int, help="Size of the region aroung the targeted position to vizualize.")
    parser_plot.set_defaults(func=plot)
    
    parser_comp.add_argument('graph', nargs=1, help='The graph to be reverse complemented.')
    parser_comp.set_defaults(func=comp_cmd)
    
    parser_finish.add_argument('reference', help='Graph or sequence to which query/contigs should be assigned.')
    parser_finish.add_argument('contigs', help='Graph or fasta that is to be reverse complemented with respect to the reference.')
    parser_finish.add_argument("-m", dest="minlength", type=int, default=100, help="Min length of maximal exact matches for considering (default 20).")
    parser_finish.add_argument("--cache", dest="cache", default=False, action="store_true", help="When specified, it caches the text, suffix and lcp array to disk after construction.")
    parser_finish.add_argument("--sa1", dest="sa1", default="", help="Specify a preconstructed suffix array for extracting matches between the two genomes in their current orientation.")
    parser_finish.add_argument("--lcp1", dest="lcp1", default="", help="Specify a preconstructed lcp array for extracting matches between the two genomes in their current orientation.")
    parser_finish.add_argument("--sa2", dest="sa2", default="", help="Specify a preconstructed suffix array for extracting matches between the two genomes in which the sequence of the contigs was reverse complemented.")
    parser_finish.add_argument("--lcp2", dest="lcp2", default="", help="Specify a preconstructed lcp array for extracting matches between the two genomes in which the sequence of the contigs was reverse complemented.")
    parser_finish.set_defaults(func=finish)

    parser_matches.add_argument('reference', help='Graph or sequence to which query/contigs should be assigned.')
    parser_matches.add_argument('contigs', help='Graph or fasta that is to be reverse complemented with respect to the reference.')    
    parser_matches.add_argument("-m", dest="minlength", type=int, default=100, help="Min length of maximal exact matches for considering (default 20).")
    parser_matches.add_argument("--cache", dest="cache", default=False, action="store_true", help="When specified, it caches the text, suffix and lcp array to disk after construction.")
    parser_matches.add_argument("--sa1", dest="sa1", default="", help="Specify a preconstructed suffix array for extracting matches between the two genomes in their current orientation.")
    parser_matches.add_argument("--lcp1", dest="lcp1", default="", help="Specify a preconstructed lcp array for extracting matches between the two genomes in their current orientation.")
    parser_matches.add_argument("--sa2", dest="sa2", default="", help="Specify a preconstructed suffix array for extracting matches between the two genomes in which the sequence of the contigs was reverse complemented.")
    parser_matches.add_argument("--lcp2", dest="lcp2", default="", help="Specify a preconstructed lcp array for extracting matches between the two genomes in which the sequence of the contigs was reverse complemented.")    
    parser_matches.set_defaults(func=matches)
    
    #TOOD: reveal TEXT; write T for any input of graphs and fastas to enable external SA and LCP construction
    
    parser_convert.add_argument('graphs', nargs='*', help='The gfa graph to convert to gml.')
    parser_convert.add_argument("-n", dest="minsamples", type=int, default=1, help="Only align nodes that occcur in this many samples (default 1).")
    parser_convert.add_argument("-x", dest="maxsamples", type=int, default=None, help="Only align nodes that have maximally this many samples (default None).")
    parser_convert.add_argument("-s", dest="targetsample", type=str, default=None, help="Only align nodes in which this sample occurs.")
    parser_convert.add_argument("--gml-max", dest="hwm", default=4000, help="Max number of nodes per graph in gml output.")
    parser_convert.add_argument("--gfa",  action="store_true", dest="gfa", default=False, help="Rewrite gfa file.")
    parser_convert.add_argument("--partition",  action="store_true", dest="partition", default=False, help="Output graph as multiple subgraphs if possible.")
    parser_convert.set_defaults(func=convert)
    
    parser_subgraph.add_argument('inputfiles', nargs='*', help='The gfa graph followed by node ids that make up the subgraph.')
    parser_subgraph.add_argument("-o", dest="outfile", type=str, default="~tmp", help="Prefix of the file to which subgraph will be written.")
    parser_subgraph.add_argument("--gml", dest="gml", action="store_true", default=False, help="Produce a gml graph instead of gfa.")
    parser_subgraph.set_defaults(func=subgraph)
    
    parser_bubbles.add_argument("graph", nargs=1, help='Graph in gfa format from which bubbles are to be extracted.')
    parser_bubbles.add_argument("-r", dest="reference", type=str, default=None, help="Name of the sequence that, if possible, should be used as a coordinate system or reference.")
    parser_bubbles.set_defaults(func=bubbles_cmd)
    
    parser_realign.add_argument("graph", nargs=1, help='Graph in gfa format for which a bubble should be realigned.') 
    parser_realign.add_argument("source", nargs='?', type=int, help='Source node.')
    parser_realign.add_argument("sink", nargs='?', type=int, help='Sink node.')
    parser_realign.add_argument("-m", dest="minlength", type=int, default=20, help="Min length of an exact match (default 20).")
    parser_realign.add_argument("-c", dest="minscore", type=int, default=0, help="Min score of an exact match (default 0), exact maches are scored by their length and penalized by the indel they create with respect to previously accepted exact matches.")
    parser_realign.add_argument("-n", dest="minn", type=int, default=2, help="Only align graph on exact matches that occur in at least this many samples.")
    parser_realign.add_argument("-o", dest="outfile", type=str, default=None, help="File to which realigned graph is to be written.")
    parser_realign.add_argument("--all", action="store_true", dest="all", default=False, help="Trigger realignment for all complex bubbles.")
    parser_realign.add_argument("--maxlen", dest="maxlen", type=int, default=10000000, help="Maximum length of the cumulative sum of all paths that run through the complex bubble.")
    parser_realign.add_argument("--maxsize", dest="maxsize", type=int, default=500, help="Maximum allowed number of nodes that are contained in a complex bubble.")
    parser_realign.set_defaults(func=realign_bubble_cmd)
    
    parser_merge.add_argument("graphs", nargs='*', help='Graphs in gfa format that should be merged.')
    parser_merge.add_argument("-o", dest="outprefix", type=str, default=None, help="Prefix of the file to which merged graph is written.")
    parser_merge.set_defaults(func=merge_cmd)
    
    args = parser.parse_args()
    logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p', level=args.loglevel)
    
    global sa64
    sa64=args.sa64
    
    args.func(args)

def merge_cmd(args):
    if len(args.graphs)<2:
        logging.fatal("Specify multiple gfa files to merge them.")
        return
    
    G=nx.DiGraph()
    for graph in args.graphs:
        logging.info("Adding %s ..." %graph)
        read_gfa(graph,None,"",G)
    
    if args.outprefix!=None:
        write_gfa(G,"",outputfile=args.outprefix+".gfa")
    else:
        write_gfa(G,"",outputfile="_".join([os.path.basename(f)[:os.path.basename(f).rfind('.')] for f in args.graphs])+".gfa")

def bubbles_cmd(args):
    if len(args.graph)<1:
        logging.fatal("Specify a gfa file to extract bubbles.")
        return
    
    reference=args.reference
    G=nx.DiGraph()
    read_gfa(args.graph[0],None,"",G)
    
    if 'samples' in G.graph:
        gori=sorted(G.graph['samples'])
    else:
        gori=[]
    
    sys.stdout.write("#source\tsink\tsubgraph\tref\tpos\tvariant")
    for sample in gori:
        sys.stdout.write("\t%s"%sample)
    sys.stdout.write("\n")
    
    for pair,bubblenodes,size,ordD in bubbles(G):
        gt=G.successors(pair[0])
        gt.sort(key=lambda l: ordD[l]) #topological sort neighbors
        sourcesamples=G.node[pair[0]]['sample']
        sucs=set(gt)
        pres=set(G.predecessors(pair[1]))
        sucs.discard(pair[1])
        pres.discard(pair[0])
        
        simple=True
        for suc in sucs:
            if len(G.successors(suc))!=1:
                simple=False
        for pre in pres:
            if len(G.predecessors(pre))!=1:
                simple=False

        #determine genotypes, associate genotype with numerical value
        if len(gt)<size or not simple:
            #complex bubble; not possible to output the actual alleles, instead we can make a call based on branching at source node and output the nodes that make up the complex bubble so we can export and vizualize the subgraph
            genotypes=['N']*len(gt)
        else:
            genotypes=[None]*len(gt)
            for i,node in enumerate(gt):
                genotypes[i]=G.node[node]['seq']
                if node==pair[1]:
                    genotypes[i]="-"
        
        #make call, associate sample with the numerical value of the genotype
        tmpsource=sourcesamples.copy()
        calls={}
        for i,node in enumerate(gt):
            n=G.node[node]
            for sample in n['offsets'].keys():
                if sample not in tmpsource:
                    continue
                calls[sample]=i
                tmpsource.discard(sample)
        
        d=G.node[pair[0]]
        
        if reference==None or reference not in d['offsets'].keys():
            ref=d['offsets'].keys()[0]
            pos=d['offsets'][ref]+len(d['seq'])+1
        else:
            ref=reference
            pos=d['offsets'][reference]+len(d['seq'])+1
        
        sys.stdout.write("%d\t%d\t%s\t%s\t%s\t%s"%(pair[0],pair[1],",".join([str(x) for x in bubblenodes]),ref,pos,",".join(genotypes)))
        
        for sample in gori:
            if sample in calls:
                sys.stdout.write("\t%s"%calls[sample])
            else:
                sys.stdout.write("\t-")
        sys.stdout.write("\n")

def bubbles(G):
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
        if op==None:
            return -1
        elif entrance(G, ordD_[op]):
            return ordD_[op]
        elif previousEntrance[ordD_[op]]==None: #
            return -1
        else:
            return ordD_[previousEntrance[ordD_[op]]]
        return startVertex
    
    ordD,ordD_,sspairs=superbubble(G)
    
    sg=set()
    for pair in sspairs:
        size=(ordD[pair[1]]-ordD[pair[0]])-1
        bubblenodes=ordD_[ordD[pair[0]]:ordD[pair[1]]+1]
        sourcenode=G.node[pair[0]]
        sourcesamples=sourcenode['sample']
        sinknode=G.node[pair[1]]
        sinksamples=sinknode['sample']
        
        #TODO: supbub algorithm detects invalid bubbles; filter them out here
        if sourcesamples!=sinksamples:
            #print "ERROR, skipping invalid bubble"
            continue
        yield pair,bubblenodes,size,ordD

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
        nodes.add(int(node))
    sg=G.subgraph(nodes)
    if args.gml:
        write_gml(sg,"",outputfile=args.outfile)
    else:
        write_gfa(sg,"",outputfile=args.outfile)

def align_seq(s1,s2,minlength=1,minscore=0,minn=2):
    global t,G,reference,o

    t=IntervalTree()
    
    if sa64:
        idx=reveallib64.index()
    else:
        idx=reveallib.index()
    
    G=nx.DiGraph()
    G.graph['samples']=[]
    reference=None
    o=0
    
    idx.addsample("s1")
    intv=idx.addsequence(s1.upper())
    Intv=Interval(intv[0],intv[1])
    t.add(Intv)
    G.add_node(Intv,sample={"s1"},offsets={"s1":0},aligned=0)

    idx.addsample("s2")
    intv=idx.addsequence(s2.upper())
    Intv=Interval(intv[0],intv[1])
    t.add(Intv)
    G.add_node(Intv,sample={"s2"},offsets={"s2":0},aligned=0)
    
    schemes.ts=t
    schemes.minlength=minlength
    schemes.minscore=minscore
    schemes.minn=minn

    idx.construct()
    
    idx.align(None,graphalign)
    
    alignedbases=0
    totbases=min([idx.n-idx.nsep[0],idx.nsep[0]])
    for node,data in G.nodes(data=True):
        if data['aligned']!=0:
            alignedbases+=(node.end-node.begin)
    
    return alignedbases/float(totbases)

def align_cmd(args):
    if len(args.inputfiles)<=1:
        logging.fatal("Specify at least 2 (g)fa files for creating a reference graph.")
        return
    
    G,idx=align_genomes(args)
    
    if args.output==None:
        args.output="_".join([os.path.basename(f)[:os.path.basename(f).rfind('.')] for f in args.inputfiles])
    
    logging.info("Merging nodes...")
    T=idx.T
    
    if len(G.graph['samples'])>2:
        prune_nodes(G,T) #TODO: do this after every alignment step to reduce memory usage of graph
    
    logging.info("Done.")
    
    alignedbases=0
    alignednodes=0
    totnodes=nx.number_of_nodes(G)

    if idx.nsamples>2: #was multi-alignment
        totbases=idx.n-totnodes
        for node,data in G.nodes(data=True):
            if data['aligned']!=0:
                alignedbases+=(node.end-node.begin)*len(data['sample'])
                alignednodes+=1
    else: #assume seq to graph
        totbases=min([(idx.n-1)-(idx.nsep[0]+1),idx.nsep[0]])
        
        if args.mumplot:
            from matplotlib import pyplot as plt
            if len(G.graph['samples'])>2:
                logging.info("Can't make a mumplot for more than 2 genomes.")
                args.mumplot=False
            else:
                s1=G.graph['samples'][0]
                plt.xlabel(s1)
                s2=G.graph['samples'][1]
                plt.ylabel(s2)
                plt.title("REVEAL "+" ".join(sys.argv[1:]))
        
        for node,data in G.nodes(data=True):
            if data['aligned']!=0:
                l=node.end-node.begin
                alignedbases+=l
                alignednodes+=1
                if args.mumplot:
                    plt.plot([data['offsets'][s1], data['offsets'][s1]+l], [data['offsets'][s2], data['offsets'][s2]+l], 'r-')

    logging.info("%s (%.2f%% identity, %d bases out of %d aligned, %d nodes out of %d aligned)."%("-".join([os.path.basename(f) for f in args.inputfiles]), (alignedbases/float(totbases))*100,alignedbases,totbases,alignednodes,totnodes))
    logging.info("Writing graph...")
    if args.gml:
        graph=write_gml(G,T, hwm=args.hwm, outputfile=args.output)
    else:
        write_gfa(G,T,nometa=args.nometa, outputfile=args.output+'.gfa', paths=args.paths)
        graph=args.output+'.gfa'
    logging.info("Done.")
    logging.info("Alignment graph written to: %s"%graph)
    
    if args.mumplot and idx.nsamples==2:
        plt.plot(0,0,'bx')
        plt.plot(idx.nsep[0],idx.n-idx.nsep[0],'bx')
        logging.info("Storing mumplot as %s.png..."%args.output)
        plt.savefig(args.output+".png")
        logging.info("Done.")
        if args.interactive:
            logging.info("Showing interactive mumplot...")
            plt.show()
            logging.info("Done.")

def align_genomes(args):
    logging.info("Loading input...")
    #global variables to simplify callbacks from c extension
    global t,G,reference,o
    
    reference=args.reference
    
    t=IntervalTree()

    if sa64:
        idx=reveallib64.index(sa=args.sa, lcp=args.lcp, cache=args.cache)
    else:
        idx=reveallib.index(sa=args.sa, lcp=args.lcp, cache=args.cache)
    
    G=nx.DiGraph()
    G.graph['samples']=[]
    o=0
    #schemes.pcutoff=args.pcutoff
    schemes.minlength=args.minlength
    schemes.minscore=args.minscore
    schemes.minn=args.minn
    graph=False
    
    for i,sample in enumerate(args.inputfiles):
        idx.addsample(os.path.basename(sample))
        if sample.endswith(".gfa"):
            graph=True
            #TODO: now applies to all graphs! probably want to have this graph specific if at all...
            logging.info("Reading graph: %s ..." % sample)
            if i==0:
                read_gfa(sample,idx,t,G,minsamples=args.minsamples,
                                        maxsamples=args.maxsamples,
                                        targetsample=args.targetsample)
            else:
                read_gfa(sample,idx,t,G)
            
            if len(G.graph['samples'])==0: #if not from reveal, might not have a header
                G.graph['samples'].append(os.path.basename(sample))

            logging.info("Done.")
        else: #consider it to be a fasta file
            logging.info("Reading fasta: %s ..." % sample)
            G.graph['samples'].append(os.path.basename(sample))
            for name,seq in fasta_reader(sample):
                intv=idx.addsequence(seq.upper())
                logging.debug("Adding interval: %s"%str(intv))
                Intv=Interval(intv[0],intv[1])
                t.add(Intv)
                G.add_node(Intv,sample={os.path.basename(sample)},offsets={os.path.basename(sample):0},aligned=0)
    
    if not nx.is_directed_acyclic_graph(G):
        logging.error("*** Input is not a DAG! ...")
        return
    
    schemes.ts=t
    schemes.G=G
    
    logging.info("Constructing index...")
    idx.construct()
    
    logging.info("Done.")
    
    if len(args.inputfiles)>2:
        logging.info("Constructing multi-alignment...")
        schemes.wscore=args.wscore
        schemes.wpen=args.wpen
        idx.align(schemes.multimumpicker,graphalign,threads=args.threads)
    else:
        if graph:
            logging.info("Constructing graph-alignment...")
            schemes.wscore=args.wscore
            schemes.wpen=args.wpen
            idx.align(schemes.graphmumpicker,graphalign,threads=args.threads)
        else:
            logging.info("Constructing pairwise-alignment...")
            idx.align(None,graphalign,threads=args.threads,wpen=args.wpen,wscore=args.wscore)
    
    return G,idx

def align_contigs(args):
    global t,G,reference,o
    
    reference=args.reference
    schemes.minlength=args.minlength
    schemes.minscore=args.minscore
    schemes.minn=args.minn

    G=nx.DiGraph()
    G.graph['samples']=[]
    
    ref=args.inputfiles[0]
    contigs=args.inputfiles[1]
    
    t=IntervalTree()
    if sa64:
        idx=reveallib64.index()
    else:
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
            G.add_node(Intv,sample={ref},contig={chromname.replace(";","")},coordcontig=chromname.replace(";",""),offsets={ref:0},aligned=0)
    
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
        G.add_node(Intv,sample={contigs},contig={contigname.replace(";","")},coordcontig=contigname.replace(";",""),offsets={contigs:0},aligned=0)

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

    if sa64:
        idx=reveallib64.index()
    else:
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

def align(aobjs,ref=None,minlength=15,minscore=0,minn=2,threads=0,global_align=False,targetsample=None,maxsamples=None,wpen=1,wscore=3):
    #seq should be a list of objects that can be (multi-) aligned by reveal, following possibilities:
    #   - fasta filename
    #   - gfa filename
    #   - tuple of the form (name,seq)
    
    #global variables to simplify callbacks from c extension
    global t,G,reference,o
    reference=ref
    t=IntervalTree()

    if sa64:
        idx=reveallib64.index()
    else:
        idx=reveallib.index()
    
    G=nx.DiGraph()
    H=G
    G.graph['samples']=[]
    o=0
    schemes.minlength=minlength
    schemes.minscore=minscore
    schemes.minn=minn
    graph=False
    
    for aobj in aobjs:
        if isinstance(aobj,tuple):
            assert(len(aobj)==2)
            name,seq=aobj
            idx.addsample(name)
            intv=idx.addsequence(seq.upper())
            if intv[1]-intv[0]>0:
                Intv=Interval(intv[0],intv[1])
                t.add(Intv)
                G.graph['samples'].append(name)
                G.add_node(Intv,sample={name},offsets={name:0},aligned=0)
        elif isinstance(aobj,str):
            if not os.path.isfile(aobj):
                logging.fatal("Not a file, expecting fasta or gfa file.")
                return
            idx.addsample(os.path.basename(aobj))
            if aobj.endswith(".gfa"):
                read_gfa(aobj,idx,t,G,targetsample=targetsample,maxsamples=maxsamples)
                graph=True
            else: #assume a file in fastaformat
                for name,seq in fasta_reader(sample):
                    intv=idx.addsequence(seq.upper())
                    if intv[1]-intv[0]>0:
                        Intv=Interval(intv[0],intv[1])
                        t.add(Intv)
                        G.graph['samples'].append(os.path.basename(sample))
                        G.add_node(Intv,sample={os.path.basename(sample)},offsets={os.path.basename(sample):0},aligned=0)
    
    #connect all sequence at end and beginning
    if global_align:
        G.add_node('start',aligned=1)
        G.add_node('end',aligned=1)
        for node in G.nodes():
            if not isinstance(node,Interval):
                continue
            if len(G.successors(node))==0:
                G.add_edge(node,'end')
            if len(G.predecessors(node))==0:
                G.add_edge('start',node)
    
    if not nx.is_directed_acyclic_graph(G):
        logging.error("*** Input is not a DAG! Not supported.")
        return
    
    schemes.ts=t
    schemes.G=G
    
    idx.construct()
    
    if len(aobjs)>2:
        idx.align(schemes.multimumpicker,graphalign,threads=threads)
    else:
        if graph:
            idx.align(schemes.graphmumpicker,graphalign,threads=threads,wpen=wpen,wscore=wscore)
        else:
            idx.align(None,graphalign,threads=threads,wpen=wpen,wscore=wscore)
    
    prune_nodes(G,idx.T)

    return G,idx







def matches(args): 
    
    if sa64:
        idx=reveallib64.index(sa=args.sa1, lcp=args.lcp1, cache=args.cache) #enable preconstruction of first SA and LCP array
    else:
        idx=reveallib.index(sa=args.sa1, lcp=args.lcp1, cache=args.cache) #enable preconstruction of first SA and LCP array
    
    G=nx.DiGraph()
    G.graph['samples']=[]
    t=IntervalTree()
    
    reffile=os.path.basename(args.reference)
    ctgfile=os.path.basename(args.contigs)
    
    idx.addsample(reffile)
    if args.reference.endswith(".gfa"):
        read_gfa(args.reference,idx,t,G)
    else:
        G.graph['samples'].append(reffile)
        for name,seq in fasta_reader(args.reference):
            intv=idx.addsequence(seq)
            intv=Interval(intv[0],intv[1],name)
            t.add(intv)
            G.add_node(intv,sample={reffile},offsets={reffile:0})
    
    contig2length=dict()
    idx.addsample(ctgfile)
    if args.contigs.endswith(".gfa"):
        read_gfa(args.contigs,idx,t,G)
    else:
        G.graph['samples'].append(ctgfile)
        for name,seq in fasta_reader(args.contigs):
            contig2length[name]=len(seq)
            intv=idx.addsequence(seq)
            intv=Interval(intv[0],intv[1],name)
            t.add(intv)
            G.add_node(intv,sample={ctgfile},offsets={ctgfile:0})
    
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
        if reffile in G.node[nodes[0]]['sample']:
            for node in nodes:
                assert(reffile in G.node[node]['sample']) #check the graph is valid
                refnode2component[node]=ri
                component2refnode[ri]=node
            ri+=1
            refcomponents.append(nodes)
        else:
            for node in nodes:
                assert(ctgfile in G.node[node]['sample']) #check the graph is valid
                ctgnode2component[node]=ci
                component2ctgnode[ci]=node
            ci+=1
            ctgcomponents.append(nodes)
    
    #for each contig, print the length
    for name in contig2length:
        print "#%s\t%d"%(name,contig2length[name])
    
    idx.construct()
    
    print "##refname\trefstart\tctgname\tctgstart\tlength\tn\tunique\torient"
    for mem in idx.getmems(args.minlength) :
        refstart=mem[2][0]
        ctgstart=mem[2][1]
        rnode=t[refstart].pop() #start position on match to node in graph
        cnode=t[ctgstart].pop()
        print "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (rnode[2], refstart-rnode[0], cnode[2], ctgstart-cnode[0], mem[0], mem[1], mem[3], 0)
    
    logging.info("Indexing reverse complement...\n")
    
    ### index reverse complement
    if sa64:
        idx=reveallib64.index(sa=args.sa2, lcp=args.lcp2) #enable preconstruction of second SA and LCP array
    else:
        idx=reveallib.index(sa=args.sa2, lcp=args.lcp2) #enable preconstruction of second SA and LCP array
    
    rcG=nx.DiGraph()
    t=IntervalTree()
    
    idx.addsample(reffile)
    if args.reference.endswith(".gfa"):
        read_gfa(args.reference,idx,t,rcG)
    else:
        rcG.graph['samples']=set([reffile])
        for name,seq in fasta_reader(args.reference):
            intv=idx.addsequence(seq)
            intv=Interval(intv[0],intv[1],name)
            t.add(intv)
            rcG.add_node(intv,sample={reffile},offsets={reffile:0},aligned=0)
            refseq=seq
    
    idx.addsample(ctgfile)
    if args.contigs.endswith(".gfa"):
        read_gfa(args.contigs,idx,t,rcG,revcomp=True)
    else:
        rcG.graph['samples']=set([ctgfile])
        for name,seq in fasta_reader(args.contigs):
            intv=idx.addsequence(rc(seq))
            intv=Interval(intv[0],intv[1],name)
            t.add(intv)
            rcG.add_node(intv,sample={ctgfile},offsets={ctgfile:0},aligned=0)
    
    idx.construct()
    
    for mem in idx.getmems(args.minlength):
        refstart=mem[2][0]
        ctgstart=mem[2][1]
        rnode=t[refstart].pop() #start position on match to node in graph
        cnode=t[ctgstart].pop()
        l=cnode[1]-cnode[0]
        print "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (rnode[2], refstart-rnode[0], cnode[2], l-((ctgstart-cnode[0])+mem[0]), mem[0], mem[1], mem[3], 1)


def finish(args):
    if sa64:
        idx=reveallib64.index(sa=args.sa1, lcp=args.lcp1, cache=args.cache) #enable preconstruction of first SA and LCP array
    else:
        idx=reveallib.index(sa=args.sa1, lcp=args.lcp1, cache=args.cache) #enable preconstruction of first SA and LCP array
    
    G=nx.DiGraph()
    G.graph['samples']=[]
    t=IntervalTree()
    
    reffile=os.path.basename(args.reference)
    ctgfile=os.path.basename(args.contigs)
    
    idx.addsample(reffile)
    if args.reference.endswith(".gfa"):
        read_gfa(args.reference,idx,t,G)
    else:
        G.graph['samples'].append(reffile)
        for name,seq in fasta_reader(args.reference):
            intv=idx.addsequence(seq)
            intv=Interval(intv[0],intv[1],name)
            t.add(intv)
            G.add_node(intv,sample={reffile},offsets={reffile:0})
    
    contig2length=dict()
    contig2seq=dict()
    
    idx.addsample(ctgfile)
    if args.contigs.endswith(".gfa"):
        read_gfa(args.contigs,idx,t,G)
    else:
        G.graph['samples'].append(ctgfile)
        for name,seq in fasta_reader(args.contigs):
            contig2length[name]=len(seq)
            contig2seq[name]=seq
            intv=idx.addsequence(seq)
            intv=Interval(intv[0],intv[1],name)
            t.add(intv)
            G.add_node(intv,sample={ctgfile},offsets={ctgfile:0})
    
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
        if reffile in G.node[nodes[0]]['sample']:
            for node in nodes:
                assert(reffile in G.node[node]['sample']) #check the graph is valid
                refnode2component[node]=ri
                component2refnode[ri]=node
            ri+=1
            refcomponents.append(nodes)
        else:
            for node in nodes:
                assert(ctgfile in G.node[node]['sample']) #check the graph is valid
                ctgnode2component[node]=ci
                component2ctgnode[ci]=node
            ci+=1
            ctgcomponents.append(nodes)
    
    idx.construct()
    mems=[]

    for mem in idx.getmems(args.minlength) :
        refstart=mem[2][0]
        ctgstart=mem[2][1]
        rnode=t[refstart].pop() #start position on match to node in graph
        cnode=t[ctgstart].pop()
        mems.append((rnode[2], refstart-rnode[0], cnode[2], ctgstart-cnode[0], mem[0], mem[1], mem[3], 0))
    
    logging.info("Indexing reverse complement...\n")
    
    ### index reverse complement
    if sa64:
        idx=reveallib64.index(sa=args.sa2, lcp=args.lcp2) #enable preconstruction of second SA and LCP array
    else:
        idx=reveallib.index(sa=args.sa2, lcp=args.lcp2) #enable preconstruction of second SA and LCP array
    
    rcG=nx.DiGraph()
    t=IntervalTree()
    
    idx.addsample(reffile)
    if args.reference.endswith(".gfa"):
        read_gfa(args.reference,idx,t,rcG)
    else:
        rcG.graph['samples']=set([reffile])
        for name,seq in fasta_reader(args.reference):
            intv=idx.addsequence(seq)
            intv=Interval(intv[0],intv[1],name)
            t.add(intv)
            rcG.add_node(intv,sample={reffile},offsets={reffile:0},aligned=0)
            refseq=seq
    
    idx.addsample(ctgfile)
    if args.contigs.endswith(".gfa"):
        read_gfa(args.contigs,idx,t,rcG,revcomp=True)
    else:
        rcG.graph['samples']=set([ctgfile])
        for name,seq in fasta_reader(args.contigs):
            intv=idx.addsequence(rc(seq))
            intv=Interval(intv[0],intv[1],name)
            t.add(intv)
            rcG.add_node(intv,sample={ctgfile},offsets={ctgfile:0},aligned=0)
    
    idx.construct()
    
    for mem in idx.getmems(args.minlength):
        refstart=mem[2][0]
        ctgstart=mem[2][1]
        rnode=t[refstart].pop() #start position on match to node in graph
        cnode=t[ctgstart].pop()
        l=cnode[1]-cnode[0]
        mems.append((rnode[2], refstart-rnode[0], cnode[2], l-((ctgstart-cnode[0])+mem[0]), mem[0], mem[1], mem[3], 1))
    
    #determine best scoring chain per contig
    ctg2mums=dict()
    
    logging.info("Relating exact matches to contigs...")
    for mem in mems:
        refchrom, refstart, ctg, ctgstart, l, n, u, o = mem
        refstart=int(refstart)
        ctgstart=int(ctgstart)
        l=int(l)
        n=int(n)
        u=int(u)
        o=int(o)
        if ctg in ctg2mums:
            if refchrom in ctg2mums[ctg]:
                ctg2mums[ctg][refchrom].append((refstart,ctgstart,l,n,u,o))
            else:
                ctg2mums[ctg][refchrom]=[(refstart,ctgstart,l,n,u,o)]
        else:
            ctg2mums[ctg]=dict({refchrom : [(refstart,ctgstart,l,n,u,o)]})

    logging.info("Done.")
    
    ref2ctg=dict()
    
    logging.info("Mapping contigs to reference.")    
    #for each contig determine the best location on the reference
    for ctg in ctg2mums:
        bestp=[]
        bestscore=0
        bestref=None
        for ref in ctg2mums[ctg]:
            mems=ctg2mums[ctg][ref]
            #determine optimal chain of mems
            path,score=bestpath(mems)
            if score>bestscore:
                bestscore=score
                bestp=path
                bestref=ref
                revcomp=False
                bestmems=mems
            
            #determine optimal chain of mems in reverse complement
            rpath,rscore=bestpath(mems,rc=True)
            
            if rscore>bestscore:
                bestscore=rscore
                bestp=rpath
                bestref=ref
                revcomp=True
                bestmems=mems
            
        if bestp==[]:
            continue
        
        begin=bestp[-1][0] #bestp[-1][0][0]
        end=bestp[0][0]+bestp[0][2] #bestp[0][0][0]
        alength=end-begin

        score=bestscore
        chaincov=sum([m[2] for m in bestp])
        chainlength=len(bestp)
         
        if bestref in ref2ctg:
            ref2ctg[bestref].append((ctg,revcomp,bestp,begin,end,contig2length[ctg]))
        else:
            ref2ctg[bestref]=[(ctg,revcomp,bestp,begin,end,contig2length[ctg])]
    
    
    resbase=reffile[:reffile.rfind('.')]+"_"+ctgfile[:ctgfile.rfind('.')]
    with open(resbase+".fasta",'w') as finished:
        #for each reference chromosome, order the assigned contigs
        for ref in ref2ctg:
            logging.info("Determining contig order for: %s"%ref)
            ref2ctg[ref].sort(key= lambda c: c[3]) #sort by start position of chains
            ctgs=ref2ctg[ref]
            frp=0
            plt.clf()
            plt.figure(0,figsize=(20,20))
            ax = plt.axes()
            
            plt.title(ref)
            offset=0
            yticks=[]
            yticklabels=[]
            
            finished.write(">%s_%s\n"%(ref,ctg))

            for ctg,revcomp,path,begin,end,ctglength in ctgs:

                if end<frp:
                    logging.info("skipping contained contig %s with length %d"%(ctg,ctglength))
                    continue
                else:
                    frp=end
                
                for mem in ctg2mums[ctg][ref]:
                    if revcomp:
                        ax.plot([mem[0],mem[0]+mem[2]],[offset+mem[1]+mem[2],offset+mem[1]],'g-')
                    else:
                        ax.plot([mem[0],mem[0]+mem[2]],[offset+mem[1],offset+mem[1]+mem[2]],'r-')
                
                if revcomp:
                    for mem in path:
                        ax.plot([mem[0],mem[0]+mem[2]],[offset+mem[1]+mem[2],offset+mem[1]],'g-',linewidth=2)
                else:
                    for mem in path:
                        ax.plot([mem[0],mem[0]+mem[2]],[offset+mem[1],offset+mem[1]+mem[2]],'r-',linewidth=2)
                ax.axhline(offset+ctglength)
                offset=offset+ctglength
                yticks.append(offset)
                yticklabels.append(ctg[0:20])
                
                if revcomp:
                    finished.write(rc(contig2seq[ctg]))
                else:
                    finished.write(contig2seq[ctg])
                finished.write("N")
            
            finished.write("\n")
            ax.set_yticks(yticks)
            ax.set_yticklabels(yticklabels)
            plt.savefig(resbase+"_"+ref.split()[0]+".png")

def bestpath(mems,rc=False):
    mems.sort(key=lambda mem: mem[0]) #sort by reference position
    init=(None, None, 0, 0, 0, 0)
    link=dict()
    score=dict({init:0})
    active=[init]
    processed=[]
    start=init
    
    #(889051, 1124491, 1113, 2, 1, 1)
    
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
        
        best=None
        w=None

        for amem in active:
            #calculate score of connecting to active point
            if rc:
                if w==None:
                    w=score[amem]+mem[2]
                    best=amem
                elif amem[1] >= mem[1]+mem[2]:
                    penalty=abs( (amem[1]-(mem[1]+mem[2]) ) - (mem[0]-(amem[0]+amem[2])  ) )
                    tmpw=score[amem]+mem[2]-penalty
                    if tmpw>w:
                        w=tmpw
                        best=amem
            else:
                if w==None:
                    w=score[amem]+mem[2]
                    best=amem
                elif amem[1] <= mem[1]:
                    penalty=abs( ((amem[1]+amem[2])-mem[1]) - ((amem[0]+amem[2])-mem[0]) )
                    tmpw=score[amem]+mem[2]-penalty
                    if tmpw>w:
                        w=tmpw
                        best=amem
        
        assert(best!=None)
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
    
    #if start!=path[-1]:
    #    print "subset of mums scores higher path!",path
    
    if score[start]<0:
        print "Negative score at start",start,score[start]
        maxscore-=score[start]
    
    return path,maxscore


def plot(args):
    #from matplotlib import pyplot as plt
    vertgaps=[]
    horzgaps=[]
    if len(args.fastas)==2:
        #get mmems for forward orientation
        if sa64:
            idx=reveallib64.index()
        else:
            idx=reveallib.index()

        for sample in args.fastas:
            idx.addsample(sample)
            for name,seq in fasta_reader(sample,truncN=False):
                intv=idx.addsequence(seq.upper())
                break #expect only one sequence per fasta for now
        idx.construct()
        
        if args.uniq:
            print "Extracting mums..."
            mmems=[(mem[0],mem[1],mem[2],0) for mem in idx.getmums(args.minlength)]
        else:
            print "Extracting mems..."
            mmems=[(mem[0],mem[1],mem[2],0,mem[3]) for mem in idx.getmems(args.minlength)]
        print "done."
        
        sep=idx.nsep[0]
        
        #get mmems for reverse orientation
        if sa64:
            idx=reveallib64.index()
        else:
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
            pc=None
            for i,c in enumerate(seq):
                if c=='N' and (pc!='N' and pc!=None):
                    vertgaps.append(i)
                elif pc=='N' and c!='N':
                    vertgaps.append(i)
                pc=c
            #vertgaps=[i for i,c in enumerate(seq) if c=='N']
            ls+=len(seq)
            intv=idx.addsequence(rc(seq.upper()))
            break #expect only one sequence per fasta for now
        idx.construct()

        if args.uniq:
            print "Extracting RC mums..."            
            mmems+=[(mem[0],mem[1],mem[2],1) for mem in idx.getmums(args.minlength)]
        else:
            print "Extracting RC mems..."            
            mmems+=[(mem[0],mem[1],mem[2],1,mem[3]) for mem in idx.getmems(args.minlength)]
        print "done."
    
    elif len(args.fastas)==1:

        if sa64:
            idx=reveallib64.index()
        else:
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

        if args.uniq:
            mmems=[(mem[0],mem[1],mem[2],0) for mem in idx.getmums(args.minlength) if mem[0]>args.minlength]
        else:
            mmems=[(mem[0],mem[1],mem[2],0,mem[3]) for mem in idx.getmems(args.minlength) if mem[0]>args.minlength]

    else:
        logging.fatal("Can only create mumplot for 2 sequences or self plot for 1 sequence.")
        return
    
    print "Drawing",len(mmems),"matches."
    
    #pos=args.pos
    #dist=args.env
    
    if args.region!=None:
        start,end=args.region.split(":")
        start=int(start)
        end=int(end)
    else:
        start=0
        end=idx.nsep[0]
    
    for mem in mmems:
        sps=sorted(mem[2])
        l=mem[0]
        sp1=sps[0]
        sp2=sps[1]-sep
        ep1=sp1+l
        ep2=sp2+l
        
        if sp1>start and ep1<end:
            if mem[3]==0:
                if args.uniq:
                    plt.plot([sp1,ep1],[sp2,ep2],'r-')
                else:
                    if mem[4]==0: #non-uniq
                        plt.plot([sp1,ep1],[sp2,ep2],'y-')
                    else:
                        plt.plot([sp1,ep1],[sp2,ep2],'r-')
            else:
                if args.uniq: #only uniq matches in the list
                    plt.plot([sp1,ep1],[ls-sp2,ls-ep2],'g-')
                else:
                    if mem[4]==0: #non-uniq
                        plt.plot([sp1,ep1],[ls-sp2,ls-ep2],'y-')
                    else:
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
    
    if args.region!=None:
        start,end=args.region.split(":")
        plt.axvline(x=int(start),linewidth=3,color='b',linestyle='dashed')
        plt.axvline(x=int(end),linewidth=3,color='b',linestyle='dashed')
    #plt.plot([1040103,1065571],[170000,170000],'g-',linewidth=2.0)
    
    if args.interactive:
        plt.show()
    else:
        fn1=args.fastas[0][0:args.fastas[0].rfind('.')] if args.fastas[0].find('.')!=-1 else args.fastas[0]
        fn2=args.fastas[1][0:args.fastas[1].rfind('.')] if args.fastas[1].find('.')!=-1 else args.fastas[1]
        plt.savefig(fn1+"_"+fn2+".png")


def realign_bubble_cmd(args):
    if len(args.graph)<1:
        logging.fatal("Specify a gfa file for which bubbles should be realigned.")
        return
    
    G=nx.DiGraph()
    read_gfa(args.graph[0],None,"",G)
    
    if args.all:
        G=realign_all(G,minscore=args.minscore,minlength=args.minlength,minn=args.minn,maxsize=args.maxsize,maxlen=args.maxlen)
    else:
        G=realign_bubble(G,args.source,args.sink,minscore=args.minscore,minlength=args.minlength,minn=args.minn,maxsize=args.maxsize,maxlen=args.maxlen)
    
    if args.outfile==None:
        fn=args.graph[0].replace(".gfa",".realigned.gfa")
    else:
        fn=args.outfile
    
    write_gfa(G,"",outputfile=fn)

def realign_bubble(G,source,sink,minscore=0,minlength=20,minn=2,maxsize=100,maxlen=10000000):
    #print "Realigning graph between %s and %s"%(source,sink)
    nn=max(G.nodes())+1
    bubblenodes=[]

    assert(source in G)
    assert(sink in G)
    sourcesamples=G.node[source]['sample']
    sinksamples=G.node[sink]['sample']
    assert(sourcesamples==sinksamples)
    
    bubblesnodes=[]
    add=False
    for node in nx.topological_sort(G):
        if node==source:
            add=True
        if add:
            bubblenodes.append(node)
        if node==sink:
            add=False
    
    sg=G.subgraph(bubblenodes)
    d={}
    aobjs=[]
    cumsum=0
    
    #extract all paths
    for sample in sourcesamples.intersection(sinksamples): #should be equal..
        #TODO: can also be that input was an assembly graph! What to do...
        for i,seq in enumerate(extract(sg,sample)): #has to be just one component
            pass
        assert(i==0)
        cumsum+=len(seq)
        if len(seq)>0:
            aobjs.append((sample,seq))
        
        if cumsum>maxlen:
            print "Bubble (%s,%s) is too big. Increase --maxlen."%(source,sink)
            return G
        #TODO: check that bubbles dont that get too big for multi alignment, put some upper limit!
    
    ng,idx=align(aobjs,global_align=False,minscore=minscore,minlength=minlength,minn=minn)
    T=idx.T
    prune_nodes(ng,T)
    
    seq2node(ng,T)
    mapping={}
    startnodes=set()
    endnodes=set()
    
    for node,data in ng.nodes(data=True): 
        if len(ng.predecessors(node))==0:
            startnodes.add(nn)
        if len(ng.successors(node))==0:
            endnodes.add(nn)
        corrected=dict()
        for sample in data['offsets']:
            corrected[sample]=data['offsets'][sample]+G.node[source]['offsets'][sample]
        ng.node[node]['offsets']=corrected
        mapping[node]=nn
        nn+=1
    
    ng=nx.relabel_nodes(ng,mapping) #relabel nodes to get rid of Intervals
    
    #store edges that have to be reconnected after removing
    pre=G.predecessors(source)
    suc=G.successors(sink)
    
    for node in bubblenodes: #remove all bubblenodes from the original graph except for source and sink
        G.remove_node(node)
    for node,data in ng.nodes(data=True): #add nodes from newly aligned graph to original graph, except for source and sink, update those with seq attribute
        G.add_node(node,data)
    for edge in ng.edges():
        G.add_edge(*edge)
    
    #reconnect nodes
    for p in pre:
        for start in startnodes:
            G.add_edge(p,start)
    for s in suc:
        for end in endnodes:
            G.add_edge(end,s)
    
    return G

def realign_all(G,minscore=0,minlength=20,minn=2,maxlen=10000000,maxsize=100):
    complexbubbles=dict()
    source2sink=dict()
    sink2source=dict()
    
    #detect all complex bubbles
    for pair,bubblenodes,size,ordD in bubbles(G):
        sourcesamples=G.node[pair[0]]['sample']
        sinksamples=G.node[pair[1]]['sample']

        if not isinstance(sourcesamples,set):
            print pair,"skipping, not set"
            continue
        if sourcesamples!=sinksamples:
            print pair,sourcesample,sinksamples
            continue
        
        sucs=set(G.successors(pair[0]))
        pres=set(G.predecessors(pair[1]))
        sucs.discard(pair[1])
        pres.discard(pair[0])
        
        simple=True
        if len(sucs)!=len(pres) or len(sucs)<size or len(pres)<size:
            simple=False
        else:
            for suc in sucs:
                if len(G.successors(suc))!=1:
                    simple=False
            for pre in pres:
                if len(G.predecessors(pre))!=1:
                    simple=False

        if not simple and size<maxsize:
            assert(len(bubblenodes)<maxsize)
            complexbubbles[pair]=bubblenodes
            source2sink[pair[0]]=pair[1]
            sink2source[pair[1]]=pair[0]
     
    #join complex bubbles that share a source/sink nodes
    converged=False
    while not converged:
        converged=True
        for source,sink in complexbubbles:
            if sink in source2sink:
                #sink is also a source, combine the two
                #print "sink is also source, merging",(source,sink),"with",(sink,source2sink[sink])
                complexbubbles[(source,source2sink[sink])]=complexbubbles[(source,sink)]+complexbubbles[(sink,source2sink[sink])]
                del complexbubbles[(sink, source2sink[sink])]
                del complexbubbles[(source,sink)]
                source2sink[source]=source2sink[sink]
                sink2source[source2sink[sink]]=source
                del sink2source[sink]
                del source2sink[sink]
                converged=False
                break
            elif source in sink2source:
                #source is also sink, combine
                #print "source is also sink, merging",(sink2source[source],source),"with",(source,sink)
                complexbubbles[(sink2source[source],sink)]=complexbubbles[(source,sink)]+complexbubbles[(sink2source[source],source)]
                del complexbubbles[(sink2source[source],source)]
                del complexbubbles[(source,sink)]
                source2sink[sink2source[source]]=sink
                sink2source[sink]=sink2source[source]
                del sink2source[source]
                del source2sink[source]
                converged=False
                break
    
    #filter out any bubbles that are a subset of another complex bubble
    distinctbubbles=dict()
    values=complexbubbles.values()
    for source,sink in complexbubbles:
        nodes=complexbubbles[(source,sink)]
        for bubble in values:
            if set(bubble).issuperset(set(nodes)) and not set(bubble)==set(nodes):
                #contained, dont add
                #print source,sink,"is contained"
                break
        else:
            distinctbubbles[(source,sink)]=nodes
    
    for bubble in distinctbubbles:
        assert(len(distinctbubbles[bubble])<500)

    #print "Total of %d bubbles to realign"%len(distinctbubbles)
    i=1
    for source,sink in distinctbubbles:
        #print i,"realigning",source,sink,len(distinctbubbles[(source,sink)])
        G=realign_bubble(G,source,sink,minscore=minscore,minlength=minlength,minn=minn)
        i+=1
    
    return G

def seq2node(G,T,toupper=True):
    for node in G:
        if isinstance(node,Interval):
            if toupper:
                G.node[node]['seq']=T[node.begin:node.end].upper()
            else:
                G.node[node]['seq']=T[node.begin:node.end]

def comp(G):
    for node in G.node:
        G.node[node]['seq']=rc(G.node[node]['seq'])
    
    genome2length=dict()
    #relabel the offsets, determine the length of all genomes in the graph, then l-pos
    for sample in G.graph['samples']:
        maxp=0
        for node,data in G.nodes(data=True):
            if sample in data['offsets']:
                if data['offsets'][sample]+len(data['seq'])>maxp:
                    maxp=data['offsets'][sample]+len(data['seq'])
        genome2length[sample]=maxp
    
    for sample in G.graph['samples']:
        for node,data in G.nodes(data=True):
            if sample in data['offsets']:
                G.node[node]['offsets'][sample]=genome2length[sample]-(G.node[node]['offsets'][sample]+len(data['seq']))
    
    G.reverse(copy=False)
    return G

def comp_cmd(args):
    g=nx.DiGraph()
    g.graph['samples']=[]
    read_gfa(args.graph[0],None,None,g,targetsample=None)
    g=comp(g)
    write_gfa(g,"",outputfile=args.graph[0].replace('.gfa','.rc.gfa'), nometa=False)

def convert(args):
    for graph in args.graphs:
        g=nx.DiGraph()
        g.graph['samples']=[]
        if graph.endswith(".gfa"): #gfa to gml/gfa
            read_gfa(graph,None,None,g,minsamples=args.minsamples,
                                 maxsamples=args.maxsamples,
                                 targetsample=args.targetsample)
            if args.gfa:
                fn=graph.replace(".gfa",".rewrite.gfa")
                graph=write_gfa(g,"", outputfile=fn)
                logging.info("gfa graph written to: %s"%fn)
            else:
                fn=write_gml(g,"", hwm=args.hwm, outputfile=graph.replace(".gfa",""), partition=args.partition)
                logging.info("gml graph written to: %s"%fn)
        elif graph.endswith(".fa") or graph.endswith(".fasta") or graph.endswith(".fna"): #assume fasta to gfa
            i=0
            for name,seq in fasta_reader(graph):
                g.graph['samples'].append(os.path.basename(graph))
                g.add_node(i,sample={os.path.basename(graph)},offsets={os.path.basename(graph):0},seq=seq)
                i+=1
            filename=graph[:graph.rfind(".")]+".gfa"
            write_gfa(g,"", outputfile=filename)
            logging.info("gfa graph written to: %s"%filename)
        else:
            logging.fatal("Unknown filetype, need gfa or fasta extension.")
            return

def extract_cmd(args):
    if not args.graph[0].endswith(".gfa"):
        logging.fatal("Invalid gfa file.")
        return
    width=args.width
    G=nx.DiGraph()
    read_gfa(args.graph[0], None, None, G)
    try:
        i=0
        for sample in args.samples:
            for seq in extract(G,sample):
                sys.stdout.write(">"+sample+" "+str(i)+"\n")
                f=0
                for i in xrange(width,len(seq),width):
                    sys.stdout.write(seq[f:i]+'\n')
                    f=i
                sys.stdout.write(seq[f:]+'\n')
                i+=1
    except IOError:
        try:
            sys.stdout.close()
        except IOError:
            pass
        try:
            sys.stderr.close()
        except IOError:
            pass

def extract(G,sample):
    sg=[]
    for node,data in G.nodes(data=True):
        if sample in data['sample']:
            sg.append(node)
    for i,contig in enumerate(nx.connected_components(G.subgraph(sg).to_undirected())):
        seq=""
        nodecount=0
        for node in nx.topological_sort(G.subgraph(contig)):
            seq+=G.node[node]['seq']
        yield seq
