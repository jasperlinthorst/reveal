#!/usr/bin/env python

from intervaltree import Interval, IntervalTree
import networkx as nx
from collections import defaultdict, deque
import threading
import reveallib
from matplotlib import pyplot as plt
import argparse
import logging

def printSA(index,maxline=100,start=0,end=100):
    sa=index.SA
    lcp=index.LCP
    t=idx.T
    sai=idx.SAi
    print len(sa), len(lcp)
    assert(len(sa)==len(lcp))
    for s,l in zip(sa[start:end],lcp[start:end]):
        print str(s).zfill(8), str(l).zfill(6), t[s:s+l].ljust(maxline) if lcp<=maxline else t[s:s+maxline].ljust(maxline), so[s]

def assertSA(idx):
    print "T"
    t=idx.T
    print "SA"
    sa=idx.SA
    print "LCP"
    lcp=idx.LCP
    for i,(s,l) in enumerate(zip(sa,lcp)):
        if i==0:
            continue
        t1=t[s:s+l]
        t2=t[sa[i-1]:sa[i-1]+l]
        if t1!=t2:
            printSA(idx,start=i-3,end=i+2)
            print i,l,t[s:s+l]
            print i-1,l,t[sa[i-1]:sa[i-1]+l]
        if s+l!=len(t)-1 and sa[i-1]+l!=len(t)-1:
            if not(t[s+l]!=t[sa[i-1]+l] or t[s+l].islower() or t[s+l]=='N'):
                print "printing SA for ",i
                if i>2:
                    printSA(idx,start=i-3,end=i+2,maxline=150)
                else:
                    printSA(idx,start=0,end=i+2,maxline=150)

                print "Done"
                assert(t[s+l]!=t[sa[i-1]+l] or t[s+l].islower() or t[s+l]=='N')
        assert(t[s:s+l]==t[sa[i-1]:sa[i-1]+l])

def fasta_reader(fn):
    seq=""
    with open(fn,'r') as ff:
        for line in ff:
            if line.startswith(">"):
                name=line.rstrip().replace(">","")
                if seq!="":
                    yield name,seq
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
            if G.has_edge(refnode,e[0]):
                print mns, refnode
                write_gml(G,T,outputfile="test.gml")
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

def graphalign(l,index,n,sp):
    nodes=index.nodes
    if l<minlength or len(nodes)==0:
        return
    #assert(n==len(sp))
    mns=[]
    topop=[]
    T=index.T
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

def mumpicker(multimums,idx):
    maxscore=0
    maxl=0
    nomum=(0,idx,0,[])
    bestmum=nomum
    for multimum in multimums:
        l,n,sp=multimum
        if n==idx.nsamples:
            if l>maxl:
                maxl=l
                bestmum=(l,idx,n,sp)
    p=(1-((1-(0.25**bestmum[0]))**(((idx.n/float(idx.nsamples))**2))))**(idx.nsamples-1)
    if p>pcutoff:
        bestmum=(0,idx,0,[])
        multimums.sort(key=lambda l: l[0]*l[1],reverse=True)
        for multimum in multimums:
            l,n,sp=multimum
            p=(1-((1-(0.25**l))**(((idx.n/idx.nsamples)**2))))**(idx.nsamples-1)
            if (p<=pcutoff):
                return (l,idx,n,sp)
        #not a single significant multi-mum
        return nomum
    else:
        #one multimum for all samples index
        return bestmum

def mumpicker2(multimums,idx):
    bestmum_by_n={}
    bestmum=(0,idx,0,[])
    for multimum in multimums:
        l,n,sp=multimum
        if n in bestmum_by_n:
            if l>bestmum_by_n[n][0]:
                bestmum_by_n[n]=(l,idx,n,sp)
        else:
            bestmum_by_n[n]=(l,idx,n,sp)
    for key in sorted(bestmum_by_n,reverse=True):
        if bestmum_by_n[key][0]>=minlength:
            bestmum=bestmum_by_n[key]
            break
    return bestmum

def prune_nodes(G):
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

def read_gfa(gfafile, index, tree, graph, minsamples=2, targetsample=None):
    f=open(gfafile,'r')
    sep=";"
    mapping={} #temp mapping object
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
                if len(ann['ORI'])<minsamples or (targetsample!=None and targetsample in ann['ORI']): #dont index these nodes, just add them to the graph
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
    i=1
    mapping={}
    for node,data in G.nodes_iter(data=True):
        f.write('S\t'+str(i)+'\t'+T[node.begin:node.end].upper())
        if not(vg):
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
        mapping[node]=i
        i+=1
    for node1,node2 in G.edges_iter():
        f.write('L\t'+str(mapping[node1])+'\t+\t'+str(mapping[node2])+"\t+\t0M\n")
    f.close()

def write_gml(G,T,outputfile="reference.gml"):
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
    nx.write_gml(G,outputfile)

def pan_core_dist(G):
    d=[0]*idx.nsamples
    for node,data in G.nodes_iter(data=True):
        if len(G.in_edges(node)+G.out_edges(node))<2:
            #print "skipping", node
            continue
        d[len(data['sample'])-1]+=node.end-node.begin
    return d

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
    parser.add_argument("-r", dest="reference", type=str, default=None, help="Name of the sequence that should be used as a coordinate system or reference.")
    parser.add_argument("-s", dest="targetsample", type=str, default=None, help="Only align nodes in which this sample occurs.")
    parser.add_argument("--gml", dest="gml", action="store_true", default=False, help="Produce a gml graph instead gfa.")
    parser.add_argument("--vg", dest="vg", action="store_true", default=False, help="Produce a gfa graph without node annotations, to ensure it's parseable by vg.")
    args = parser.parse_args()
    logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p', level=args.loglevel)
    
    if len(args.inputfiles)==0:
        logging.fatal("Specify at least 2 (g)fa files for creating a reference graph.")
        return
    
    if args.output==None:
        args.output="_".join([f[:f.rfind('.')] for f in args.inputfiles])
    
    #globale variables to simplify callbacks from c extension
    global idx,t,G,T,so,pcutoff,minlength,reference
    
    reference=args.reference
    
    t=IntervalTree()
    idx=reveallib.index()
    G=nx.DiGraph()
    G.graph['samples']=[]
    pcutoff=args.pcutoff
    minlength=args.minlength
    
    for sample in args.inputfiles:
        idx.addsample(sample)
        if sample.endswith(".gfa"):
            read_gfa(sample,idx,t,G,minsamples=args.minsamples,
                                    targetsample=args.targetsample)
        else: #consider it to be a fasta file
            G.graph['samples'].append(sample)
            for name,seq in fasta_reader(sample):
                intv=idx.addsequence(seq.upper())
                i=Interval(intv[0],intv[1])
                t.add(i)
                G.add_node(i,sample={sample},contig={name},coordsample=sample,coordcontig=name,start=0,aligned=0)
    
    for sep in idx.nsep:
        assert(idx.T[sep]=='$')
    
    idx.construct()
    sa=idx.SA
    sai=idx.SAi
    lcp=idx.LCP

    T=idx.T
    so=idx.SO
    
    if len(args.inputfiles)>2:
        idx.align(mumpicker2,graphalign,threads=args.threads)
    else:
        idx.align(None,graphalign,threads=args.threads)
    
    T=idx.T
    prune_nodes(G)
    
    if args.gml:
        write_gml(G,T, outputfile=args.output+'.gml')
    else:
        write_gfa(G,T,vg=args.vg, outputfile=args.output+'.gfa')