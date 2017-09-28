import logging
from utils import *
from intervaltree import Interval, IntervalTree
from collections import defaultdict, deque
import os

import reveallib
import reveallib64
import schemes

try:
    from matplotlib import pyplot as plt
    from matplotlib import patches as patches
except:
    pass

def breaknode(node,pos,l):
    logging.debug("Breaking node: %s"%str(node))
    att=G.node[node]
    in_edges=G.in_edges(node,data=True)
    out_edges=G.out_edges(node,data=True)
    mn=Interval(pos,pos+l)
    other=set()
    
    if mn==node: #no breaking needed
        t.remove(node)
        return node,other
    
    paths=set()
    
    moffsets=dict()
    for s in att['offsets']:
        moffsets[s]=att['offsets'][s]+(pos-node.begin)
        paths.add(s)
    
    soffsets=dict()
    for s in att['offsets']:
        soffsets[s]=att['offsets'][s]+((pos+l)-node.begin)
    
    G.add_node(mn,offsets=moffsets,aligned=0)#create merge node
    
    if (node[0]!=pos):
        pn=Interval(node[0],pos)
        G.add_node(pn,offsets=att['offsets'],aligned=0)#create prefix node
        assert(not G.has_edge(pn,mn))
        assert(not G.has_edge(mn,pn))
        G.add_edge(pn,mn,paths=paths)
        t.add(pn)
        other.add(pn)
    else:
        pn=mn
    if (node[1]!=pos+l):
        sn=Interval(pos+l,node[1])
        G.add_node(sn,offsets=soffsets,aligned=0)#create suffix node
        assert(not G.has_edge(mn,sn))
        assert(not G.has_edge(sn,mn))
        G.add_edge(mn,sn,paths=paths)
        t.add(sn)
        other.add(sn)
    else:
        sn=mn
    G.remove_node(node)                     #update Graph
    t.remove(node)                          #update intervaltree
    for fro,to,d in in_edges:
        assert(not G.has_edge(fro,pn))
        assert(not G.has_edge(pn,fro))
        G.add_edge(fro,pn,**d)
    for fro,to,d in out_edges:
        assert(not G.has_edge(sn,to))
        assert(not G.has_edge(to,sn))
        G.add_edge(sn,to,**d)
    
    logging.debug("Leading/Trailing node(s): %s"%str(other))
    logging.debug("Matching node: %s"%str(mn))

    return mn,other #return merge node

def mergenodes(mns,mark=True):
    logging.debug("Merging nodes %s"%str(mns))
    global o
    ri=0
    if reference!=None:
        for i,node in enumerate(mns):
            if reference in G.node[node]['offsets'].keys():
                refnode=node
                ri=i
                break
        else:
            refnode=mns[ri]
    else:
        refnode=mns[ri]
    
    newoffsets=dict()
    for node in mns:
        d=G.node[node]
        for sample in d['offsets']:
            if sample in newoffsets:
                logging.warn("WARNING: merging nodes that originate from the same sample: %s in %s."%(sample,str(newoffsets.keys())))
            #assert(sample not in newoffsets)
            newoffsets[sample]=d['offsets'][sample]
    
    G.node[refnode]['offsets']=newoffsets
    assert(len(G.node[refnode]['offsets'])==len(newoffsets))
    
    if mark:
        o+=1 #increment counter, to be able to keep track of when nodes were aligned 
        G.node[refnode]['aligned']=o
    
    tmp=mns.pop(ri)
    assert(tmp==refnode)
    
    for mn in mns: #leave the first node, merge the rest
        for e0,e1,d in G.in_edges(mn,data=True):
            assert(not G.has_edge(refnode,e0))
            if G.has_edge(e0,refnode):
                for p in d['paths']:
                    G[e0][refnode]['paths'].add(p)
            else:
                G.add_edge(e0,refnode,**d)

        for e0,e1,d in G.out_edges(mn,data=True):
            assert(not G.has_edge(e1,refnode))
            if G.has_edge(refnode,e1):
                for p in d['paths']:
                    G[refnode][e1]['paths'].add(p)
            else:
                G.add_edge(refnode,e1,**d)
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
    logging.debug("In graphalign.")

    nodes=index.nodes
    isize=index.n
    
    if len(nodes)==0:
        logging.debug("Invalid set of nodes (length=%d, samples=%d, score=%d, penalty=%d, sp=%s, indexsize=%d)"%(l,n,score,penalty,sp,index.n))
        return
    
    if l==0:
        logging.debug("Invalid length (length=%d, samples=%d, score=%d, penalty=%d, sp=%s, indexsize=%d)"%(l,n,score,penalty,sp,index.n))
        return

    if schemes.minscore!=None:
        if score<schemes.minscore:
            if schemes.pcutoff!=None: #when too much deviation from diagonal, check validity of a match, by approximating a pvalue based on a uniform random distribution of bases, use threshold
                asl=((isize/len(nodes))-l) #assume equal length sequence
                if asl<0:
                    asl=1
                p=float((.25**(n-1))**l)
                npos=(asl**n)
                cp=(npos*(p*((1-p)**(npos-1))))# / ((1-p)**(npos))
                logging.error("cp=%.4g -- p=%.4g (l=%d,score=%d,n=%d,samples=%d,npos=%d,asl=%d,indexsize=%d,cutoff=%.4g)"%(cp,p,l,score,n,len(nodes),npos,asl,isize,schemes.pcutoff))
                if cp>schemes.pcutoff:
                    logging.error("Reject MUM, pvalue=%.4g (l=%d,n=%d,samples=%d,indexsize=%d,cutoff=%.4g)"%(cp,l,n,len(nodes),isize,schemes.pcutoff))
                    return
                else:
                    logging.error("Accept MUM, pvalue=%.4g (l=%d,n=%d,samples=%d,indexsize=%d,cutoff=%.4g)"%(cp,l,n,len(nodes),isize,schemes.pcutoff))
            else:
                logging.debug("Reject MUM, score too low (length=%d, samples=%d, score=%d, penalty=%d, sp=%s, indexsize=%d)"%(l,n,score,penalty,sp,index.n))
                return
    
    if l<schemes.minlength:
        logging.debug("Reject MUM, too short (length=%d, samples=%d, score=%d, penalty=%d, sp=%s, indexsize=%d)"%(l,n,score,penalty,sp,index.n))
        return
    
    logging.debug("Align graph to MUM of length %d (samples=%d, score=%d, penalty=%d, sp=%s, indexsize=%d)"%(l,n,score,penalty,sp,isize))
    
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
    
    logging.debug("Merging samples: %s"%str(msamples))

    intervals=segmentgraph(mn,nodes)
    
    leading,trailing,rest, merged=intervals
    
    logging.debug("Merged interval: %s"%str(merged))
    logging.debug("Number of leading intervals: %d"%len(leading))
    logging.debug("Number of trailing intervals: %d"%len(trailing))
    logging.debug("Number of parallel intervals: %d"%len(rest))
    logging.debug("Number of nodes in the entire graph: %d"%G.number_of_nodes())
    
    #G.node[mn]['penalty']=penalty
    #G.node[mn]['score']=score
    
    newleft=mn
    newright=mn
    
    for intv in leading:
        if not set(G.node[Interval(intv[0],intv[1])]['offsets'].keys()).issubset(msamples): #no clean dissection of all paths on the left
            newright=index.right
            break
    
    for intv in trailing:
        if not set(G.node[Interval(intv[0],intv[1])]['offsets'].keys()).issubset(msamples): #no clean dissection of all paths on the right
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
    converged=False
    while not(converged):
        converged=True
        for node,data in G.nodes_iter(data=True):
            if node not in G:
                continue
            
            for run in [0,1]:
                if 'aligned' in data:
                    if data['aligned']!=0:
                        if run==0:
                            neis=G.successors(node)
                        else:
                            neis=G.predecessors(node)
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
                                    mergenodes(group,mark=True)
                                    converged=False

def align_seq(s1,s2,minlength=1,minscore=None,minn=2,sa64=False):
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
    G.add_node(Intv,offsets={"s1":0},aligned=0)

    idx.addsample("s2")
    intv=idx.addsequence(s2.upper())
    Intv=Interval(intv[0],intv[1])
    t.add(Intv)
    G.add_node(Intv,offsets={"s2":0},aligned=0)
    
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
        pref=[]
        for f in args.inputfiles:
            bn=os.path.basename(f)
            if '.' in bn:
                pref.append(bn[:bn.find('.')])
            else:
                pref.append(bn)
        args.output="_".join(pref)
    
    logging.info("Merging nodes...")
    T=idx.T

    if len(G.graph['samples'])>2:
        prune_nodes(G,T) #TODO: do this after every alignment step to reduce memory usage of graph
    
    logging.info("Done.")
    
    alignedbases=0
    alignednodes=0
    
    totnodes=G.number_of_nodes()
    if idx.nsamples>2: #was multi-alignment
        totbases=idx.n-T.count('$') #TODO: inefficient, need a count of the number of nodes before the alignment
        for node,data in G.nodes(data=True):
            if data['aligned']!=0:
                alignedbases+=(node.end-node.begin)*len(data['offsets'])
                alignednodes+=1
    else: #assume seq to graph
        totbases=min([(idx.n-1)-(idx.nsep[0]+1),idx.nsep[0]])
        
        for node,data in G.nodes(data=True):
            if data['aligned']!=0:
                l=node.end-node.begin
                alignedbases+=l
                alignednodes+=1

    logging.info("%s (%.2f%% identity, %d bases out of %d aligned, %d nodes out of %d aligned)."%("-".join([os.path.basename(f) for f in args.inputfiles]), (alignedbases/float(totbases))*100,alignedbases,totbases,alignednodes,totnodes))
    logging.info("Writing graph...")
    if args.gml:
        graph=write_gml(G,T, hwm=args.hwm, outputfile=args.output)
    else:
        write_gfa(G,T,nometa=args.nometa, outputfile=args.output+'.gfa', paths=args.paths)
        graph=args.output+'.gfa'
    logging.info("Done.")
    logging.info("Alignment graph written to: %s"%graph)
    
    if args.mumplot:
        if len(G.graph['samples']==2):
            plotgraph(sg,G.graph['samples'][0],G.graph['samples'][1],interactive=args.interactive)
        else:
            logging.info("Unable to make plot for graphs with more than 2 paths.")
    
    #plotgraph(G,G.graph['samples'][0],G.graph['samples'][1],interactive=args.interactive)

    
def align_genomes(args):
    logging.info("Loading input...")
    #global variables to simplify callbacks from c extension
    global t,G,reference,o
    
    reference=args.reference
    
    t=IntervalTree()

    if args.sa64:
        idx=reveallib64.index(sa=args.sa, lcp=args.lcp, cache=args.cache)
    else:
        idx=reveallib.index(sa=args.sa, lcp=args.lcp, cache=args.cache)
    
    G=nx.DiGraph()
    G.graph['samples']=[]
    G.graph['sample2id']=dict()
    G.graph['id2sample']=dict()

    o=0
    schemes.pcutoff=args.pcutoff
    schemes.minlength=args.minlength
    schemes.minscore=args.minscore
    schemes.minn=args.minn
    schemes.wscore=args.wscore
    schemes.wpen=args.wpen
    
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
            for name,seq in fasta_reader(sample):
                sid=len(G.graph['samples'])
                name=name.replace(":","").replace(";","")
                G.graph['samples'].append(name)
                G.graph['sample2id'][name]=sid
                G.graph['id2sample'][sid]=name
                intv=idx.addsequence(seq.upper())
                logging.debug("Adding interval: %s"%str(intv))
                Intv=Interval(intv[0],intv[1])
                t.add(Intv)
                G.add_node(Intv,offsets={sid:0},aligned=0)
    
    if not nx.is_directed_acyclic_graph(G):
        logging.error("*** Input is not a DAG! ...")
        return
    
    schemes.exp=args.exp
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

def align(aobjs,ref=None,minlength=15,minscore=None,minn=2,threads=0,global_align=False,targetsample=None,maxsamples=None,wpen=1,wscore=3,sa64=False):
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
    G.graph['sample2id']=dict()
    G.graph['id2sample']=dict()
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
                G.graph['sample2id'][name]=len(G.graph['samples'])
                G.graph['id2sample'][len(G.graph['samples'])]=name
                G.graph['samples'].append(name)
                G.add_node(Intv,offsets={G.graph['sample2id'][name]:0},aligned=0)
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
                        G.add_node(Intv,offsets={os.path.basename(sample):0},aligned=0)
    
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

