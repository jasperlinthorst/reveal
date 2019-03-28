import logging
import traceback

from utils import *

from collections import defaultdict, deque
import os

import reveallib
import reveallib64
import schemes
import bubbles

def breaknode(node,pos,l):
    att=G.node[node]

    in_edges=[(e[0],e[2]) for e in G.in_edges(node,data=True)]
    out_edges=[(e[1],e[2]) for e in G.out_edges(node,data=True)]

    mn=Interval(pos,pos+l)
    other=set()
    
    if mn==node: #no breaking needed
        logging.log(1,"Node %s does not need to be broken."%str(node))
        t.remove(node)
        return node,other
    
    logging.log(1,"Breaking node: %s into: %s"%(str(node),str(mn)))

    allpaths=set()
    
    moffsets=dict()
    for s in att['offsets']:
        moffsets[s]=att['offsets'][s]+(pos-node.begin)
        allpaths.add(s)
    
    logging.log(1,"Offsets after break: %s"%str(moffsets))

    soffsets=dict()
    for s in att['offsets']:
        soffsets[s]=att['offsets'][s]+((pos+l)-node.begin)
    
    #if node is traversed via the other strand, add reverse edges
    negstrand=False
    negpaths=set()
    pospaths=set()
    
    if len(in_edges)==0 and len(out_edges)==0:
        pospaths=allpaths
    else:
        if len(in_edges)>0:
            for fro,d in in_edges:
                if d['oto']=='-':
                    negstrand=True
                    for p in d['paths']:
                        negpaths.add(p)
                else:
                    assert(d['oto']=='+')
                    for p in d['paths']:
                        pospaths.add(p)
        if len(out_edges)>0:
            for to,d in out_edges:
                if d['ofrom']=='-':
                    negstrand=True
                    for p in d['paths']:
                        negpaths.add(p)
                else:
                    assert(d['ofrom']=='+')
                    for p in d['paths']:
                        pospaths.add(p)

    if pospaths.intersection(negpaths)!=set():
        logging.error("Unable to properly separate paths through node: %s [%s,%s,%s,%s]"%(node, G.graph['id2path'], allpaths, pospaths, negpaths))

    assert(pospaths.intersection(negpaths)==set())

    G.add_node(mn,offsets=moffsets,aligned=0)#create merge node
    
    if (node[0]!=pos):
        pn=Interval(node[0],pos)
        logging.log(1,"Creating prefix node: %s"%str(pn))
        G.add_node(pn,offsets=att['offsets'],aligned=0)#create prefix node
        assert(not G.has_edge(pn,mn))
        assert(not G.has_edge(mn,pn))
        G.add_edge(pn,mn,paths=pospaths.copy(),ofrom='+',oto='+')
        assert(pospaths!=set())
        if negstrand:
            G.add_edge(mn,pn,paths=negpaths.copy(),ofrom='-',oto='-')
            assert(negpaths!=set())
        t.add(pn)
        other.add(pn)
    else:
        pn=mn

    if (node[1]!=pos+l):
        sn=Interval(pos+l,node[1])
        logging.log(1,"Creating suffix node: %s"%str(sn))
        G.add_node(sn,offsets=soffsets,aligned=0)#create suffix node
        assert(not G.has_edge(mn,sn))
        assert(not G.has_edge(sn,mn))
        G.add_edge(mn,sn,paths=pospaths.copy(),ofrom='+',oto='+')
        assert(pospaths!=set())
        if negstrand:
            G.add_edge(sn,mn,paths=negpaths.copy(),ofrom='-',oto='-')
            assert(negpaths!=set())
        t.add(sn)
        other.add(sn)
    else:
        sn=mn

    G.remove_node(node)                     #update Graph
    t.remove(node)                          #update intervaltree

    for fro,d in in_edges:
        if d['oto']=="+":
            G.add_edge(fro,pn,**d)
        else:
            G.add_edge(fro,sn,**d)

    for to,d in out_edges:
        if d['ofrom']=="+":
            G.add_edge(sn,to,**d)
        else:
            G.add_edge(pn,to,**d)
    
    logging.log(1,"Leading/Trailing node(s): %s"%str(other))
    logging.log(1,"Matching node: %s"%str(mn))

    return mn,other #return merge node

def mergenodes(G,mns):
    # logging.trace("Merging nodes %s"%str(mns))
    
    # global o
    ri=0
    # if reference!=None:
    #     for i,node in enumerate(mns):
    #         if reference in G.node[node]['offsets'].keys():
    #             refnode=node
    #             ri=i
    #             break
    #     else:
    #         refnode=mns[ri]
    # else:
    refnode=mns[ri]

    #merge the offset dictionaries
    newoffsets=dict()
    for node in mns:
        d=G.node[node]
        for sampleid in d['offsets']:
            if sampleid in newoffsets:
                logging.warn("WARNING: merging nodes that originate from the same sample: %s in %s."%(sample,str(newoffsets.keys())))
            #assert(sample not in newoffsets)
            newoffsets[sampleid]=d['offsets'][sampleid]
    
    G.node[refnode]['offsets']=newoffsets
    assert(len(G.node[refnode]['offsets'])==len(newoffsets))
    
    G.node[refnode]['aligned']=1
    
    tmp=mns.pop(ri)
    assert(tmp==refnode)

    for mn in mns: #leave the first node, merge the rest

        if type(G)==nx.MultiDiGraph:
            for e0,e1,k,d in G.in_edges(mn,keys=True,data=True):
                for _e0,_e1,_k,_d in G.in_edges(refnode,keys=True,data=True):
                    if type(_e0)==type(e0) and _e0==e0 and _d['oto']==d['oto'] and _d['ofrom']==d['ofrom']: #edge already exists, merge paths
                        for p in d['paths']:
                            G[_e0][_e1][_k]['paths'].add(p)
                        break
                else:
                    G.add_edge(e0,refnode,**d)

            for e0,e1,k,d in G.out_edges(mn,keys=True,data=True):
                for _e0,_e1,_k,_d in G.out_edges(refnode,keys=True,data=True):
                    if type(_e1)==type(e1) and _e1==e1 and _d['oto']==d['oto'] and _d['ofrom']==d['ofrom']: #edge already exists, merge paths
                        for p in d['paths']:
                            G[_e0][_e1][_k]['paths'].add(p)
                        break
                else:
                    G.add_edge(refnode,e1,**d)
        else:
            for e0,e1,d in G.in_edges(mn,data=True):
                if G.has_edge(e0,refnode):
                    for p in d['paths']:
                        G[e0][refnode]['paths'].add(p)
                else:
                    G.add_edge(e0,refnode,**d)
            for e0,e1,d in G.out_edges(mn,data=True):
                if G.has_edge(refnode,e1):
                    for p in d['paths']:
                        G[refnode][e1]['paths'].add(p)
                else:
                    G.add_edge(refnode,e1,**d)

        G.remove_node(mn)
    
    return refnode

def predecessorsfilter_iter(G,node):
    if type(G)==nx.MultiDiGraph:
        for pre in G.predecessors(node):
            for i in G[pre][node]:
                for p in G[pre][node][i]['paths']:
                    if not G.graph['id2path'][p].startswith("*"):
                        yield pre
                        break
    else:
        for pre in G.predecessors(node):
            for p in G[pre][node]['paths']:
                if not G.graph['id2path'][p].startswith("*"):
                    yield pre
                    break

def successorsfilter_iter(G,node):
    if type(G)==nx.MultiDiGraph:
        for suc in G.successors(node):
            for i in G[node][suc]:
                for p in G[node][suc][i]['paths']:
                    if not G.graph['id2path'][p].startswith("*"):
                        yield suc
                        break
    else:
        for suc in G.successors(node):
            for p in G[node][suc]['paths']:
                if not G.graph['id2path'][p].startswith("*"):
                    yield suc
                    break

def bfs(G, source, reverse=False, ignore=set()):
    if reverse:
        neighbors = predecessorsfilter_iter
    else:
        neighbors = successorsfilter_iter
    visited = set([source])
    queue = deque([(source, neighbors(G,source))])
    while queue:
        parent, children = queue[0]
        try:
            child = next(children)
            if child not in visited:
                visited.add(child)
                if 'aligned' not in G.node[child]:
                    assert(type(child)==str) #has to be start or end node
                    yield child,2
                elif (G.node[child]['aligned']==0):
                    queue.append((child, neighbors(G,child)))
                    yield child,0
                elif (G.node[child]['aligned']!=0 and child in ignore): #keep searching
                    queue.append((child, neighbors(G,child)))
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
    fwdmerge=set()
    for c,t in bfs(G,node):
        if t==0:
            trailing.add(c)
        elif t==1 or t==2:
            endpoints.add(c)
            if t==2:
                fwdmerge.add(c)
        else:
            logging.error("Node traversal failed, encountered: %s"%str((c,t)))
            sys.exit(1)
    
    #reverse search for each endpoint
    if len(endpoints)>1:
        for endpoint in endpoints:
            for c,t in bfs(G,endpoint,reverse=True,ignore=endpoints):
                if t==0:
                    reverse_trailing.add(c)
        trailing=trailing.intersection(reverse_trailing)
    
    #backward search
    endpoints=set()
    bwdmerge=set()
    for c,t in bfs(G,node,reverse=True):
        if t==0:
            leading.add(c)
        elif t==1 or t==2:
            endpoints.add(c)
            if t==2:
                bwdmerge.add(c)
        else:
            logging.error("Node traversal failed, encountered: %s"%str((c,t)))
            sys.exit(1)
    
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
    
    return leading, trailing, rest

def graphalign(index,mum):
    try:
        logging.debug("In graphalign with %s"%str(mum))
        l,n,spd=mum
        nodes=index.nodes
        isize=index.n
        mns=[]
        topop=[]
        logging.debug("Nodes in subgraph:")
        for node in nodes:
            logging.debug("%s"%str(node))

        sp=[sp for gid,sp in spd]
        matching=set()
        for pos in sp:
            matching.add((pos,pos+l))
            # logging.debug("Lookup node for sp=%d"%pos)
            old=t[pos].pop()
            assert(old.end-old.begin>=l)
            mn,other=breaknode(old,pos,l)
            mns.append(mn)
            if isinstance(old,Interval):
                nodes.remove((old.begin,old.end))
            for node in other:
                if isinstance(node,Interval):
                    nodes.add((node.begin,node.end))
        
        mn=mergenodes(G,mns)
        msamples=set(G.node[Interval(mn[0],mn[1])]['offsets'].keys())
        # logging.trace("Merging samples: %s"%str(msamples))
        # logging.trace("Nodes before segmenting: %s"%nodes)

        intervals=segmentgraph(mn,nodes)
        leading,trailing,rest=intervals

        # logging.trace("Leading nodes after segmenting: %s"%leading)
        # logging.trace("Trailing nodes after segmenting: %s"%trailing)
        # logging.trace("Parallel nodes after segmenting: %s"%rest)

        logging.debug("Merged interval: %s"%str(mn))
        logging.debug("Number of leading intervals: %d"%len(leading))
        logging.debug("Number of trailing intervals: %d"%len(trailing))
        logging.debug("Number of parallel intervals: %d"%len(rest))
        # logging.trace("Number of nodes in the entire graph: %d"%G.number_of_nodes())
        newleftnode=mn
        newrightnode=mn

        for intv in leading:
            if not set(G.node[Interval(intv[0],intv[1])]['offsets'].keys()).issubset(msamples): #no clean dissection of all paths on the left
                newrightnode=index.rightnode
                break
        
        for intv in trailing:
            if not set(G.node[Interval(intv[0],intv[1])]['offsets'].keys()).issubset(msamples): #no clean dissection of all paths on the right
                newleftnode=index.leftnode
                break

        return leading,trailing,matching,rest,mn,newleftnode,newrightnode

    except Exception:
        print traceback.format_exc()
        raise Exception
        return

#TODO: rewrite so this uses bubble definition code
def prune_nodes(G,T=""):
    converged=False
    while not(converged):
        converged=True
        for node,data in G.nodes(data=True):
            if node not in G:
                continue

            for run in [0,1]:
                # if type(node)==str or data['aligned']!=0:
                if run==0:
                    if type(G)==nx.MultiDiGraph:
                        neis=[n2 for n1,n2,k,d in G.out_edges(node,keys=True,data=True) if d['ofrom']=='+' and d['oto']=='+']
                    else:
                        neis=[n2 for n1,n2,d in G.out_edges(node,data=True) if d['ofrom']=='+' and d['oto']=='+']
                else:
                    if type(G)==nx.MultiDiGraph:
                        neis=[n1 for n1,n2,k,d in G.in_edges(node,keys=True,data=True) if d['ofrom']=='+' and d['oto']=='+']
                    else:
                        neis=[n1 for n1,n2,d in G.in_edges(node,data=True) if d['ofrom']=='+' and d['oto']=='+']
                
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

                for key in seqs:
                    group=seqs[key]
                    if len(group)>1:
                        merge=True
                        for v in group:
                            if run==0:
                                if type(G)==nx.MultiDiGraph:
                                    if len([n1 for n1,n2,k,d in G.in_edges(v,keys=True,data=True) if d['ofrom']=='+' and d['oto']=='+'])>1:
                                        merge=False
                                        break
                                else:
                                    if len([n1 for n1,n2,d in G.in_edges(v,data=True) if d['ofrom']=='+' and d['oto']=='+'])>1:
                                        merge=False
                                        break
                            else:
                                if type(G)==nx.MultiDiGraph:
                                    if len( [n2 for n1,n2,k,d in G.out_edges(v,keys=True,data=True) if d['ofrom']=='+' and d['oto']=='+'] )>1:
                                        merge=False
                                        break
                                else:
                                    if len( [n2 for n1,n2,d in G.out_edges(v,data=True) if d['ofrom']=='+' and d['oto']=='+'] )>1:
                                        merge=False
                                        break
                        if merge:
                            mergenodes(G,group)
                            converged=False

def align_cmd(args):
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
        args.output+=".gfa"
    
    logging.info("Merging nodes...")
    T=idx.T

    if len(G.graph['paths'])>2:
        prune_nodes(G,T=T)

    seq2node(G,T,remap=False)
    
    logging.info("Done.")
    
    alignedbases=0
    alignednodes=0
    
    totnodes=G.number_of_nodes()

    #TODO: start, report identity per sample that was aligned, this does not make sense...
    if idx.nsamples>2: #was multi-alignment
        totbases=idx.n-T.count('$')-T.count('N')
        for node,data in G.nodes(data=True):
            if 'aligned' in data and data['aligned']!=0 and type(node)!=str:
                alignedbases+=(node.end-node.begin)*len([k for k in data['offsets'] if not G.graph['id2path'][k].startswith("*")])
                alignednodes+=1
    else: #assume seq to graph
        totbases=idx.n-T.count('$')-T.count('N') # min([(idx.n-1)-(idx.nsep[0]+1),idx.nsep[0]])
        for node,data in G.nodes(data=True):
            if 'aligned' in data and data['aligned']!=0 and type(node)!=str:
                l=node.end-node.begin
                alignedbases+=(l*2)
                alignednodes+=1
    
    logging.info("%s (%.2f%% identity, %d bases out of %d aligned, %d nodes out of %d aligned)."%("-".join([os.path.basename(f) for f in args.inputfiles]), (alignedbases/(float(totbases)))*100,alignedbases,totbases,alignednodes,totnodes))
    #TODO: end

    logging.info("Writing graph...")
    if args.gml:
        graph=write_gml(G,T, hwm=args.hwm, outputfile=args.output, partition=False)
    else:
        write_gfa(G,T, outputfile=args.output)
        graph=args.output

    logging.info("Done.")
    logging.info("Graph written to: %s"%graph)
    
    if args.mumplot:
        if len(G.graph['paths'])==2:
            plotgraph(G,G.graph['paths'][0],G.graph['paths'][1],interactive=args.interactive)
        else:
            logging.info("Unable to make plot for graphs with more than 2 paths.")
    
    #plotgraph(G,G.graph['paths'][0],G.graph['paths'][1],interactive=args.interactive)

def align_genomes(args):
    logging.info("Loading input...")
    #global variables to simplify callbacks from c extension
    global t,G
    
    # global reference
    # reference=args.reference
    
    t=IntervalTree()

    if args.sa64:
        idx=reveallib64.index(sa=args.sa, lcp=args.lcp, cache=args.cache)
    else:
        idx=reveallib.index(sa=args.sa, lcp=args.lcp, cache=args.cache)
    
    #G=nx.DiGraph()
    G=nx.MultiDiGraph()

    o=0
    schemes.args=args
    
    graph=False
    
    for i,sample in enumerate(args.inputfiles):
        
        if sample.endswith(".gfa"):
            idx.addsample(os.path.basename(sample))
            graph=True

            logging.info("Reading graph: %s ..." % sample)
            if i==0:
                read_gfa(sample,idx,t,G,minsamples=args.minsamples,
                                        maxsamples=args.maxsamples,
                                        targetsample=args.targetsample,
                                        remap=True)
            else:
                read_gfa(sample,idx,t,G,remap=True)

        else: #consider it to be a fasta file
            read_fasta(sample,idx,t,G,contigs=args.contigs,toupper=args.toupper)
    
    logging.debug("Graph contains the following paths: %s"%G.graph['paths'])

    logging.debug("Index contains the following samples: %s"%idx.samples)

    if len(idx.samples)<=1:
        logging.fatal("Specify at least 2 targets to construct alignment. In case of multi-fasta, consider the --nocontigs flag.")
        sys.exit(1)

    if not nx.is_directed_acyclic_graph(G):
        logging.info("*** Input is not a DAG! ...")

    for n1,n2,data in G.edges(data=True):
        assert('paths' in data)

    schemes.ts=t
    schemes.G=G
    
    logging.info("Constructing index...")
    idx.construct()
    
    logging.info("Done.")
    
    if len(args.inputfiles)==2 and not graph:
        logging.info("Constructing pairwise-alignment...")
        idx.align(schemes.graphmumpicker,graphalign,threads=args.threads,wpen=args.wpen,wscore=args.wscore,minl=args.minlength,minn=args.minn)
    else:
        logging.info("Constructing graph-based multi-alignment...")
        idx.align(schemes.graphmumpicker,graphalign,threads=args.threads,wpen=args.wpen,wscore=args.wscore,minl=args.minlength,minn=args.minn)
    
    # from multiprocessing import Process
    # from Queue import Queue
    # main=idx #make sure we keep the main ref count, since it has the reference to T
    # q=Queue()
    # q.put(idx)
    # while not q.empty():
    #     idx=q.get()
    #     if len(args.inputfiles)>2:
    #         multimums=idx.getmultimums(minlength=args.minlength, minn=args.minn)
    #     else:
    #         multimums=idx.mums(args.minlength)
    #     if len(multimums)==0:
    #         continue
    #     ret=schemes.graphmumpicker(multimums,idx)
    #     if ret==None:
    #         continue
    #     else:
    #         splitmum,skipleft,skipright=ret
    #     ret=graphalign(idx,splitmum)
    #     if ret==None:
    #         continue
    #     else:
    #         leading,trailing,matching,rest,merged,newleftnode,newrightnode=ret
    #     ilead,itrail,ipar=idx.splitindex(leading,trailing,matching,rest,merged,newleftnode,newrightnode,skipleft,skipright)
    #     if ilead!=None and ilead.n>1:
    #         q.put(ilead)
    #     if itrail!=None and itrail.n>1:
    #         q.put(itrail)
    #     if ipar!=None and ipar.n>1:
    #         q.put(ipar)

    return G,idx


#seq should be a list of objects that can be (multi-) aligned by reveal:
#   - tuple of the form (name,seq)
def align(aobjs,ref=None,minlength=20,minn=2,seedsize=None,threads=0,targetsample=None,maxsamples=None,\
                maxmums=10000,wpen=1,wscore=1,sa64=False,pcutoff=1e-8,gcmodel="sumofpairs",maxsize=None,\
                trim=True):
    
    kwargs = dict(locals()) #hack the kwargs into a dict so we can pass it to schemes as if it were the argparsed args object
    class dict2class(object):
        def __init__(self, d):
            self.__dict__ = d
    args=dict2class(kwargs)
    schemes.args=args
    
    #global variables to simplify callbacks from c extension
    global t,G

    t=IntervalTree()

    if sa64:
        idx=reveallib64.index()
    else:
        idx=reveallib.index()
    
    G=nx.DiGraph()

    G.graph['paths']=[]
    G.graph['path2id']=dict()
    G.graph['id2path']=dict()
    G.graph['id2end']=dict()
    o=0

    graph=False
    
    startnode=uuid.uuid4().hex
    G.add_node(startnode)
    endnode=uuid.uuid4().hex
    G.add_node(endnode)

    for aobj in aobjs:
        if isinstance(aobj,tuple):
            name,seq=aobj
            idx.addsample(name)
            intv=idx.addsequence(seq.upper())
            if intv[1]-intv[0]>0:
                Intv=Interval(intv[0],intv[1])
                t.add(Intv)
                sid=len(G.graph['paths'])
                G.graph['path2id'][name]=len(G.graph['paths'])
                G.graph['id2path'][sid]=name
                G.graph['id2end'][sid]=len(seq)
                
                # G.node[endnode]['offsets'][sid]=len(seq)
                # G.node[startnode]['offsets'][sid]=0

                G.graph['paths'].append(name)
                G.add_node(Intv,offsets={sid:0},aligned=0)
                G.add_edge(startnode,Intv,paths={sid},ofrom='+',oto='+')
                G.add_edge(Intv,endnode,paths={sid},ofrom='+',oto='+')

        # elif isinstance(aobj,str):
        #     if not os.path.isfile(aobj):
        #         logging.fatal("Not a file, expecting fasta or gfa file.")
        #         return
        #     idx.addsample(os.path.basename(aobj))
        #     if aobj.endswith(".gfa"):
        #         read_gfa(aobj,idx,t,G,targetsample=targetsample,maxsamples=maxsamples)
        #         graph=True
        #     else: #assume a file in fastaformat
        #         for name,seq in fasta_reader(sample):
        #             intv=idx.addsequence(seq.upper())
        #             if intv[1]-intv[0]>0:
        #                 Intv=Interval(intv[0],intv[1])
        #                 t.add(Intv)
        #                 sid=len(G.graph['paths'])
        #                 G.graph['path2id'][name]=len(G.graph['paths'])
        #                 G.graph['id2path'][sid]=name
        #                 G.graph['id2end'][sid]=len(seq)
        #                 G.graph['paths'].append(name)
        #                 G.add_node(Intv,offsets={sid:0},aligned=0)
        #                 G.add_edge(startnode,Intv,paths={sid})
        #                 G.add_edge(endnode,Intv,paths={sid})
    
    if not nx.is_directed_acyclic_graph(G):
        logging.error("*** Input is not a DAG! Not supported.")
        return
    
    schemes.ts=t
    schemes.G=G
    
    idx.construct()
    
    idx.align(schemes.graphmumpicker,graphalign,threads=threads,wpen=wpen,wscore=wscore,minl=minlength,minn=minn)

    prune_nodes(G,T=idx.T)

    G.remove_node(startnode)
    G.remove_node(endnode)

    return G,idx