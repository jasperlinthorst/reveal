import reveallib
import reveallib64
from utils import *
from intervaltree import IntervalTree
import networkx as nx
import uuid

def chain_cmd(args):
    fastas=args.fastas
    idx=reveallib.index()
    minn=args.minn
    
    tree=IntervalTree()
    
    for fasta in fastas:
        sample=os.path.basename(fasta)
        idx.addsample(sample)
        for i,t in enumerate(fasta_reader(fasta)):
            name,seq=t
            f,t=idx.addsequence(seq)
            tree[f:t]=sample
            if i==1:
                logging.error("Can't handle multi-fasta input. Use single fasta file per sequence.")
                sys.exit(1)
    
    idx.construct()
    
    G=nx.DiGraph()
    G.graph['paths']=idx.samples
    G.graph['path2id']=dict()
    G.graph['id2path']=dict()
    G.graph['startnodes']=[]
    G.graph['endnodes']=[]

    for sid,sample in enumerate(G.graph['paths']):
        G.graph['path2id'][sample]=sid
        G.graph['id2path'][sid]=sample
    
    k=len(idx.samples)
    
    T=idx.T
    
    istart=tuple([-1]+[sep for sep in idx.nsep]) #no matches possible at these loci
    iend=tuple([sep for sep in idx.nsep]+[idx.n-1]) #loci of sentinels, also no matches possible
    startcoords=tuple([0]+[sep+1 for sep in idx.nsep])
    G.add_node(istart,l=0)
    G.add_node(iend,l=0)
    G.add_edge(istart,iend)
    
    G.graph['startnodes'].append(istart)
    G.graph['endnodes'].append(iend)
    
    idc=range(idx.nsamples)

    stack=[(idx,idc,istart,iend,startcoords,0,False)]
    
    while len(stack)!=0:
        idx,idc,p1,p2,startcoords,depth,keepedge=stack.pop()
        subg,pp1,pp2,nodepath=chain(idx,startcoords,args.minlength,depth,args.maxmums,recurse=args.recurse,uniq=True,gcmodel=args.gcmodel,wpen=args.wpen,wscore=args.wscore)

        if len(nodepath)==2: #no more chain, output variant sequence
            localstart=tuple([-1]+[sep for sep in idx.nsep])
            localend=tuple([sep-1 for sep in idx.nsep]+[idx.n-2])
            lengths=tuple([e-s for s,e in zip(localstart,localend)])
            outputVariantNodes(G,T,p1,p2,startcoords,lengths)
            if not keepedge:
                G.remove_edge(p1,p2)
            continue
         
        #replace the edge (start,end) in G with the chain in subg
        insertSubgraph(G,p1,p2,subg,pp1,pp2,keepedge)
        
        coordpath=list(nodepath)
        coordpath[0]=tuple([d+1 for d in nodepath[0]])
        nodepath[0]=p1
        nodepath[-1]=p2

        fromcoord=coordpath[0]
        fromnode=nodepath[0]
        l=0
        
        #for every edge in subg construct idx and add to stack
        for node,pos in zip(nodepath[1:],coordpath[1:]):
            seq=[]
            idc_=[]
            keepedge=False

            for i in idc:
                f=fromcoord[i]
                t=pos[i]
                assert(f>=0)
                assert(t>=0)
                if f+l<t:
                    seq.append(T[f+l:t])
                    idc_.append(i)
                elif f+l==t:
                    keepedge=True
                else:
                    print "Error overlapping matches",f,l,t
                    sys.exit(1)
            
            if len(seq)>=minn and args.recurse==True:
                idx=reveallib.index()
                for i,s in enumerate(seq):
                    assert('$' not in s)
                    idx.addsample(str(i))
                    idx.addsequence(s)
                idx.construct()
                
                newoffsets=tuple([fromcoord[i]+l for i in idc_])
                idc_=range(len(newoffsets))
                stack.append((idx, idc_, fromnode, node, newoffsets, depth+1, keepedge))
            else:
                varnodes=[fromcoord[i]+l for i in idc_]
                lengths=[pos[i]-(fromcoord[i]+l) for i in idc_]
                outputVariantNodes(G,T,fromnode,node,varnodes,lengths)
                if not keepedge:
                    G.remove_edge(fromnode,node)
            
            fromcoord=pos
            fromnode=node
            
            if node!=nodepath[-1]:
                l=subg.node[node]['l']
    
    G.remove_node(istart)
    G.remove_node(iend)
    
    tot=0
    totn=0
    for node,data in G.nodes(data=True):
        G.node[node]['offsets']=dict()
        
        if isinstance(node,tuple):
            G.node[node]['seq']=T[node[0]:node[0]+data['l']]
            for c in node:
                intv=list(tree[c])[0]
                G.node[node]['offsets'][G.graph['path2id'][intv[2]]]=c-intv[0]
        else:
            if 'l' in data:
                G.node[node]['seq']=T[node:node+data['l']]
            intv=list(tree[node])[0]
            G.node[node]['offsets'][G.graph['path2id'][intv[2]]]=node-intv[0]
        
        if 'aligned' in data:
            if data['aligned']==1:
                tot+=data['l']
                totn+=1

    print "Aligned",tot,"bases in",totn,"nodes. Nodes total:",G.number_of_nodes(),"Edges total:",G.number_of_edges()
    
    if args.mumplot:
        plotgraph(G, G.graph['paths'][0], G.graph['paths'][1], interactive=args.interactive)
    
    if args.output==None:
        pref=[]
        for f in args.fastas:
            bn=os.path.basename(f)
            if '.' in bn:
                pref.append(bn[:bn.find('.')])
            else:
                pref.append(bn)
        args.output="_".join(pref)
    
    #add paths annotation to edges
    for sample in G.graph['paths']:
        sid=G.graph['path2id'][sample]
        sg=[]
        for node,data in G.nodes(data=True):
            if sid in data['offsets']:
                sg.append(node)
        subgraph=G.subgraph(sg)
        topsort=list(nx.topological_sort(subgraph))
        pnode=topsort[0]
        for node in topsort[1:]:
            if 'paths' in G[pnode][node]:
                G[pnode][node]['paths'].add(sid)
            else:
                G[pnode][node]['paths']={sid}
            pnode=node

    write_gfa(G,T,nometa=args.nometa,outputfile=args.output+'.gfa')

def outputVariantNodes(G,T,source,sink,varnodes,lengths,merge=True):
    if merge:
        seq=[]
        uvarseq=dict()
        gaps=[]
        for n,l in zip(varnodes,lengths):
            s=T[n:n+l]
            if 'N' not in s: #dont merge gaps, biases statistics
                if s in uvarseq:
                    uvarseq[s]+=[n]
                else:
                    uvarseq[s]=[n]
            else:
                gaps.append((n,l))
        
        for uv in uvarseq:
            G.add_node(tuple(uvarseq[uv]),l=len(uv),aligned=1 if len(uvarseq[uv])>1 else 0)
            G.add_edge(source,tuple(uvarseq[uv]))
            G.add_edge(tuple(uvarseq[uv]),sink)
        
        for gap,l in gaps:
            G.add_node(tuple([gap]),l=l,aligned=0)
            G.add_edge(source,tuple([gap]))
            G.add_edge(tuple([gap]),sink)
    else:
        for v,l in zip(varnodes,lengths):
            G.add_node(v,l=l)
            G.add_edge(source,v)
            G.add_edge(v,sink)

def chain(idx,offsets,minlength,depth,maxmums,recurse=True,uniq=True,gcmodel="sumofpairs",wpen=1,wscore=1):
    k=idx.nsamples
    
    if k>2:
        logging.debug("Extracting mmums of length %d."%minlength)
        mums=idx.getmultimums(minlength=minlength,minn=k)
        logging.debug("Extracted %d mums."%len(mums))
    else:
        logging.debug("Extracting mums of length %d."%minlength)
        mums=idx.getmums(minlength)
        logging.debug("Extracted %d mums."%len(mums))
    
    points=[]
    G=nx.DiGraph()
    localoffsets=tuple([0]+[sep+1 for sep in idx.nsep])
    localstart=tuple([-1]+[sep for sep in idx.nsep])
    localend=tuple([sep-1 for sep in idx.nsep]+[idx.n-2])
    lengths=tuple([e-s for s,e in zip(localstart,localend)])

    p1=tuple([o-1 for o in offsets])
    p2=tuple([o+l for o,l in zip(offsets,lengths)])
    
    mums=[m for m in mums if m[1]==k] #filter only mums that occur in all genomes
    
    if len(mums)>maxmums and maxmums!=0:
        logging.info("Capping the %d anchors that were detected, taking the maxmums=%d longest for chaining."%(len(mums),maxmums))
        mums=sorted(mums,key=lambda m: m[0])[-maxmums:] #take top n longest mums
    elif len(mums)==0:
        return G,p1,p2,[p1,p2]
    else:
        #print "Found %d anchors."%len(mums)
        mums=sorted(mums,key=lambda m: m[0])
    
    logging.debug("Computing chain for %d mums."%len(mums))
    #add all nodes to the graph
    for mum in mums:
        # point=sorted(mum[2])
        # point=sorted(mum[2].values())
        point=sorted([sp for gid,sp in mum[2]])
        for i,p in enumerate(point):
            point[i]=offsets[i]+(point[i]-localoffsets[i]) #map positions back to toplevel T index
        point=tuple(point)
        points.append(point)
        G.add_node(point,l=mum[0])
    
    G.add_node(p1,s=0,l=0,score=0)
    G.add_node(p2,s=0,l=0,score=0)
    
    points.append(p2)
    points=sorted(points,key=lambda p:p[0]) #sort points by first dimension
    
    #build the k-dimensional tree for fast k-dimensional range queries
    tree=kdtree(points,k)
    
    #add edges to graph
    for t in points:
        bestpoint=p1
        bestpenalty=gapcost(p1,t,model=gcmodel)
        bestscore=-1*(wpen*bestpenalty)
        for v in range_search(tree,p1,t):
            if v==t:
                continue
            l=G.node[v]['l']
            for i,d in enumerate(v): #no overlapping mums
                if d+l>t[i]:
                    break
            else:
                penalty=gapcost(v,t,model=gcmodel)
                score=G.node[v]['score']+(wscore*(l*((len(v)*(len(v)-1))/2)))-(wpen*penalty)
                if score>bestscore:
                    bestscore=score
                    bestpoint=v
                    bestpenalty=penalty

        G.node[t]['score']=bestscore
        G.add_edge(bestpoint,t,p=bestpenalty)
    
    logging.debug("backtrack")
    bestpath=[]
    #backtrack the optimal path
    v=p2
    while v!=p1:
        bestpath.append(v)
        G.node[v]['aligned']=1
        for v in G.predecessors(v):
            break

    bestpath.append(p1)
    
    logging.debug("remove nodes")
    #remove nodes that aren't part of the optimal path
    delete=[]
    bestpaths=set(bestpath)
    for node in G.nodes():
        if node not in bestpaths:
            delete.append(node)
    #remove nodes that are not contained in the bestpath
    G.remove_nodes_from(delete)
    logging.debug("done")

    return G,p1,p2,bestpath[::-1]

def insertSubgraph(G,start,end,subg,sstart,send,keepedge):
    upref=uuid.uuid4().hex
    mapping={ sstart : upref+str(sstart), send : upref+str(send) }
    nx.relabel_nodes(subg,mapping,copy=False) #relabel the start and end node in the subgraph to prevent overlap in id space
    
    for node in subg.nodes():
        if node in G:
            print "node already exists!",node
        assert(node not in G)
    
    G.add_nodes_from(subg.nodes(data=True))
    G.add_edges_from(subg.edges(data=True))

    for nei in G.successors(upref+str(sstart)):
        G.add_edge(start,nei)
    
    for nei in G.predecessors(upref+str(send)):
        G.add_edge(nei,end)
    
    if not keepedge:
        G.remove_edge(start,end)

    G.remove_node(upref+str(send))
    G.remove_node(upref+str(sstart))

