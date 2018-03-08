import utils
import bubbles
import networkx as nx
import os
import logging
import sys

def chop_cmd(args):
    if not args.graph[0].endswith(".gfa"):
        logging.fatal("Invalid gfa file.")
        return

    G=nx.DiGraph()
    utils.read_gfa(args.graph[0], None, None, G, remap=False)
    
    # assert(len(G.edges())==0)
    if args.output==None:
        fof=os.path.splitext(args.graph[0])[0]+".chopped.fasta"
        gof=os.path.splitext(args.graph[0])[0]+".chopped.gfa"
    else:
        fof=args.output+".fasta"
        gof=args.output+".gfa"

    if args.check:
        Gorg=G.copy()

    chop(G,k=args.k)
    
    for node in G.nodes():
        G.node[node]['seq']=G.node[node]['prefix']+G.node[node]['seq']+G.node[node]['suffix']

    utils.write_gfa(G,None,outputfile=gof,remap=False)

    with open(fof,'w') as ff:
        for node in G.nodes():
            name=">"+str(node)+"\n"
            seq=G.node[node]['seq']
            ff.write(name)
            for i in range( (len(seq)/args.lw)+(len(seq) % args.lw > 0)):
                ff.write(seq[i*args.lw:(i+1)*args.lw]+"\n")
    
    if args.check:
        import extract
        r="$".join([G.node[node]['seq'] for node in G])
        for path in Gorg.graph['paths']:
            logging.debug("Check: %s"%path)
            s=extract.extract(Gorg,path)
            for i in xrange(len(s)-args.k):
                if r.find(s[i:i+args.k])==-1:
                    print i,s[i:i+args.k]
                    logging.error("Flat representation does not cover all k-length substrings for %s!"%path)
                    sys.exit(1)

def contract(G,topsort):
    newtopsort=[]
    stretches=[[]]    
    pnode=topsort[0]
    newtopsort=[topsort[0]]
    for i,node in enumerate(topsort[1:]):
        pred=list(G.predecessors(node))
        suc=list(G.successors(pnode))
        if pred==[pnode] and suc==[node]:
            if len(stretches[-1])==0:
                stretches[-1].append(pnode)
            stretches[-1].append(node)
        else:
            if len(stretches[-1])!=0:
                stretches.append([])
            newtopsort.append(node)
        pnode=node

    for stretch in stretches:
        if len(stretch)>0:
            contract_nodes(G,stretch)

    assert(len(newtopsort)==len(set(newtopsort)))

    return newtopsort

def contract_nodes(G,nodes):
    logging.debug("Contract: Contracting nodes: %s"%nodes)
    G.node[nodes[0]]['seq']="".join([G.node[n]['seq'] for n in nodes])
    for n1,n2,data in G.out_edges(nodes[-1],data=True):
        G.add_edge(nodes[0],n2,**data)
    G.remove_nodes_from(nodes[1:])

def duplicate_node(G,node):
    logging.debug("Duplicate: node %s"%node)
    offsets=G.node[node]['offsets']
    seq=G.node[node]['seq']
    es=[]
    duplicates=[]
    predecessors=list(G.predecessors(node))
    successors=list(G.successors(node))
    if len(predecessors)>0 and len(successors)>0:
        for pred in predecessors:
            for suc in successors:
                i=G[pred][node]['paths'].intersection(G[node][suc]['paths'])
                if len(i)>0:
                    G.add_node(G.graph['noffset'],offsets={k:offsets[k] for k in offsets if k in i},seq=seq,prefix="",suffix="") #TODO: prevent contract, by doing so immediately
                    duplicates.append(G.graph['noffset'])
                    es.append((pred,G.graph['noffset'],{'paths':i,'ofrom':G[pred][node]['ofrom'],'oto':G[pred][node]['oto']}))
                    es.append((G.graph['noffset'],suc,{'paths':i,'ofrom':G[node][suc]['ofrom'],'oto':G[node][suc]['oto']}))
                    G.graph['noffset']+=1
    elif len(predecessors)>0:
        for pred in predecessors:
            i=G[pred][node]['paths']
            G.add_node(G.graph['noffset'],offsets={k:offsets[k] for k in offsets if k in i},seq=seq,prefix="",suffix="") #TODO: prevent contract, by doing so immediately
            duplicates.append(G.graph['noffset'])
            es.append((pred,G.graph['noffset'],{'paths':i,'ofrom':G[pred][node]['ofrom'],'oto':G[pred][node]['oto']}))
            G.graph['noffset']+=1
    elif len(successors)>0:
        for suc in successors:
            i=G[node][suc]['paths']
            G.add_node(G.graph['noffset'],offsets={k:offsets[k] for k in offsets if k in i},seq=seq,prefix="",suffix="") #TODO: prevent contract, by doing so immediately
            duplicates.append(G.graph['noffset'])
            es.append((G.graph['noffset'],suc,{'paths':i,'ofrom':G[node][suc]['ofrom'],'oto':G[node][suc]['oto']}))
            G.graph['noffset']+=1
    
    G.remove_node(node)
    G.add_edges_from(es)

    return duplicates

def checkedges(G,k=100):
    for u,v,d in G.edges(data=True):
        d['overlap']=None
    es=[]
    update=True
    while update:
        update=False
        for u,v,d in G.edges(data=True):
            if d['overlap']!=None:
                continue
            if len(G.node[u]['seq'])>=k-1 and len([e for e in G.in_edges(v,data=True) if e[2]['overlap']==None])==1:
                d['overlap']=u
                update=True
                continue #can use k-1 suffix of u as prefix of v
            if len(G.node[v]['seq'])>=k-1 and len([e for e in G.out_edges(u,data=True) if e[2]['overlap']==None])==1:
                d['overlap']=v
                update=True
                continue #can use k-1 prefix of v as suffix of u
    for u,v,d in G.edges(data=True):
        if d['overlap']==None:
            es.append((u,v))
    return es

def chop(G,k=100):
    remove=[]
    for node in G.nodes():
        if type(node)==str:
            remove.append(node)
        else: #add prefix and suffix attributes
            G.node[node]['prefix']=""
            G.node[node]['suffix']=""

    G.remove_nodes_from(remove)
    iteration=1

    es=checkedges(G,k=k)

    while len(es)!=0:
        logging.info("Running iteration %d"%iteration)
        
        #determine subgraph for duplication
        sg=nx.DiGraph(es)
        nodes=list(sg.nodes())
        nodes=[node for node in nodes if (len(sg.in_edges(node))>1 or len(sg.out_edges(node))>1)]# and len(G.node[node]['seq'])<k-1]
        nodes.sort(key=lambda n: len(G.node[n]['seq']))
        d=set()
        dups=[]
        for n in nodes:
            dup=True
            for n1,n2 in sg.in_edges(n):
                if n1 in d:
                    dup=False
                d.add(n1)
            for n1,n2 in sg.out_edges(n):
                if n2 in d:
                    dup=False
                d.add(n2)
            if dup:
                dups.append(n)
        
        for n in dups:
            duplicate_node(G,n)

        contract(G,list(nx.topological_sort(G)))
        es=checkedges(G,k=k)

        iteration+=1

    #all edges can now be extended
    for u,v,d in G.edges(data=True):
        assert(d['overlap']!=None)
        if d['overlap']==u:
            G.node[v]['prefix']=G.node[u]['seq'][-(k-1):]
        else:
            assert(d['overlap']==v)
            G.node[u]['suffix']=G.node[v]['seq'][:k-1]
        d['cigar']=str(k)+"M"

    return G