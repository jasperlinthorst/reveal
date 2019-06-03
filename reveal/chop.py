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

    chop(G,k=args.k,extend=args.extend)
    
    logging.debug("Merging node sequence...")
    for node in G.nodes():
        if type(node)==str: #skip start/end nodes
            continue
        G.node[node]['seq']=G.node[node]['prefix']+G.node[node]['seq']+G.node[node]['suffix']
    
    logging.debug("Done.")

    logging.debug("Write overlap graph...")
    utils.write_gfa(G,None,outputfile=gof,remap=False)
    logging.debug("Done.")

    if args.fasta:
        logging.debug("Write corresponding fasta file...")
        with open(fof,'w') as ff:
            for node in G.nodes():
                if type(node)==str: #skip start/end nodes
                    continue
                name=">"+str(node)+"\n"
                seq=G.node[node]['seq']
                ff.write(name)
                for i in range( (len(seq)/args.lw)+(len(seq) % args.lw > 0)):
                    ff.write(seq[i*args.lw:(i+1)*args.lw]+"\n")
        logging.debug("Done.")
    
    if args.check:
        logging.debug("Validate if all k-mers from the original graph are contained in overlap graph...")
        import extract
        r="$".join([G.node[node]['seq'] for node in G])
        for path in Gorg.graph['paths']:
            logging.debug("Check: %s"%path)
            s=extract.extract(Gorg,path)
            for i in xrange(len(s)-args.k):
                if r.find(s[i:i+args.k])==-1:
                    logging.error("Flat representation does not cover all k-length substrings for %s, could not find: %s!"%(path,s[i:i+args.k]))
                    sys.exit(1)
        logging.debug("Done.")

def duplicate_node(G,node):
    if typ(node)==str:
        logging.fatal("Attempt to duplicate end/start node, shouldn't happen. Exit.")
        sys.exit(1)
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
        
        remove=[]
        for u,v,d in G.edges(data=True):
            if d['overlap']!=None:
                continue

            if type(u)==str:
                d['overlap']=u
                continue

            if type(v)==str:
                d['overlap']=v
                continue

            if len(G.node[u]['seq'])>=k-1 and len([e for e in G.in_edges(v)])==1:
                d['overlap']=u
                update=True
                continue #can use k-1 suffix of u as prefix of v

            if len(G.node[v]['seq'])>=k-1 and len([e for e in G.out_edges(u)])==1:
                d['overlap']=v
                update=True
                continue #can use k-1 prefix of v as suffix of u

            #substitute edge by node that contains only pre- and suffix
            if len(G.node[v]['seq'])>=k-1 and len(G.node[u]['seq'])>=k-1:
                nid=G.graph['noffset']
                G.add_node(nid,seq="",offsets=G.nodes[v]['offsets'],prefix="",suffix="")
                G.graph['noffset']+=1
                G.add_edge(u,nid,**G[u][v])
                G[u][nid]['overlap']=u
                G.add_edge(nid,v,**G[u][v])
                G[nid][v]['overlap']=v
                remove.append((u,v))
                update=True

            # if len(G.node[u]['seq'])>=k-1 and len([e for e in G.in_edges(v,data=True) if e[2]['overlap']==None])==1:
            #     d['overlap']=u
            #     update=True
            #     continue #can use k-1 suffix of u as prefix of v

            # if len(G.node[v]['seq'])>=k-1 and len([e for e in G.out_edges(u,data=True) if e[2]['overlap']==None])==1:
            #     d['overlap']=v
            #     update=True
            #     continue #can use k-1 prefix of v as suffix of u

        G.remove_edges_from(remove)

    for u,v,d in G.edges(data=True):
        # if type(u)==str or type(v)==str:
        #     continue
        if d['overlap']==None:
            es.append((u,v))

    return es

def chop(G,k=100,extend=True):
    # remove=[]
    for node in G.nodes():
        if type(node)==str:
            pass
            # remove.append(node)
        else: #add prefix and suffix attributes
            G.node[node]['prefix']=""
            G.node[node]['suffix']=""

    # G.remove_nodes_from(remove)
    iteration=1

    es=checkedges(G,k=k)

    maxiter=1e22

    while len(es)!=0 and iteration<maxiter:
        logging.info("Running iteration %d"%iteration)
        
        #determine subgraph for duplication
        sg=nx.DiGraph(es)

        nodes=list(sg.nodes())

        nodes=[node for node in nodes if (len(sg.in_edges(node))>1 or len(sg.out_edges(node))>1) and type(node)!=str]# and len(G.node[node]['seq'])<k-1]
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
        
        logging.info("Duplicating nodes...")
        for n in dups:
            logging.debug("Duplicating node %d..."%n)
            for dup in duplicate_node(G,n):
                logging.debug("Generated node %d"%dup)
        logging.info("Duplicating done.")

        logging.info("Contracting nodes...")
        
        topsort=list(nx.topological_sort(G))[1:-1]
        topsort=[v for v in topsort if type(v)!=str]
        
        utils.contract(G,topsort)
        logging.info("Contracting done.")

        # es=checkedges(G,k=k)
        # sg=nx.DiGraph(es)
        # for u,v,d in sg.edges(data=True):
        #     if len(sg.in_edges(v))==len(sg.out_edges(u))==1 and len(G.nodes[u]['seq'])>=k-1 and len(G.nodes[v]['seq'])>=k-1:
        #         nid=G.graph['noffset']
        #         assert(nid not in G)
        #         G.add_node(nid,seq="",offsets=G.nodes[v]['offsets'],prefix="",suffix="")
        #         G.graph['noffset']+=1
        #         G.add_edge(u,nid,**G[u][v])
        #         G.add_edge(nid,v,**G[u][v])
        #         G.remove_edge(u,v)

        logging.debug("Checking edges...")
        es=checkedges(G,k=k)

        logging.info("Done. %d unextendable edges remain."%len(es))

        for u,v in es:
            logging.debug("Edge %s,%s can't be extended yet: %s"%(u,v,G[u][v]))

        iteration+=1

    if len(es)>0:
        logging.fatal("Error, maxiterations reached, chop did not converge!")
        sys.exit(1)

    if extend:
        logging.info("Extending nodes with prefix/suffix...")

        #all edges can now be extended
        for u,v,d in G.edges(data=True):
            if type(u)==str or type(v)==str:
                continue

            assert(d['overlap']!=None)
            
            if d['overlap']==u:
                logging.debug("Add prefix to %s"%v)
                assert(G.node[v]['prefix']=="")
                G.node[v]['prefix']=G.node[u]['seq'][-(k-1):]
                #if we give v a prefix, all other incoming edges of v are also affected, so cigar should be increased on those as well 
                # for _u,_v,_d in G.in_edges(v,data=True):
                    # _d['overlap_length']+=(k-1)
            else:
                assert(d['overlap']==v)
                logging.debug("Add suffix to %s"%u)
                assert(G.node[u]['suffix']=="")
                G.node[u]['suffix']=G.node[v]['seq'][:k-1]
                #if we give u a suffix, all other outgoing edges of u are also affected, so cigar should be increased on those as well
                # for _u,_v,_d in G.out_edges(u,data=True):
                    # _d['overlap_length']+=(k-1)

            d['cigar']=str(k-1)+"M"

        logging.info("Done.")
    
    # logging.info("Unzipping bubbles...")
    # #all edges can now be extended
    # for n in G.nodes():
    #     extend=False
    #     for u,v,d in G.out_edges(n,data=True):
    #         if d['overlap']==v: #overlap comes from v, so use prefix of v as suffix for u
    #             if extend:
    #                 logging.error("PROBLEM out!")
    #                 print n,u,v,d
    #             extend=True
    #     prefix=False
    #     for u,v,d in G.in_edges(n,data=True):
    #         if d['overlap']==u: #overlap comes from u, so use suffix of u as prefix for v
    #             if prefix:
    #                 logging.error("PROBLEM in!")
    #                 print n,u,v,d
    #             prefix=True

    return G