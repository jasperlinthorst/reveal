import networkx as nx
import utils
import sys
import logging

def extract_cmd(args):
    if not args.graph[0].endswith(".gfa"):
        logging.fatal("Invalid gfa file.")
        return
    width=args.width
    G=nx.MultiDiGraph()
    utils.read_gfa(args.graph[0], None, None, G, remap=False)

    if args.all:
        args.input=G.graph['paths']
    
    try:
        i=0
        for ins in args.input:
            if args.type=="pathname":
                seq=extract(G,ins)
            elif args.type=="path":
                seq=extract_path(G,ins.split(","))
            else:
                logging.fatal("Unknown input type")
                sys.exit(1)
            
            sys.stdout.write(">"+ins+"\n")
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
    logging.debug("Extracting path: %s"%sample)

    if sample not in G.graph['path2id']:
        logging.fatal("Unknown path: %s, graph contains: %s"%(sample, G.graph['path2id'].keys()))
        sys.exit(1)

    sid=G.graph['path2id'][sample]
    sg=[]

    for n1,n2,d in G.edges(data=True):
        if sid in d['paths']:
            sg.append((n1,n2,d))

    if len(sg)>0:
        #G can be a MultiDiGraph, but subgraph should be single edge!
        subgraph=nx.DiGraph(sg)
        seq=""
        path=list(nx.topological_sort(subgraph))
        subgraph.add_edge(0,path[0],ofrom='+',oto='+')

        pnode=0
        for node in path:
            o=subgraph[pnode][node]['oto']
            if o=="+":
                seq+=G.node[node]['seq'] # --> seq content should come from the original graph
            else:
                seq+=utils.rc(G.node[node]['seq'])
            pnode=node

    else: #has to be a single node
        seq=""
        for n,d in G.nodes(data=True):
            if sid in d['offsets']:
                seq=d['seq']
                break

    return seq

def extract_path(G,path):
    logging.debug("Extracting path: %s"%path)

    seq=""
    for n in path:
        nid,o=int(n[:-1]),n[-1:]
        assert(o=='+' or o=='-')

        if o=="+":
            seq+=G.node[nid]['seq']
        else:
            seq+=utils.rc(G.node[nid]['seq'])

    return seq