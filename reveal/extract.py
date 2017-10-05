import networkx as nx
import utils
import sys
import logging

def extract_cmd(args):
    if not args.graph[0].endswith(".gfa"):
        logging.fatal("Invalid gfa file.")
        return
    width=args.width
    #G=nx.DiGraph();
    G=nx.MultiDiGraph()
    utils.read_gfa(args.graph[0], None, None, G, remap=False)
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

    sid=G.graph['sample2id'][sample]
    sg=[]
    for n1,n2,d in G.edges_iter(data=True):
        if sid in d['paths']:
            sg.append((n1,n2,d))

    #G can be a MultiDiGraph, but subgraph should be single edge!
    subgraph=nx.DiGraph(sg)

    seq=""
    path=nx.topological_sort(subgraph)
    pnode=path[0]
    for node in path[1:]:
        o=subgraph[pnode][node]['ofrom']
        if o=="+":
            seq+=G.node[pnode]['seq'] # --> seq content should come from the original graph
        else:
            seq+=utils.rc(G.node[pnode]['seq'])
        if node==path[-1]:#last
            o=subgraph[pnode][node]['oto']
            if o=="+":
                seq+=G.node[node]['seq']
            else:
                seq+=utils.rc(G.node[node]['seq'])
        pnode=node

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