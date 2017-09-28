import networkx as nx
import utils
import sys

def extract_cmd(args):
    if not args.graph[0].endswith(".gfa"):
        logging.fatal("Invalid gfa file.")
        return
    width=args.width
    G=nx.DiGraph()
    utils.read_gfa(args.graph[0], None, None, G)
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
        if G.graph['sample2id'][sample] in data['offsets']:
            sg.append(node)
    for i,contig in enumerate(nx.connected_components(G.subgraph(sg).to_undirected())):
        seq=""
        nodecount=0
        for node in nx.topological_sort(G.subgraph(contig)):
            seq+=G.node[node]['seq']
        yield seq
