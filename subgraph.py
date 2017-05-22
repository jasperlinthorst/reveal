import networkx as nx
from utils import *

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

