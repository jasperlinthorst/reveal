import networkx as nx
from utils import *

def subgraph(args):
    if len(args.inputfiles)<=1:
        logging.fatal("Specify 1 gfa file followed by comma separated list of node ids that make up the subgraph.")
        return
    if not args.inputfiles[0].endswith('.gfa'):
        logging.fatal("Specify gfa file as first argument of subgraph subcommand.")
        return
    
    G=nx.DiGraph()
    read_gfa(args.inputfiles[0],None,"",G)
    
    nodes=set()
    
    for arg in args.inputfiles[1:]:
        for node in arg.split(','):
            if node.find('-')!=-1: #then bubble definition
                source,sink=node.split('-')
                for node in bubbles.bubble(G,source,sink).nodes:
                    nodes.add(int(node))
            else:
                nodes.add(int(node))
    
    sg=G.subgraph(nodes)
    if args.gml:
        write_gml(sg,"",outputfile=args.outfile)
    else:
        write_gfa(sg,"",outputfile=args.outfile, remap=False)

