import networkx as nx
from utils import *
from intervaltree import IntervalTree
import bubbles

def subgraph(args):
    if len(args.inputfiles)<=1:
        logging.fatal("Specify 1 gfa file followed by comma separated list of node ids that make up the subgraph.")
        return

    if not args.inputfiles[0].endswith('.gfa'):
        logging.fatal("Specify gfa file as first argument of subgraph subcommand.")
        return

    G=nx.DiGraph()

    read_gfa(args.inputfiles[0],None,"",G)
    
    topsort=list(nx.topological_sort(G))

    nodes=set()
    
    for arg in args.inputfiles[1:]:
        if arg.find(':')!=-1: #then parse as interval definition: eg chr4:230000-230010
            cds,intv=arg.split(':')
            tree=graph_to_tree(G,cds) #use interval tree for retrieval
            start,stop=intv.split('-')
            
            source=tree[int(start)].pop()[2]
            sink=tree[int(stop)].pop()[2]
            source_idx,sink_idx=0,0
            for i,v in enumerate(topsort):
                if v==source:
                    source_idx=i
                if v==sink:
                    sink_idx=i
                    break
            if sink_idx<=source_idx:
                logging.fatal("Invalid source/sink pair. Topsort indices: %d, %d"%(source_idx,sink_idx))
                return

            for node in topsort[source_idx:sink_idx+1]:
                nodes.add(int(node))

        elif arg.find('-')!=-1: #then bubble definition
            source,sink=arg.split('-')

            for i,v in enumerate(topsort):
                if v==int(source):
                    source_idx=i
                if v==int(sink):
                    sink_idx=i
                    break
            else:
                logging.fatal("Sink node could not be found.")
                sys.exit(1)

            if sink_idx<=source_idx:
                logging.fatal("Invalid source/sink pair. Topsort indices: %d, %d"%(source_idx,sink_idx))
                return
            
            for node in topsort[source_idx:sink_idx+1]:
                nodes.add(int(node))
        else:
            for node in arg.split(','): #assume a comma separated list of nodes
                nodes.add(int(node))
    
    sg=G.subgraph(nodes)
    sg.graph['startnodes']=[]

    for sid in sg.graph['id2path']:
        start=None
        for node in sg:
            if sid in sg.node[node]['offsets']:
                if start==None or sg.node[node]['offsets'][sid]<start:
                    start=sg.node[node]['offsets'][sid]
                    startnode=node
        sg.graph['startnodes'].append(startnode)

    # for sid in sg.graph['id2path']:
    #     end=None
    #     for node in sg:
    #         if sid in sg.node[node]['offsets']:
    #             if end==None or sg.node[node]['offsets'][sid]>end:
    #                 end=sg.node[node]['offsets'][sid]
    #                 endnode=node
    #     sg.graph['endnodes'].append(endnode)

    if args.gml:
        write_gml(sg,"",outputfile=args.outfile)
    else:
        write_gfa(sg,"",outputfile=args.outfile, remap=False)

def graph_to_tree(G,cds):
    tree=IntervalTree()
    cds=G.graph['path2id'][cds]
    for node in G:
        if type(node)!=str:
            if cds in G.node[node]['offsets']:
                tree[G.node[node]['offsets'][cds]:G.node[node]['offsets'][cds]+len(G.node[node]['seq'])]=node
    return tree