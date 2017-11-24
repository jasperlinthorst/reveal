import networkx as nx
import utils

def comp(G):
    for node in G.node:
        G.node[node]['seq']=utils.rc(G.node[node]['seq'])
    
    genome2length=dict()
    #relabel the offsets, determine the length of all genomes in the graph, then l-pos
    for sample in G.graph['paths']:
        maxp=0
        for node,data in G.nodes(data=True):
            if sample in data['offsets']:
                if data['offsets'][sample]+len(data['seq'])>maxp:
                    maxp=data['offsets'][sample]+len(data['seq'])
        genome2length[sample]=maxp
    
    for sample in G.graph['paths']:
        for node,data in G.nodes(data=True):
            if sample in data['offsets']:
                G.node[node]['offsets'][sample]=genome2length[sample]-(G.node[node]['offsets'][sample]+len(data['seq']))
    
    G.reverse(copy=False)
    return G

def comp_cmd(args):
    g=nx.DiGraph()
    g.graph['paths']=[]
    utils.read_gfa(args.graph[0],None,None,g,targetsample=None)
    g=comp(g)
    utils.write_gfa(g,"",outputfile=args.graph[0].replace('.gfa','.rc.gfa'), nometa=False)
