import networkx as nx
from utils import *
import sys

def split_cmd(args):
    if len(args.gfa)!=1:
        logging.fatal("Specify 1 gfa file.")
        return
    
    if not args.gfa[0].endswith('.gfa'):
        logging.fatal("Use .gfa as extension of the gfa file.")
        return
    
    if args.nocycles:
        G=nx.DiGraph()
    else:
        G=nx.MultiDiGraph()

    split(G,args.gfa[0])

def split(G,gfafile):

    logging.debug("Reading graph...")
    read_gfa(gfafile,None,"",G)
    logging.debug("Done.")

    for i,sg in enumerate(nx.weakly_connected_component_subgraphs(G)):
        
        sgpaths=[]

        sids=set()
        for node in sg.nodes():
            if type(node)!=str:
                for sid in sg.node[node]['offsets']:
                    sids.add(sid)
        
        for sid in sids:
            sgpaths.append(sg.graph['id2path'][sid])

        sg.graph['paths']=sgpaths
        sg.graph['id2path']={sid:sg.graph['id2path'][sid] for sid in sids}
        sg.graph['path2id']={path:sg.graph['path2id'][path] for path in sgpaths}

        name="_".join([p for p in sorted(sgpaths) if not p.startswith("*")]).replace("|","_").replace(" ","_")[:200]

        logging.info("Write component (%d, size=%d) to: %s"%(i,len(sg.nodes()),name))
        write_gfa(sg,None,outputfile="%s.gfa"%name)