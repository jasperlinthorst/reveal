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
    
    split(args.gfa[0])

def split(gfafile):
    G=nx.DiGraph()
    read_gfa(gfafile,None,"",G)
    
    for i,sub in enumerate(nx.connected_components(G.to_undirected())):
        sg=nx.subgraph(G,sub)
        sids=set()
        for node in sg.nodes():
            for sid in sg.node[node]['offsets']:
                sids.add(sid)
        
        #determine a mapping
        mapping=dict()
        sgsamples=[]
        for sid in sids:
            mapping[sid]=len(sgsamples)
            sg.graph['path2id'][G.graph['id2path'][sid]]=len(sgsamples)
            sgsamples.append(G.graph['id2path'][sid])
        
        sg.graph['paths']=sgsamples
        
        for e1,e2,d in sg.edges(data=True):
            np=set()
            for p in d['paths']:
                np.add(mapping[p])
            d['paths']=np
        
        for n,d in sg.nodes(data=True):
            no=dict()
            for p in d['offsets']:
                no[mapping[p]]=d['offsets'][p]
            sg.node[n]['offsets']=no
        name="_".join(sg.graph['paths'])
        write_gfa(sg,None,outputfile="%s.gfa"%name)
