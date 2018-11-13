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

    read_gfa(gfafile,None,"",G)
    
    for i,sg in enumerate(nx.weakly_connected_component_subgraphs(G)):
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
        
        sgsamples.sort()
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

        name="_".join([p for p in sg.graph['paths'] if not p.startswith("*")]).replace("|","_").replace(" ","_")[:200]
        
        logging.info("Write component (%d, size=%d) to: %s"%(i,len(sg.nodes()),name))
        write_gfa(sg,None,outputfile="%s.gfa"%name)