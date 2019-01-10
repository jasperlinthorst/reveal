import networkx as nx
import utils
import sys
import logging
import re

def extract_cmd(args):
    if not args.graph[0].endswith(".gfa"):
        logging.fatal("Invalid gfa file.")
        return
    width=args.width

    if args.nocycles:
        G=nx.DiGraph()
    else:
        G=nx.MultiDiGraph()

    utils.read_gfa(args.graph[0], None, None, G, remap=False)

    if args.all:
        args.input=sorted(G.graph['paths'])
    
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
    logging.debug("Extracting path: %s from graph (%s) of size: (%d,%d)"%(sample,type(G),G.number_of_nodes(),G.number_of_edges()))

    if sample not in G.graph['path2id']:
        logging.fatal("Unknown path: %s, graph contains: %s"%(sample, G.graph['path2id'].keys()))
        sys.exit(1)

    sid=G.graph['path2id'][sample]
    
    sg=[]
    for n1,n2,d in G.edges(data=True):
        if sid in d['paths']:
            sg.append((n1,n2,d))

    # if type(G)==nx.MultiDiGraph:
    #     sg=[(n1,n2,G[n1][n2][i]) for n1,n2,i in G.edges if sid in G[n1][n2][i]['paths']]
    # else:
    #     sg=[(n1,n2) for n1,n2 in G.edges if sid in G[n1][n2]['paths']]
        
    if len(sg)>0:
        #G can be a MultiDiGraph, but subgraph should be single edge!
        subgraph=nx.DiGraph(sg)
        seq=""
        path=list(nx.topological_sort(subgraph))

        if type(G)==nx.MultiDiGraph:
            inito=G[path[0]][path[1]][0]['ofrom']
        else:
            inito=G[path[0]][path[1]]['ofrom']

        pnode=None

        for node in path:
            offset=0
            if pnode==None:
                o=inito
            else:
                o=subgraph[pnode][node]['oto']
                if 'cigar' in subgraph[pnode][node] and subgraph[pnode][node]['cigar']!='0M':
                    cigar=subgraph[pnode][node]['cigar']
                    a=re.findall(r'(\d+)(\w)', cigar)
                    for l,t in a: #determine offset within the segment to allow for overlapping segments
                        if t=='M' or t=='I' or t=='S' or t=='P': #source of the edge (pnode) is considered the reference
                            offset+=int(l)
                
            if o=="+":
                s=G.node[node]['seq']
            else:
                s=utils.rc(G.node[node]['seq'])

            assert(len(s)>=offset)

            seq+=s[offset:]
            pnode=node

    else: #has to be a single node
        seq=""
        for n in G:
            if sid in G.node[n]['offsets']:
                seq=G.node[n]['seq']
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

    # return seq