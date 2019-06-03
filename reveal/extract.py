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

#TODO: contribute patch back to networkx
def dag_longest_path_custom(G, weight='weight', default_weight=1):
    if not G:
        return []
    dist = {}  # stores {v : (length, u)}
    for v in nx.topological_sort(G):

        if type(G)==nx.MultiDiGraph:
            us = [(dist[u][0] + max([data[k].get(weight, default_weight) for k in data]), u)
                  for u, data in G.pred[v].items()]
        else:
            us = [(dist[u][0] + data.get(weight, default_weight), u)
                  for u, data in G.pred[v].items()]

        # Use the best predecessor if there is one and its distance is
        # non-negative, otherwise terminate.
        maxu = max(us, key=lambda x: x[0]) if us else (0, v)
        dist[v] = maxu if maxu[0] >= 0 else (0, v)
    u = None
    v = max(dist, key=lambda x: dist[x][0])
    path = []
    while u != v:
        path.append(v)
        u = v
        v = dist[v][1]
    path.reverse()
    return path


def extract(G,sample):
    logging.info("Extracting path: %s from graph (%s) of size: (%d,%d)"%(sample,type(G),G.number_of_nodes(),G.number_of_edges()))
    
    if sample == "_longest_":
        #shortcut to extract the "longest" path in terms of sequence

        if type(G)==nx.MultiDiGraph:
            sv=utils.MultiGraphToDiGraph(G)
            for v,t,k in G.edges:
                G[v][t][k]['weight']=len(G.node[t]['seq'])-G.node[t]['seq'].count("N") if 'seq' in G.node[t] else 0
        else:
            for v,t in G.edges:
                G[v][t]['weight']=len(G.node[t]['seq'])-G.node[t]['seq'].count("N") if 'seq' in G.node[t] else 0

        # p=[]
        seq=""
        # e=None
        # weights=[0]
        for n in dag_longest_path_custom(G, weight='weight'):
            # p.append(n)
            # if e!=None:
            #     if 0 in G[e][n]:
            #         weights.append(G[e][n][0]['weight'])
            #     else:
            #         weights.append(G[e][n]['weight'])
            seq+=G.node[n]['seq']
            # e=n

        # with open("path.txt",'w') as f:
        #     f.write("total length: %d\n"%sum(weights))
        #     for n,w in zip(p,weights):
        #         f.write("%s-%d\n"%(n,w))

        return seq
        
    elif sample not in G.graph['path2id']:
        logging.fatal("Unknown path: %s, graph contains: %s"%(sample, G.graph['path2id'].keys()))
        sys.exit(1)

    else:
        sid=G.graph['path2id'][sample]
        
        sg=[]
        for n1,n2,d in G.edges(data=True):
            if sid in d['paths']:
                sg.append((n1,n2,d))
        
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

def extract_path(G,path,type="str"):
    logging.debug("Extracting path of length: %d"%len(path))

    seq=""
    for n in path:
        nid,o=int(n[:-1]),n[-1:]
        assert(o=='+' or o=='-')

        if o=="+":
            seq+=G.node[nid]['seq']
        else:
            seq+=utils.rc(G.node[nid]['seq'])

    # return seq