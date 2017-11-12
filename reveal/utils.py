import networkx as nx
import logging
from intervaltree import Interval, IntervalTree
import sys
import os

def fasta_reader(fn,truncN=False,toupper=True,cutN=0):
    seq=""
    gapseq=""
    sub=0
    ntract=0
    with open(fn,'r') as ff:
        for line in ff:
            if line.startswith(">"):
                if seq!="":
                    yield name,seq
                    sub=0
                name=line.rstrip().replace(">","").replace("\t","")
                seq=""
            else:
                if truncN:
                    for base in line.rstrip():
                        if base.upper()=='N':
                            if len(seq)==0:
                                continue
                            elif seq[-1]=='N':
                                continue
                            else:
                                seq+='N'
                        else:
                            if toupper:
                                seq+=base.upper()
                            else:
                                seq+=base
                elif cutN>0:
                    for base in line.rstrip():
                        if base.upper()=='N':
                            gapseq+='N'
                        else:
                            if len(gapseq)<cutN:
                                seq+=gapseq
                                gapseq=""
                            else:
                                if seq!="":
                                    yield name+"_"+str(sub),seq
                                    seq=""
                                    sub+=1
                                gapseq=""
                            if toupper:
                                seq+=base.upper()
                            else:
                                seq+=base
                else:
                    if toupper:
                        seq+=line.upper().rstrip()
                    else:
                        seq+=line.rstrip()
        if seq!="":
            if cutN>0:
                yield name+"_"+str(sub),seq
            else:
                yield name,seq

def fasta_writer(fn,name_seq,lw=100):
    seq=""
    with open(fn,'w') as ff:
        for name,seq in name_seq:
            if not name.startswith(">"):
                name=">"+name+"\n"
            ff.write(name)
            for i in range( (len(seq)/lw)+(len(seq) % lw > 0)):
                ff.write(seq[i*lw:(i+1)*lw]+"\n")

def gapcost(pointa,pointb,model="sumofpairs"): #model is either sumofpairs or star
    assert(len(pointa)==len(pointb))
    
    if model=="star-avg":
        return abs(sum([pointa[i]-pointb[i] for i in range(len(pointa))]))/len(pointa)
    elif model=="star-med":
        return sorted([abs(pointa[i]-pointb[i]) for i in range(len(pointa))])[len(pointa)/2]
    elif model=="sumofpairs":
        p=0
        if len(pointa)==1:
            return abs(pointa[0]-pointb[0])
        D=[pointa[i]-pointb[i] for i in range(len(pointa))]
        for i in range(len(D)): #all pairwise distances
            for j in range(i+1,len(D)):
                p+=abs(D[i]-D[j])
        return p
    else:
        logging.warn("Uknown penalty model: %s."%model)
        return 0

def rc(seq):
    d = {'A':'T','C':'G','G':'C','T':'A','N':'N','a':'t','c':'g',\
        'g':'c','t':'a','n':'n','Y':'R','R':'Y','K':'M','M':'K',\
        'S':'S','W':'W','B':'V','V':'B','D':'H','H':'D','N':'N',\
        'X':'X','-':'-'}
    return "".join([d[b] for b in reversed(seq)])

def plotgraph(G, s1, s2, interactive=False, region=None, minlength=1):
    try:
        from matplotlib import pyplot as plt
        from matplotlib import patches as patches
    except:
        logging.error("Install matplotlib to generate mumplot.")
        return

    plt.xlabel(s1)
    plt.ylabel(s2)
    plt.title("REVEAL "+" ".join(sys.argv[1:]))
    maxx=0
    maxy=0
    
    minx=None
    miny=None

    #map names to ids
    s1=G.graph['sample2id'][s1]
    s2=G.graph['sample2id'][s2]

    for node,data in G.nodes(data=True):
        if 'seq' in data:
            l=len(data['seq']) #either seq argument
        else:
            l=node.end-node.begin #or interval as node
        if l<minlength:
            continue
        
        s1t=False
        s2t=False
        if s1 in data['offsets']:
            s1t=True
            if minx==None:
                minx=data['offsets'][s1]
            
            if data['offsets'][s1]+l > maxx:
                maxx=data['offsets'][s1]+l
            if data['offsets'][s1] < minx:
                minx=data['offsets'][s1]
        
        if s2 in data['offsets']:
            s2t=True
            if miny==None:
                miny=data['offsets'][s2]
            
            if data['offsets'][s2]+l > maxy:
                maxy=data['offsets'][s2]+l
            if data['offsets'][s2] < miny:
                miny=data['offsets'][s2]
        
        if s1t and s2t:
            plt.plot([data['offsets'][s1], data['offsets'][s1]+l], [data['offsets'][s2], data['offsets'][s2]+l], 'r-')
    
    if minx==None:
        minx=0
    if miny==None:
        miny=0
    
    plt.plot(minx,miny,'bx')
    plt.plot(maxx,maxy,'bx')
    
    if region!=None:
        rstart,rend=region.split(":")
        plt.axvline(x=int(rstart),linewidth=3,color='b',linestyle='dashed')
        plt.axvline(x=int(rend),linewidth=3,color='b',linestyle='dashed')
    
    if interactive:
        plt.show()
    else:
        plt.savefig("%s_%s.png"%(s1[:s1.rfind('.')],s2[:s2.rfind('.')]))


def read_gfa(gfafile, index, tree, graph, minsamples=1, maxsamples=None, targetsample=None, revcomp=False, remap=True):

    f=open(gfafile,'r')
    sep=";"
    nmapping={} #temp mapping object for nodeids in gfa file
    edges=[] #tmp list for edges
    paths=[]
    
    i=0
    gnodeid=graph.number_of_nodes()+1
    
    if 'samples' not in graph.graph:
        graph.graph['samples']=list()
    
    if 'id2sample' not in graph.graph:
        graph.graph['id2sample']=dict()
    
    if 'sample2id' not in graph.graph:
        graph.graph['sample2id']=dict()
    else:
        i=len(graph.graph['sample2id'])
        assert(i not in graph.graph['id2sample'])
    
    if 'id2end' not in graph.graph:
        graph.graph['id2end']=dict()

    for line in f:
    
        if line.startswith('H'):
            pass
        
        elif line.startswith('S'):
            s=line.strip().split('\t')
            nodeid=int(s[1])

            if remap:
                if graph.has_node(gnodeid):
                    logging.fatal("Id space for nodes is larger than total number of nodes in the graph.")
                    sys.exit(1)
            
            if index!=None:
                if revcomp:
                    intv=index.addsequence(rc(s[2]).upper())
                else:
                    intv=index.addsequence(s[2].upper())
                intv=Interval(intv[0],intv[1])
                tree.add(intv)
                graph.add_node(intv,aligned=0,offsets={})
                nmapping[nodeid]=intv
            else:
                if remap:
                    nmapping[nodeid]=gnodeid
                    gnodeid+=1
                else:
                    nmapping[nodeid]=nodeid

                if revcomp:
                    graph.add_node(nmapping[nodeid],seq=rc(s[2].upper()),aligned=0,offsets={})
                else:
                    graph.add_node(nmapping[nodeid],seq=s[2].upper(),aligned=0,offsets={})        
        
        elif line.startswith('L'):
            edges.append(line)
        
        elif line.startswith('P'): #traverse paths to add offset values
            paths.append(line)
    
    for line in edges:
        e=line.strip().split("\t")
        #assert(not graph.has_edge(nmapping[int(e[1])],nmapping[int(e[3])]))
        #assert(not graph.has_edge(nmapping[int(e[3])],nmapping[int(e[1])]))
        tags=dict()
        tags['ofrom']=e[2]
        tags['oto']=e[4]
        tags['cigar']=e[5]
        if '*' in e: #there are additional tags parse them
            for tag in e[7:]:
                key,ttype,value=tag.split(':')
                tags[key.lower()]=value
        tags['paths']=set()
        graph.add_edge(nmapping[int(e[1])],nmapping[int(e[3])],**tags)
    
    if len(paths)==0:
        logging.fatal("No paths defined in GFA, exit.")
        sys.exit(1)

    for line in paths:
        cols=line.rstrip().split("\t")
        sample=cols[1]
        
        if type(graph)==nx.DiGraph:
            logging.debug("DiGraph as input, so exclude *-paths.")
            if sample.startswith("*"):
                continue

        if sample in graph.graph['samples']:
            logging.fatal("ERROR: Graph already contains path for: %s"%sample)
            sys.exit(1)
        
        graph.graph['samples'].append(sample)
        
        if sample in graph.graph['sample2id']:
            logging.fatal("ERROR: Graph already contains path for: %s"%sample)
            sys.exit(1)
        
        sid=len(graph.graph['sample2id'])
        
        if sid in graph.graph['id2sample']:
            logging.fatal("ERROR: Id %d already linked to a path in the graph."%sid)
            sys.exit(1)
        
        graph.graph['sample2id'][sample]=sid
        graph.graph['id2sample'][sid]=sample

        o=0
        path=[(nid[:-1],nid[-1:]) for nid in cols[2].split(',')]

        for pi,gfn in enumerate(path):
            nid,orientation=gfn
            node=nmapping[int(nid)]

            # graph.graph['sample2path'][sample].append((node,orientation))

            graph.node[node]['offsets'][sid]=o

            if 'seq' in graph.node[node]:
                o+=len(graph.node[node]['seq'])
            elif isinstance(node,Interval):
                o+=node[1]-node[0]
            else:
                logging.warn("Node %s has unknown sequence content."%node)
            
            if pi==0:
                pnode=node
                pnid=nid
                porientation=orientation
                continue
            else:
                if node not in graph[pnode]:
                    logging.fatal("Path %s has %s -> %s, but no edge between these nodes exists in the graph definition!"%(sample,pnid,nid))
                assert(node in graph[pnode])
                if type(graph)==nx.MultiDiGraph:
                    for i in graph[pnode][node]:
                        if graph[pnode][node][i]['oto']==orientation and graph[pnode][node][i]['ofrom']==porientation:
                            graph[pnode][node][i]['paths'].add(sid)
                            break
                    else:
                        logging.fatal("Edge missing for path %s between %s (%s) and %s (%s)"%(sample,pnode,porientation,node,orientation))
                        sys.exit(1)
                else:
                    graph[pnode][node]['paths'].add(sid)

            pnode=node
            pnid=nid
            porientation=orientation
        
        graph.graph['id2end'][sid]=o
    
    #remove nodes and edges that are not associated to any path
    remove=[]
    for n1,n2,d in graph.edges_iter(data=True):
        if d['paths']==set(): #edge that is not traversed by any of the paths
            remove.append((n1,n2))
    graph.remove_edges_from(remove)

    remove=[]
    for n,d in graph.nodes_iter(data=True):
        if graph.node[n]['offsets']=={}: #edge that is not traversed by any of the paths
            remove.append(n)
    graph.remove_nodes_from(remove)

    if revcomp:
        genome2length=dict()
        #relabel the offsets, determine the length of all genomes in the graph, then l-pos
        for sample in graph.graph['samples']:
            maxp=0
            for node,data in graph.nodes(data=True):
                if graph.graph['sample2id'][sample] in data['offsets']:
                    if data['offsets'][graph.graph['sample2id'][sample]]+ (node[1]-node[0]) >maxp:
                        maxp=data['offsets'][graph.graph['sample2id'][sample]]+(node[1]-node[0])
            genome2length[sample]=maxp
        
        for sample in graph.graph['samples']:
            for node,data in graph.nodes(data=True):
                if graph.graph['sample2id'][sample] in data['offsets']:
                    graph.node[node]['offsets'][graph.graph['sample2id'][sample]]=genome2length[sample]-(graph.node[node]['offsets'][graph.graph['sample2id'][sample]]+(node[1]-node[0]))
        
        graph.reverse(copy=False)

#simply write sequence without the graph topology
def write_fasta(G,T,outputfile="reference.fasta"):
    f=open(outputfile,'wb')
    for i,node in enumerate(nx.topological_sort(G)):
        if isinstance(node,str):
            if node=='start' or node=='end':
                continue
        i+=1
        data=G.node[node]
        seq=""
        if len(node)==3:
            nodename=node[2]
        else:
            nodename=str(node)
        if 'seq' in data:
            f.write(">%s\n"%nodename)
            f.write(data['seq'].upper()+"\n")
        else:
            if isinstance(node,Interval):
                f.write(">%s\n"%nodename)
                f.write(T[node.begin:node.end].upper()+"\n")
            else:
                f.write(">%s\n"%nodename)
                logging.warn("No sequence for node: %s"%nodename)
    f.close()

def write_gfa(G,T,outputfile="reference.gfa",nometa=False, paths=True, remap=True):
    
    if not outputfile.endswith(".gfa"):
        outputfile+=".gfa"
    
    f=open(outputfile,'wb')
    sep=';'
    f.write('H\tVN:Z:1.0\tCL:Z:%s\n'%" ".join(sys.argv))
    
    sample2id=dict()
    sample2id=G.graph['sample2id']

    mapping={}
    
    if type(G)==nx.DiGraph:
        iterator=nx.topological_sort(G)
        logging.debug("Writing gfa in topological order (%d)."%len(iterator))
    elif type(G)==nx.MultiDiGraph:
        iterator=G.nodes()
        logging.debug("Writing gfa in random order (%d)."%len(iterator))
    else:
        logging.fatal("Unsupported graph type: %s"%type(G))
        sys.exit(1)

    for i,node in enumerate(iterator): #iterate once to get a mapping of ids to intervals
        mapping[node]=i+1 #one-based for vg

    for i,node in enumerate(iterator):
        i+=1
        data=G.node[node]
        seq=""

        if 'seq' in data:
            f.write('S\t'+str(mapping[node])+'\t'+data['seq'].upper())
            seq=data['seq']
        else:
            if isinstance(node,Interval):
                seq=T[node.begin:node.end].upper()
                f.write('S\t'+str(mapping[node])+'\t'+seq)
            elif isinstance(node,tuple):
                seq=T[node[0]:node[0]+G.node[node]['l']].upper()
                f.write('S\t'+str(mapping[node])+'\t'+seq)
            else:
                f.write('S\t'+str(mapping[node])+'\t')
        
        f.write("\n")
        
        for to in G[node]:
            if type(G)==nx.MultiDiGraph:
                for edgeid in G[node][to]:
                    f.write("L\t"+str(mapping[node])+"\t"+G[node][to][edgeid]['ofrom']+"\t"+str(mapping[to])+"\t"+G[node][to][edgeid]['oto']+"\t0M\n")
            else:
                if 'ofrom' in G[node][to] and 'oto' in G[node][to]:
                    f.write("L\t"+str(mapping[node])+"\t"+G[node][to]['ofrom']+"\t"+str(mapping[to])+"\t"+G[node][to]['oto']+"\t0M\n")
                else: #if not there, assume same orientation
                    f.write("L\t"+str(mapping[node])+"\t+\t"+str(mapping[to])+"\t+\t0M\n")
    
    #write paths
    for sample in G.graph['samples']:
        logging.debug("Writing path: %s"%sample)

        sid=G.graph['sample2id'][sample]
        subgraph=[]
        for e1,e2,d in G.edges(data=True):
            if sid in d['paths']:
                subgraph.append((e1,e2,d))

        path=[]
        if len(subgraph)>0:
            sg=nx.DiGraph(subgraph)

            if (len([c for c in nx.connected_components(sg.to_undirected())] )!=1):
                write_gml(sg,None,outputfile="%s_subgraph.gml"%sample)

            nodepath=nx.topological_sort(sg)
            pn=nodepath[0]
            for n in nodepath[1:]:
                #assert(n in G[pn]) #path is unconnected in graph! Something went wrong..

                if n not in G[pn]:
                    logging.error("Path %s spells a path that is not supported by the graph %s->%s!"%(sample,str(mapping[n]),str(mapping[pn])))
                    write_gml(sg,None,outputfile="%s_subgraph.gml"%sample)
                    break
                else:
                    path.append("%d%s"% (mapping[pn], sg[pn][n]['ofrom'] if 'ofrom' in sg[pn][n] else '+') )
                
                if n==nodepath[-1]:
                    path.append("%d%s"% (mapping[n], sg[pn][n]['oto'] if 'oto' in sg[pn][n] else '+') )
                pn=n
        else: #Maybe just a single node, strictly not a path..
            for node,data in G.nodes_iter(data=True):
                if sid in data['offsets']:
                    path.append("%s+"%mapping[node]) #TODO: have no way of knowing what the original orientation was now...
            if len(path)==0: #not just a single node, sample not part of the graph, dont write path
                continue

        f.write("P\t"+sample+"\t"+",".join(path)+"\t"+",".join(["0M"]*len(path))+"\n")
    
    f.close()


def write_gml(G,T,outputfile="reference",partition=False,hwm=4000):
    G=G.copy()
    mapping={}

    if 'samples' in G.graph:
        totn=len(G.graph['samples'])
        logging.debug("Graph contains %d samples"%totn)
    else:
        totn=0
    
    for key in G.graph:
        G.graph[key]=str(G.graph[key])
    
    if type(G)==nx.MultiDiGraph:
        for n1,n2,k,d in G.edges(keys=True,data=True):
            for key in d:
                v=d[key]
                if not isinstance(v,str) and not isinstance(v,int):
                    G[n1][n2][k][key]=str(v)
    else:
        for n1,n2,d in G.edges(data=True):
            for key in d:
                v=d[key]
                if not isinstance(v,str) and not isinstance(v,int):
                    G[n1][n2][key]=str(v)

    for n,d in G.nodes(data=True):
        mapping[n]=str(n)
        d=G.node[n]
        
        if 'offsets' in d:
            G.node[n]['n']=len(d['offsets'])
        
        for key in d:
            v=d[key]
            if not isinstance(v,str) and not isinstance(v,int):
                G.node[n][key]=str(v)
        
        if 'seq' not in G.node[n]:
            if isinstance(n,Interval):
                G.node[n]['seq']=T[n.begin:n.end].upper()
            else:
                G.node[n]['seq']=""
        G.node[n]['l']=len(G.node[n]['seq'])
    
    G=nx.relabel_nodes(G,mapping)

    outputfiles=[]
    
    if partition:
        logging.debug("Trying to partion graph into subgraphs of size %d."%hwm)
        i=0
        for sgi,subset in enumerate(nx.connected_components(G.to_undirected())):
            logging.debug("Partitioning connected component: %d"%sgi)
            sgn=[]
            g=G.subgraph(subset)
            gn=G.graph['samples']
            for n in nx.topological_sort(g):
                sgn.append(n)
                if G.node[n]['n']==totn: #join/split node
                    logging.debug("Can split graph at node: %s."%n)
                    if len(sgn)>=hwm:
                        logging.debug("Splitting graph at node: %s"%n)
                        sg=G.subgraph(sgn)
                        fn=outputfile+'.'+str(i)+'.gml'
                        nx.write_gml(sg,fn)
                        outputfiles.append(fn)
                        sgn=[n]
                        i+=1
            
            if len(sgn)>1:
                sg=G.subgraph(sgn)
                fn=outputfile+'.'+str(i)+'.gml'
                nx.write_gml(sg,fn)
                i+=1
                outputfiles.append(fn)
    else:
        if not outputfile.endswith(".gml"):
            outputfile=outputfile+'.gml'
        nx.write_gml(G,outputfile)
        outputfiles.append(outputfile)
    
    return outputfiles

def kdtree(points, k, depth=0):
    n=len(points)    
    if n==0:
        return None
    if n==1:
        return points[0]
    splitdim=depth % k
    spoints=sorted(points,key=lambda p: p[splitdim])
    splitvalue=spoints[n/2][splitdim] #take median for splitting
    if splitvalue==spoints[0][splitdim]:
        splitvalue+=1
    left=[p for p in spoints if p[splitdim] < splitvalue]
    right=[p for p in spoints if p[splitdim] >= splitvalue]
    return { 'left': kdtree(left, k, depth=depth+1) , 'split' : splitvalue, 'right': kdtree(right, k, depth=depth+1) }

#return all points within the query range
def range_search(kdtree, qstart, qend):
    k=len(qstart)
    points=[]
    stack=[(kdtree,0)]
    while len(stack)!=0:
        tree,depth=stack.pop()
        splitdim=depth%k

        if tree==None: #something is wrong
            continue

        if isinstance(tree,tuple): #reached leaf, tree==point
            if tree[splitdim]>=qstart[splitdim] and tree[splitdim]<=qend[splitdim]:
                #check ik point is contained in the range query
                for d in range(k):
                    if tree[d]>=qstart[d] and tree[d]<=qend[d]:
                        continue
                    else:
                        break
                else:
                    points.append(tree)
            continue
        
        splitvalue=tree['split']

        if splitvalue>=qstart[splitdim] and splitvalue<=qend[splitdim]: #intersect
            if splitvalue != qstart[splitdim]: #equal values go right
                stack.append((tree['left'],depth+1))
            stack.append((tree['right'],depth+1))
        elif splitvalue < qstart[splitdim]:
            stack.append((tree['right'],depth+1))
        else:
            stack.append((tree['left'],depth+1))

    return points
