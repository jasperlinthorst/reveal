import networkx as nx
import logging

#OVERRIDE intervaltree to bug in incorrect __eq__ function
import intervaltree
class IntervalPatched(intervaltree.Interval):    
    def __eq__(self,other):
        if type(self)==type(other)==Interval:
            return super(intervaltree.Interval,self).__eq__(other)
        else:
            return False
intervaltree.Interval=IntervalPatched
from intervaltree import Interval, IntervalTree


import sys
import os
from math import log
import uuid
import subprocess

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

def gapcost(pointa,pointb,model="sumofpairs",convex=True): #model is either sumofpairs or star
    assert(len(pointa)==len(pointb))
    
    if model=="star-avg":
        return abs(sum([pointa[i]-pointb[i] for i in range(len(pointa))]))/len(pointa)
    elif model=="star-med":
        return sorted([abs(pointa[i]-pointb[i]) for i in range(len(pointa))])[len(pointa)/2]
    elif model=="sumofpairs":
        p=0
        D=[abs(pointa[i]-pointb[i]) for i in range(len(pointa))]
        for i in range(len(D)): #all pairwise distances
            for j in range(i+1,len(D)):
                if convex:
                    p+=log(abs(D[i]-D[j])+1)
                else:
                    p+=abs(D[i]-D[j])
        return p
    else:
        logging.warn("Unknown penalty model: %s."%model)
        return 0

def rc(seq):
    d = {'A':'T','C':'G','G':'C','T':'A','N':'N','a':'t','c':'g',\
        'g':'c','t':'a','n':'n','Y':'R','R':'Y','K':'M','M':'K',\
        'S':'S','W':'W','B':'V','V':'B','D':'H','H':'D','N':'N',\
        'X':'X','-':'-'}
    return "".join([d[b] for b in reversed(seq)])

#extract all combinations for 'non-unique max exact matches' and return list of mems as if they were unique
def mem2mums(mem):
    import itertools
    l,n,spd=mem
    spd=sorted(spd, key=lambda m:m[1])
    pos=[[spd[0]]]
    for i in range(1,len(spd)):
        if spd[i-1][0]==spd[i][0]:
            pos[-1].append(spd[i])
        else:
            pos.append([spd[i]])
    mums=[]
    for t in itertools.product(*pos):
        yield (l,n,t)

def plotgraph(G, s1, s2, interactive=False, region=None, minlength=1):
    try:
        from matplotlib import pyplot as plt
        from matplotlib import patches as patches
    except:
        logging.error("Install matplotlib to generate mumplot.")
        return

    logging.debug("Generating plot for %s and %s."%(s1,s2))
    plt.xlabel(s1)
    plt.ylabel(s2)
    plt.title("REVEAL "+" ".join(sys.argv[1:]))
    maxx=0
    maxy=0
    
    minx=None
    miny=None

    #map names to ids
    s1=G.graph['path2id'][s1]
    s2=G.graph['path2id'][s2]

    logging.debug("Samples in graph: %s"%G.graph['path2id'])
    logging.debug("Generating plot for %s and %s, with minlength=%d."%(s1,s2,minlength))

    for node,data in G.nodes(data=True):
        if type(node)==str:
            continue

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
        else:
            continue

        if s2 in data['offsets']:
            s2t=True
            if miny==None:
                miny=data['offsets'][s2]
            if data['offsets'][s2]+l > maxy:
                maxy=data['offsets'][s2]+l
            if data['offsets'][s2] < miny:
                miny=data['offsets'][s2]
        else:
            continue

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
        plt.savefig("%s_%s.png"%(s1,s2))


def read_fasta(fasta, index, tree, graph, contigs=True):
    logging.info("Reading fasta: %s ..." % fasta)
    
    if 'paths' not in graph.graph:
        graph.graph['paths']=list()
    
    if 'id2path' not in graph.graph:
        graph.graph['id2path']=dict()
    
    if 'path2id' not in graph.graph:
        graph.graph['path2id']=dict()
    
    if 'id2end' not in graph.graph:
        graph.graph['id2end']=dict()

    if contigs:
        index.addsample(os.path.basename(fasta))
        for name,seq in fasta_reader(fasta):
            sid=len(graph.graph['paths'])
            name=name.replace(":","").replace(";","")
            if name in graph.graph['paths']:
                logging.fatal("Fasta with this name: \"%s\" is already contained in the graph."%name)
                sys.exit(1)
            graph.graph['paths'].append(name)
            graph.graph['path2id'][name]=sid
            graph.graph['id2path'][sid]=name
            graph.graph['id2end'][sid]=len(seq)
            intv=index.addsequence(seq.upper())
            logging.debug("Adding interval: %s"%str(intv))
            Intv=Interval(intv[0],intv[1])
            tree.add(Intv)
            startnode=uuid.uuid4().hex
            endnode=uuid.uuid4().hex
            graph.add_node(startnode,offsets={sid:0},endpoint=True)
            graph.add_node(Intv,offsets={sid:0},aligned=0)
            graph.add_node(endnode,offsets={sid:len(seq)},endpoint=True)
            graph.add_edge(startnode,Intv,paths=set([sid]),ofrom='+',oto='+')
            graph.add_edge(Intv,endnode,paths=set([sid]),ofrom='+',oto='+')
    else: #treat every sequence in the multifasta as a target
        for name,seq in fasta_reader(fasta):
            index.addsample(name)
            sid=len(graph.graph['paths'])
            name=name.replace(":","").replace(";","")
            if name in graph.graph['paths']:
                logging.fatal("Fasta with this name: \"%s\" is already contained in the graph."%name)
                sys.exit(1)
            graph.graph['paths'].append(name)
            graph.graph['path2id'][name]=sid
            graph.graph['id2path'][sid]=name
            graph.graph['id2end'][sid]=len(seq)
            intv=index.addsequence(seq.upper())
            logging.debug("Adding interval: %s"%str(intv))
            Intv=Interval(intv[0],intv[1])
            tree.add(Intv)
            startnode=uuid.uuid4().hex
            endnode=uuid.uuid4().hex
            graph.add_node(startnode,offsets={sid:0},endpoint=True)
            graph.add_node(Intv,offsets={sid:0},aligned=0)
            graph.add_node(endnode,offsets={sid:len(seq)},endpoint=True)
            graph.add_edge(startnode,Intv,paths=set([sid]),ofrom='+',oto='+')
            graph.add_edge(Intv,endnode,paths=set([sid]),ofrom='+',oto='+')

def read_gfa(gfafile, index, tree, graph, minsamples=1, maxsamples=None, targetsample=None, revcomp=False, remap=True):
    f=open(gfafile,'r')
    sep=";"
    nmapping={} #temp mapping object for nodeids in gfa file
    edges=[] #tmp list for edges
    paths=[]
    
    i=0
    gnodeid=graph.number_of_nodes()+1
    
    if 'paths' not in graph.graph:
        graph.graph['paths']=list()
    
    if 'id2path' not in graph.graph:
        graph.graph['id2path']=dict()
    
    if 'path2id' not in graph.graph:
        graph.graph['path2id']=dict()
    else:
        i=len(graph.graph['path2id'])
        assert(i not in graph.graph['id2path'])
    
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

        if type(graph)==nx.DiGraph and (e[2]!='+' or e[4]!='+'):
            continue #skip these edges if we only want a directed acyclic graph

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

    if index==None:
        graph.graph['noffset']=max([v for v in nmapping.values() if type(v)==int])+1

    startnodes=set()
    endnodes=set()

    for line in paths:
        cols=line.rstrip().split("\t")
        sample=cols[1]
        
        if type(graph)==nx.DiGraph:
            if sample.startswith("*"):
                logging.debug("DiGraph as input, so exclude path: %s"%sample)
                continue

        if sample in graph.graph['paths']:
            logging.fatal("ERROR: Graph already contains path for: %s"%sample)
            sys.exit(1)
        
        graph.graph['paths'].append(sample)
        
        if sample in graph.graph['path2id']:
            logging.fatal("ERROR: Graph already contains path for: %s"%sample)
            sys.exit(1)
        
        sid=len(graph.graph['path2id'])
        
        if sid in graph.graph['id2path']:
            logging.fatal("ERROR: Id %d already linked to a path in the graph."%sid)
            sys.exit(1)
        
        graph.graph['path2id'][sample]=sid
        graph.graph['id2path'][sid]=sample

        o=0
        
        path=[(nid[:-1],nid[-1:]) for nid in cols[2].split(',')]
        
        for pi,gfn in enumerate(path):
            nid,orientation=gfn
            node=nmapping[int(nid)]
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
        
        start=uuid.uuid4().hex
        graph.add_node(start,offsets={sid:0},endpoint=True)
        graph.add_edge(start,nmapping[int(path[0][0])],paths={sid},ofrom='+',oto=path[0][1])
        startnodes.add(start)

        end=uuid.uuid4().hex
        graph.add_node(end,offsets={sid:o},endpoint=True)
        graph.add_edge(nmapping[int(path[-1][0])],end,paths={sid},ofrom=path[-1][1],oto='+')
        endnodes.add(end)
        
        graph.graph['id2end'][sid]=o
    
    #remove nodes and edges that are not associated to any path
    remove=[]
    for n1,n2,d in graph.edges(data=True):
        if d['paths']==set(): #edge that is not traversed by any of the paths
            remove.append((n1,n2))
    if len(remove)>0:
        logging.debug("Removing %d edges from the graph as they are not traversed..."%len(remove))
        graph.remove_edges_from(remove)
        logging.debug("Done.")

    remove=[]
    for n,d in graph.nodes(data=True):
        if graph.node[n]['offsets']=={}: #node that is not traversed by any of the paths
            remove.append(n)
    if len(remove)>0:
        logging.debug("Removing %d nodes from the graph as they are not traversed..."%len(remove))
        graph.remove_nodes_from(remove)
        logging.debug("Done.")

    logging.debug("Converting to undirected graph...")
    Gu=graph.to_undirected()
    logging.debug("Done.")

    logging.debug("Extracting connected components...")
    conncomp=nx.connected_components(Gu)
    logging.debug("Done.")

    #merge start/end nodes per connected component in the graph
    for i,comp in enumerate(conncomp):
        logging.debug("Inspecting connected component: %d (%d)"%(i,len(comp)))

        startmerge=set()
        endmerge=set()
        for node in comp:
            if node in startnodes:
                startmerge.add(node)
            if node in endnodes:
                endmerge.add(node)

        if len(endmerge)>0:
            endnode=uuid.uuid4().hex            
            graph.add_node(endnode,offsets={},seq="",endpoint=True) #add dummy node for end of each sequence in the subgraph
            
            for node in endmerge: #copy offset values
                for k in graph.node[node]['offsets']:
                    graph.node[endnode]['offsets'][k]=graph.node[node]['offsets'][k]
                #reconnect
                pred=set()
                predids=set()
                for pnode in graph.predecessors(node):
                    if type(graph)==nx.MultiDiGraph:
                        graph.add_edge(pnode,endnode,paths=graph[pnode][node][0]['paths'],ofrom=graph[pnode][node][0]['ofrom'],oto=graph[pnode][node][0]['oto'])
                    else:
                        if not graph.has_edge(pnode,endnode):
                            graph.add_edge(pnode,endnode,paths=graph[pnode][node]['paths'],ofrom=graph[pnode][node]['ofrom'],oto=graph[pnode][node]['oto'])
                        else:
                            for p in graph[pnode][node]['paths']:
                                graph[pnode][endnode]['paths'].add(p)
        
        if len(startmerge)>0:
            startnode=uuid.uuid4().hex
            graph.add_node(startnode,offsets={},seq="",endpoint=True) #add dummy node for start of each sequence in the subgraph

            for node in startmerge: #copy offset values
                for k in graph.node[node]['offsets']:
                    graph.node[startnode]['offsets'][k]=graph.node[node]['offsets'][k]
                #reconnect
                pred=set()
                predids=set()
                for nnode in graph.successors(node):
                    if type(graph)==nx.MultiDiGraph:
                        graph.add_edge(startnode,nnode,paths=graph[node][nnode][0]['paths'],ofrom=graph[node][nnode][0]['ofrom'],oto=graph[node][nnode][0]['oto'])
                    else:
                        if not graph.has_edge(startnode,nnode):
                            graph.add_edge(startnode,nnode,paths=graph[node][nnode]['paths'],ofrom=graph[node][nnode]['ofrom'],oto=graph[node][nnode]['oto'])
                        else:
                            for p in graph[node][nnode]['paths']:
                                graph[startnode][nnode]['paths'].add(p)

        graph.remove_nodes_from(list(startmerge)+list(endmerge))
        logging.debug("Done.")

    if revcomp:
        genome2length=dict()
        #relabel the offsets, determine the length of all genomes in the graph, then l-pos
        for sample in graph.graph['paths']:
            maxp=0
            for node,data in graph.nodes(data=True):
                if graph.graph['path2id'][sample] in data['offsets']:
                    if data['offsets'][graph.graph['path2id'][sample]]+ (node[1]-node[0]) >maxp:
                        maxp=data['offsets'][graph.graph['path2id'][sample]]+(node[1]-node[0])
            genome2length[sample]=maxp
        
        for sample in graph.graph['paths']:
            for node,data in graph.nodes(data=True):
                if graph.graph['path2id'][sample] in data['offsets']:
                    graph.node[node]['offsets'][graph.graph['path2id'][sample]]=genome2length[sample]-(graph.node[node]['offsets'][graph.graph['path2id'][sample]]+(node[1]-node[0]))
        
        graph.reverse(copy=False)

#simply write sequence without the graph topology
def write_fasta(G,T,outputfile="reference.fasta"):
    f=open(outputfile,'wb')
    for i,node in enumerate(nx.topological_sort(G)):
        if isinstance(node,str):
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
    sample2id=G.graph['path2id']

    mapping={}
    
    if type(G)==nx.DiGraph or type(G)==nx.classes.graphviews.SubDiGraph:
        iterator=nx.topological_sort(G)
        logging.debug("Writing gfa in topological order.")
    elif type(G)==nx.MultiDiGraph or type(G)==nx.classes.graphviews.SubMultiDiGraph:
        iterator=G.nodes()
        logging.debug("Writing gfa in random order.")
    else:
        logging.fatal("Unsupported graph type: %s"%type(G))
        sys.exit(1)

    iterator=[node for node in iterator if type(node)!=str] #exclude start/end node

    if remap:
        for i,node in enumerate(iterator): #iterate once to get a mapping of ids to intervals
            mapping[node]=i+1
    else:
        for node in iterator: #quick and dirty..
            mapping[node]=node

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
            if 'endpoint' in G.node[to]:
                continue
            if type(G)==nx.MultiDiGraph:
                for edgeid in G[node][to]:
                    if 'cigar' in G[node][to][edgeid]:
                        cigar=G[node][to][edgeid]['cigar']
                    f.write("L\t"+str(mapping[node])+"\t"+G[node][to][edgeid]['ofrom']+"\t"+str(mapping[to])+"\t"+G[node][to][edgeid]['oto']+"\t"+(G[node][to][edgeid]['cigar'] if 'cigar' in G[node][to][edgeid] else "0M")+"\n")
            else:
                if 'ofrom' in G[node][to] and 'oto' in G[node][to]:
                    f.write("L\t"+str(mapping[node])+"\t"+G[node][to]['ofrom']+"\t"+str(mapping[to])+"\t"+G[node][to]['oto']+"\t"+(G[node][to]['cigar'] if 'cigar' in G[node][to] else "0M")+"\n")
                else: #if not there, assume same orientation
                    f.write("L\t"+str(mapping[node])+"\t+\t"+str(mapping[to])+"\t+\t"+(G[node][to]['cigar'] if 'cigar' in G[node][to] else "0M")+"\n")
    
    #write paths
    for sample in G.graph['paths']:
        sid=G.graph['path2id'][sample]
        logging.debug("Writing path: %s (sid=%d)"%(sample,sid))

        subgraph=[]
        for e1,e2,d in G.edges(data=True):
            if sid in d['paths']:
                if type(e1)!=str and type(e2)!=str: #string nodes are dummy start/stop nodes, ignore
                    subgraph.append((e1,e2,d))

        path=[]
        cigarpath=[]
        if len(subgraph)>0:

            sg=nx.DiGraph(subgraph)

            for node in sg.nodes():
                if len(sg[node])>1:
                    print "writing gml",sample
                    write_gml(sg,None,outputfile="%s_subgraph.gml"%sample)
                    print "done"
                assert(len(sg[node])<2) #There should only one path!

            # if (len([c for c in nx.connected_components(sg.to_undirected())] )!=1):
            #     write_gml(sg,None,outputfile="%s_subgraph.gml"%sample)

            nodepath=list(nx.topological_sort(sg))

            if type(nodepath[0])==str and type(nodepath[-1])==str:
                nodepath=nodepath[1:-1]

            pn=nodepath[0]
            for n in nodepath[1:]:
                #assert(n in G[pn]) #path is unconnected in graph! Something went wrong..
                if n not in G[pn]:
                    logging.error("Path %s spells a path that is not supported by the graph %s->%s (%s->%s)!"%(sample,str(mapping[n]),str(mapping[pn]),n,pn))
                    write_gml(sg,None,outputfile="%s_subgraph.gml"%sample)
                    break
                else:
                    path.append("%d%s"% (mapping[pn], sg[pn][n]['ofrom'] if 'ofrom' in sg[pn][n] else '+') )
                    cigarpath.append("%s"% (sg[pn][n]['cigar'] if 'cigar' in sg[pn][n] else '0M') )
                
                if n==nodepath[-1]: #if last node
                    path.append("%d%s"% (mapping[n], sg[pn][n]['oto'] if 'oto' in sg[pn][n] else '+') )
                    cigarpath.append("%s"% (sg[pn][n]['cigar'] if 'cigar' in sg[pn][n] else '0M') )
                pn=n
        elif len(subgraph)==1: #Maybe just a single node, strictly not a path..
            for node,data in G.nodes(data=True):
                if sid in data['offsets'] and type(node)!=str:
                    path.append("%s+"%mapping[node]) #TODO: have no way of knowing what the original orientation was now...
        else:
            #if len(path)==0: #not just a single node, sample not part of the graph, dont write path
            logging.debug("Sample %s not part of the graph."%sid)
            continue

        f.write("P\t"+sample+"\t"+",".join(path)+"\t"+",".join(cigarpath)+"\n")
    
    f.close()


def write_gml(G,T,outputfile="reference",partition=False,hwm=4000):
    G=G.copy()
    mapping={}

    if 'paths' in G.graph:
        totn=len(G.graph['paths'])
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

            if type(v)!=str and type(v)!=int:
                G.node[n][key]=str(v)
        
        if 'seq' not in G.node[n]:
            if isinstance(n,Interval):
                G.node[n]['seq']=T[n.begin:n.end].upper()
            else:
                G.node[n]['seq']=""
        G.node[n]['l']=len(G.node[n]['seq'])
        # G.node[n]['seqstart']=G.node[n]['seq'][:20]
        G.node[n]['seqend']=G.node[n]['seq'][-20:]
    
    G=nx.relabel_nodes(G,mapping)

    outputfiles=[]
    
    if partition:
        logging.debug("Trying to partion graph into subgraphs of size %d."%hwm)
        i=0
        for sgi,subset in enumerate(nx.connected_components(G.to_undirected())):
            logging.debug("Partitioning connected component: %d"%sgi)
            sgn=[]
            g=G.subgraph(subset)
            gn=G.graph['paths']
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

#copy interval based nodes to sequence attribute
def seq2node(G,T,toupper=True,remap=False):
    i=1
    mapping=dict()
    for node in G:
        if isinstance(node,Interval):
            if toupper:
                G.node[node]['seq']=T[node.begin:node.end].upper()
            else:
                G.node[node]['seq']=T[node.begin:node.end]
            if remap:
                mapping[node]=i
                i+=1
    if remap: #get rid of interval objects
        G=nx.relabel_nodes(G,mapping,copy=False)
