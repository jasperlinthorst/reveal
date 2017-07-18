import networkx as nx
import logging
from intervaltree import Interval, IntervalTree
import sys
import os

def fasta_reader(fn,truncN=False,toupper=True):
    seq=""
    with open(fn,'r') as ff:
        for line in ff:
            if line.startswith(">"):
                if seq!="":
                    yield name,seq
                name=line.rstrip().replace(">","")
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
                else:
                    if toupper:
                        seq+=line.upper().rstrip()
                    else:
                        seq+=line.rstrip()
        if seq!="":
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

def rc(seq):
    d = {'A':'T','C':'G','G':'C','T':'A','N':'N','a':'t','c':'g',\
        'g':'c','t':'a','n':'n','Y':'R','R':'Y','K':'M','M':'K',\
        'S':'S','W':'W','B':'V','V':'B','D':'H','H':'D','N':'N',\
        'X':'X','-':'-'}
    return "".join([d[b] for b in reversed(seq)])

def read_gfa(gfafile, index, tree, graph, minsamples=1, maxsamples=None, targetsample=None, revcomp=False):
    f=open(gfafile,'r')
    sep=";"
    mapping={} #temp mapping object
    edges=[] #tmp list for edges
    id2sample={}
    sample2id={}
    i=0
    for line in f:
        if line.startswith('H'):
            h=line.split()
            tag=h[1].split(':')
            if tag[0]=="ORI":
                for sample in tag[2].rstrip(sep).split(sep):
                    id2sample[i]=sample
                    sample2id[sample]=i
                    i+=1
                    if 'samples' in graph.graph:
                        graph.graph['samples'].append(sample)
                    else:
                        graph.graph['samples']=[sample]
        
        if line.startswith('S'):
            s=line.strip().split('\t')
            nodeid=int(s[1])
            gnodeid=graph.number_of_nodes()+1
            
            if graph.has_node(gnodeid):
                logging.fatal("Id space for nodes is larger than total number of nodes in the graph.")
                return
            
            ann={}
            
            if len(s)>3:
                for v in s[4:]:
                    v=v.split(':')
                    if v[2].find(sep)!=-1: #multi-valued
                        ann[v[0]]=v[2].rstrip(sep).split(sep)
                    else:
                        ann[v[0]]=v[2]
            
            if "ORI" not in ann: #not a reveal graph, so no metadata on nodes, just index all
                if index!=None:
                    if len(s)<3:
                        intv=None
                    else:
                        if revcomp:
                            intv=index.addsequence(rc(s[2]).upper())
                        else:
                            intv=index.addsequence(s[2].upper())
                        intv=Interval(intv[0],intv[1])
                        tree.add(intv)
                    graph.add_node(intv,offsets={gfafile:0},aligned=0)
                    mapping[nodeid]=intv
                else:
                    if len(s)<3:
                        s.append("")
                    graph.add_node(gnodeid,offsets={gfafile:0},seq=s[2].upper(),aligned=0)
                    mapping[nodeid]=gnodeid

            else: #there are annotations on the nodes, so use them
                if not(isinstance(ann['ORI'], list)):
                    ann['ORI']=[ann['ORI']]
                if not(isinstance(ann['OFFSETS'], list)):
                    ann['OFFSETS']=[ann['OFFSETS']]
                offsets=dict()
                if 'OFFSETS' in ann: #create dictionary by sample name
                    for sampleid,offset in zip(ann['ORI'],ann['OFFSETS']):
                        offsets[id2sample[int(sampleid)]]=int(offset)
                tmp=set()
                for ori in ann['ORI']: #map sample ids back to sample names
                    tmp.add(id2sample[int(ori)])
                ann['ORI']=tmp
                
                assert(isinstance(ann['ORI'],set))
                
                if len(ann['ORI'])<minsamples: #dont index these nodes, just add them to the graph
                    graph.add_node(gnodeid,seq=s[2].upper(),aligned=0,offsets=offsets)
                    mapping[nodeid]=gnodeid
                elif maxsamples!=None and len(ann['ORI'])>maxsamples: #dont index these nodes, just add them to the graph
                    graph.add_node(gnodeid,seq=s[2].upper(),aligned=0,offsets=offsets)
                    mapping[nodeid]=gnodeid
                elif (targetsample!=None and targetsample not in ann['ORI']): #dont index these nodes, just add them to the graph
                    graph.add_node(gnodeid,seq=s[2].upper(),aligned=0,offsets=offsets)
                    mapping[nodeid]=gnodeid
                else:
                    if index!=None:
                        if revcomp:
                            intv=index.addsequence(rc(s[2]).upper())
                        else:
                            intv=index.addsequence(s[2].upper())
                        intv=Interval(intv[0],intv[1])
                        tree.add(intv)
                        graph.add_node(intv,aligned=0,offsets=offsets)
                        mapping[nodeid]=intv
                    else:
                        graph.add_node(gnodeid,seq=s[2].upper(),aligned=0,offsets=offsets)
                        mapping[nodeid]=gnodeid
        
        #L      206     +       155     +       0M
        if line.startswith('L'):
            edges.append(line)

    for line in edges:
        e=line.strip().split()
        assert(not graph.has_edge(mapping[int(e[1])],mapping[int(e[3])]))
        assert(not graph.has_edge(mapping[int(e[3])],mapping[int(e[1])]))
        graph.add_edge(mapping[int(e[1])],mapping[int(e[3])],ofrom=e[2],oto=e[4],cigar=e[5])
    
    if revcomp:
        genome2length=dict()
        #relabel the offsets, determine the length of all genomes in the graph, then l-pos
        for sample in graph.graph['samples']:
            maxp=0
            for node,data in graph.nodes(data=True):
                if sample2id[sample] in data['offsets']:
                    if data['offsets'][sample2id[sample]]+ (node[1]-node[0]) >maxp:
                        maxp=data['offsets'][sample2id[sample]]+(node[1]-node[0])
            genome2length[sample]=maxp
        
        for sample in graph.graph['samples']:
            for node,data in graph.nodes(data=True):
                if sample2id[sample] in data['offsets']:
                    graph.node[node]['offsets'][sample2id[sample]]=genome2length[sample]-(graph.node[node]['offsets'][sample2id[sample]]+(node[1]-node[0]))
        
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

def write_gfa(G,T,outputfile="reference.gfa", nometa=False, paths=False, remap=True):
    
    if not outputfile.endswith(".gfa"):
        outputfile+=".gfa"
    
    f=open(outputfile,'wb')
    sep=';'
    f.write('H\tVN:Z:1.0\n')
    f.write('H\tORI:Z:')
    
    sample2id=dict()
    i=0
    
    for sample in G.graph['samples']:
        for subsample in sample.rstrip().rstrip(sep).split(sep):
            f.write(subsample+sep)
            sample2id[subsample]=i
            i+=1
    
    f.write('\n')
    
    mapping={}
    
    if remap:
        for i,node in enumerate(nx.topological_sort(G)): #iterate once to get a mapping of ids to intervals
            mapping[node]=i+1 #one-based for vg
    else:
        for i,node in enumerate(nx.topological_sort(G)): #iterate once to get a mapping of ids to intervals
            mapping[node]=node

    for i,node in enumerate(nx.topological_sort(G)):
        if isinstance(node,str):
            if node=='start' or node=='end':
                continue

        i+=1
        data=G.node[node]
        seq=""
        if 'seq' in data:
            #f.write('S\t'+str(i)+'\t'+data['seq'].upper())
            f.write('S\t'+str(mapping[node])+'\t'+data['seq'].upper())
            seq=data['seq']
        else:
            if isinstance(node,Interval):
                seq=T[node.begin:node.end].upper()
                #f.write('S\t'+str(i)+'\t'+seq)
                f.write('S\t'+str(mapping[node])+'\t'+seq)
            elif isinstance(node,tuple):
                seq=T[node[0]:node[0]+G.node[node]['l']].upper()
                #f.write('S\t'+str(i)+'\t'+seq)
                f.write('S\t'+str(mapping[node])+'\t'+seq)
            else:
                #f.write('S\t'+str(i)+'\t')
                f.write('S\t'+str(mapping[node])+'\t')
        
        if nometa or 'offsets' not in data:
            f.write("\n")
        else:
            tmp=data['offsets'].keys()
            tmp.sort(key=lambda s: sample2id[s])
            data['sample']=tmp
            
            if not isinstance(data['offsets'],dict):
                print "ERROR offset data:", data['offsets']
            
            f.write("\t*\tORI:Z:%s\tOFFSETS:Z:%s\tRC:i:%s\n" % 
                    (
                    sep.join([str(sample2id[s]) for s in data['offsets']]),
                    sep.join([str(data['offsets'][s]) for s in data['offsets']]),
                    len(data['offsets'])*len(seq)
                    )
                )
         
        for to in G[node]:
            f.write('L\t'+str(mapping[node])+'\t+\t'+str(mapping[to])+"\t+\t0M\n")
    
    if paths:
        for sample in G.graph['samples']:
            path=[]
            for node in nx.topological_sort(G):
                if sample in G.node[node]['offsets']:
                    i=mapping[node]
                    path.append(str(i)+"+")
            f.write("P\t"+sample+"\t"+",".join(path)+"\t"+",".join(["0M"]*len(path))+"\n")
    
    f.close()

def write_gml(G,T,outputfile="reference",partition=True,hwm=4000):
    G=G.copy()
    mapping={}
    totn=len(G.graph['samples'])

    for key in G.graph:
        G.graph[key]=str(G.graph[key])
    
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
        i=0
        for subset in nx.connected_components(G.to_undirected()):
            sgn=[]
            g=G.subgraph(subset)
            gn=G.graph['samples']
            for n in nx.topological_sort(g):
                if sgn==[]:
                    fr=n
                sgn.append(n)
                if 'offsets' in G.node[n]:
                    if len(G.node[n]['offsets'])==totn: #join/split node
                        if len(sgn)>=hwm:
                            sg=G.subgraph(sgn)
                            fn=outputfile+'.'+str(i)+'.gml'
                            nx.write_gml(sg,fn)
                            outputfiles.append(fn)
                            sgn=[]
                            i+=1
            if len(sgn)>0:
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
