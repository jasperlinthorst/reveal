import networkx as nx
import utils
import sys
import logging
import os
import uuid

def convert(args):
    for graph in args.graphs:
        
        if args.nocycles:
            g=nx.DiGraph()
        else:
            g=nx.MultiDiGraph()

        g.graph['paths']=[]
        g.graph['path2id']=dict()
        g.graph['id2path']=dict()

        if graph.endswith(".gfa"): #gfa to gml/gfa
            utils.read_gfa(graph,None,None,g,minsamples=args.minsamples,
                                 maxsamples=args.maxsamples,
                                 targetsample=args.targetsample,
                                 remap=False)
            if args.type=="gfa":
                fn=graph.replace(".gfa",".rewrite.gfa")
                graph=utils.write_gfa(g,"", outputfile=fn)
                logging.info("gfa graph written to: %s"%fn)
            elif args.type=="gml":
                fn=utils.write_gml(g,"", hwm=args.hwm, outputfile=graph.replace(".gfa",""), partition=args.partition)
                logging.info("gml graph written to: %s"%fn)
            elif args.type=="maf":
                logging.info("Converting graph to maf..")
                graph2maf(g,graph.replace(".gfa",".maf"))
        
        elif graph.endswith(".maf"): #multiple alignment format, convert to graph
            g=maf2graph(graph)
            filename=graph[:graph.rfind(".")]+".gml"
            utils.write_gml(g,"", outputfile=filename)

            filename=graph[:graph.rfind(".")]+".gfa"
            utils.write_gfa(g,"", outputfile=filename)
            logging.debug("gfa graph written to: %s"%filename)

        elif graph.endswith(".fa") or graph.endswith(".fasta") or graph.endswith(".fna"): #assume fasta to gfa
            if args.aligned:
                seqs=[]
                names=[]
                for name,seq in utils.fasta_reader(graph,keepdash=True):
                    seqs.append(seq)
                    names.append(name)
                g,nid=utils.aln2graph(seqs,names)
            else:
                i=0
                start=uuid.uuid4().hex
                end=uuid.uuid4().hex
                g.graph['startnodes']=[start]
                g.graph['endnodes']=[end]
                g.add_node(start,offsets=dict())
                g.add_node(end,offsets=dict())
                for i,v in enumerate(utils.fasta_reader(graph)):
                    name,seq=v
                    g.graph['paths'].append(name)
                    g.graph['path2id'][name]=i
                    g.graph['id2path'][i]=name
                    g.node[start]['offsets'][i]=0
                    g.node[end]['offsets'][i]=len(seq)
                    g.add_node(i,offsets={i:0},seq=seq)
                    g.add_edge(start,i,paths=set([i]))
                    g.add_edge(i,end,paths=set([i]))

            filename=graph[:graph.rfind(".")]+".gfa"
            utils.write_gfa(g,"", outputfile=filename)
            logging.debug("gfa graph written to: %s"%filename)
        else:
            logging.fatal("Unknown filetype, need gfa or fasta extension.")
            return

#converts a multiple alignment format file to a graph
def maf2graph(maffile):
    files=set()
    G=nx.MultiDiGraph()

    startnode=uuid.uuid4().hex
    endnode=uuid.uuid4().hex

    G.graph['startnodes']=set([startnode])
    G.graph['endnodes']=set([endnode])
    G.graph['path2id']=dict()

    G.add_node(startnode,offsets=dict())
    G.add_node(endnode,offsets=dict())

    nid=0
    with open(maffile,"r") as maf:
        for line in maf:
            if line.startswith("#"):
                continue
            elif line.startswith("a"): #start of an aligned segment
                nid+=1
                G.add_node(nid,data=dict())
            elif line.startswith("s"):
                cols=line.rstrip().split()
                if '.' in cols[1]: #TODO: use db parameter to specificy a single mfa file with all sequence
                    file,name=cols[1][:cols[1].find('.')],cols[1][cols[1].find('.')+1:]
                    files.add(file)
                else:
                    file=None #args.db?
                    name=cols[1]

                if name not in G.graph['path2id']:
                    G.graph['path2id'][name]=len(G.graph['path2id'])
                    G.node[startnode]['offsets'][G.graph['path2id'][name]]=0

                G.node[nid]['data'][(file,name)]={'start':int(cols[2]),
                                                  'end':int(cols[2])+int(cols[3]),
                                                  'orientation':cols[4],
                                                  'aln':cols[6]
                                                  }
    nid+=1

    remove=[]
    for node,d in G.nodes(data=True):
        if 'data' in d and len(d['data'])==1: #multiplicity of 1, strictly not an alignment
            remove.append(node)

    G.remove_nodes_from(remove)

    db=dict() #map name to sequence
    for file in files:
        for name,seq in utils.fasta_reader(file+".fasta"): #guess that the original file has a ".fasta" extension
            name=name.split()[0]
            key=(file,name)
            if key in db:
                logging.fatal("Non unique contig-name: %s. quit."%name)
                sys.exit(1)
            else:
                db[key]=seq

    remove=[]

    #for every sequence, check that none of the alignments overlap, otherwise assignment is not 1-1
    for file,name in db:
        seq=db[(file,name)]

        intvs=[]
        for node in G:
            if 'data' in G.node[node]: #does the node represent an aligned segment?
                if (file,name) in G.node[node]['data']:
                    intvs.append((G.node[node]['data'][(file,name)]['start'] , G.node[node]['data'][(file,name)]['end'], node))
        
        intvs.sort() #sort by start position
        pstart=0
        pend=0
        pnode=startnode
        unaligned=[]

        for start,end,node in intvs:
            if start>pend:
                unaligned.append((pend,start))
                G.add_node(nid,intv=(pend,start),seq=seq[pend:start])
                G.add_edge(pnode,nid,paths=set([G.graph['path2id'][name]]),ofrom="+",oto="+")
                G.add_edge(nid,node,paths=set([G.graph['path2id'][name]]),ofrom="+",oto="+")
                nid+=1
            elif start<pend:
                logging.fatal("Overlapping alignments for sequence: %s.%s --> (%d,%d) and (%d,%d)."%(file,name,pstart,pend,start,end))
                remove.append(node)
                # sys.exit(1)
            else: #no gap, just connect subsequent intervals
                G.add_edge(pnode,node,paths=set([G.graph['path2id'][name]]),ofrom="+",oto="+")

            pstart,pend,pnode=start,end,node
        
        if len(seq)!=pend:
            unaligned.append((pend,len(seq)))
            G.add_node(nid,intv=((pend,len(seq))),seq=seq[pend:len(seq)])
            G.add_edge(pnode,nid,paths=set([G.graph['path2id'][name]]),ofrom="+",oto="+")
            G.add_edge(nid,endnode,paths=set([G.graph['path2id'][name]]),ofrom="+",oto="+")
            nid+=1
        else:
            G.add_edge(pnode,endnode,paths=set([G.graph['path2id'][name]]),ofrom="+",oto="+")

    G.remove_nodes_from(remove)

    # print "Unaligned segments",unaligned

    alignments=[node for node in G if 'data' in G.node[node]]
    
    for node in alignments: #expand all alignments in the graph

        if 'data' in G.node[node]:
            seqs=[]
            names=[]
            offsets={}

            for file,name in G.node[node]['data']:
                seqs.append(G.node[node]['data'][(file,name)]['aln'])
                offsets[G.graph['path2id'][name]]=G.node[node]['data'][(file,name)]['start']
                names.append(name)

            sg,nid=utils.aln2graph(seqs,names,idoffset=nid,path2id=G.graph['path2id'],offsets=offsets)

            nid+=1

            G.add_nodes_from(sg.nodes(data=True))
            G.add_edges_from(sg.edges(data=True))

            assert(len(sg.graph['startnodes'])==1)
            assert(len(sg.graph['endnodes'])==1)

            sgstart=sg.graph['startnodes'][0]
            sgend=sg.graph['endnodes'][0]

            for v,t,d in G.in_edges(node,data=True):
                G.add_edge(v,sgstart,paths=d['paths'],ofrom="+",oto="+")

            for v,t,d in G.out_edges(node,data=True):
                G.add_edge(sgend,t,paths=d['paths'],ofrom="+",oto="+")

            #hack this in here so we can continue
            G.node[sgstart]['seq']=""
            G.node[sgend]['seq']=""
            nx.relabel_nodes(G,{sgstart: nid, sgend: nid+1},copy=False)

            nid+=2

            G.remove_node(node)

    return G

def graph2maf(G,filename):

    if isinstance(G,nx.MultiDiGraph):
        #TODO: decompose global alignment into local alignments by deconnecting structure edges
        #determine set of structure edges
        orgpaths=set([G.graph['path2id'][p] for p in G.graph['paths'] if p.startswith('*')])
        refpaths=set([G.graph['path2id'][p] for p in G.graph['paths'] if not p.startswith('*')])
        
        es=[]
        for e0,e1,d in G.edges(data=True):
            if len(d['paths'] & refpaths)==0: #edge that exclusively represents structural event
                es.append((e0,e1))

        toremove=es
        G.remove_edges_from(toremove)

    sizes={sid:0 for sid in G.graph['id2path']}
    
    with open(filename,'w') as maf:

        for g in nx.weakly_connected_component_subgraphs(G):
            
            longest=0
            sids=set()
            for node in nx.topological_sort(g):
                if type(node)!=str:
                    go=max([0]+[G.node[pred]['graphoffset']+len(G.node[pred]['seq']) for pred in G.predecessors(node) if type(pred)!=str])
                    G.node[node]['graphoffset']=go
                    
                    if go+len(G.node[node]['seq'])>longest:
                        longest=go+len(G.node[node]['seq'])

                    for k in G.node[node]['offsets']:
                        sids.add(k)
                        if G.node[node]['offsets'][k]+len(G.node[node]['seq'])>sizes[k]:
                            sizes[k]=G.node[node]['offsets'][k]+len(G.node[node]['seq'])
            
            ml=max([len(p) for p in G.graph['paths']])

            maf.write("##maf version=1\n")
            maf.write("a\n")
            for sid in sids:
                path=G.graph['id2path'][sid]
                o=0
                sl=0
                maf.write("s %s %d %d + %-10d "%(path.ljust(ml), 0, sizes[G.graph['path2id'][path]], sizes[G.graph['path2id'][path]]) )
                for node in nx.topological_sort(g):
                    if type(node)!=str and sid in G.node[node]['offsets']:
                        while o<G.node[node]['graphoffset']:
                            maf.write("-")
                            o+=1
                        sl+=len(G.node[node]['seq'].replace("-",""))
                        maf.write("%s"%G.node[node]['seq'])
                        o+=len(G.node[node]['seq'])
                
                maf.write("-"*(longest-o)) #pad with dash so all lines are equally long
                maf.write("\n")
            maf.write("\n")



