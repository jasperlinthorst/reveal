import networkx as nx
import utils
import sys
import logging
import os

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
                                 targetsample=args.targetsample)
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

        elif graph.endswith(".fa") or graph.endswith(".fasta") or graph.endswith(".fna"): #assume fasta to gfa
            i=0
            for name,seq in utils.fasta_reader(graph):
                g.graph['paths'].append(os.path.basename(graph))
                g.graph['path2id'][os.path.basename(graph)]=0
                g.graph['id2path'][0]=os.path.basename(graph)
                g.add_node(i,offsets={0:0},seq=seq)
                i+=1
            filename=graph[:graph.rfind(".")]+".gfa"
            utils.write_gfa(g,"", outputfile=filename)
            logging.info("gfa graph written to: %s"%filename)
        else:
            logging.fatal("Unknown filetype, need gfa or fasta extension.")
            return

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
        for node in nx.topological_sort(G):
            if type(node)!=str:
                go=max([0]+[G.node[pred]['graphoffset']+len(G.node[pred]['seq']) for pred in G.predecessors(node) if type(pred)!=str])
                G.node[node]['graphoffset']=go
                for k in G.node[node]['offsets']:
                    if G.node[node]['offsets'][k]+len(G.node[node]['seq'])>sizes[k]:
                        sizes[k]=G.node[node]['offsets'][k]+len(G.node[node]['seq'])

        maf.write("a score=?\n")
        for path in G.graph['paths']:
            o=0
            maf.write("s %s %d %d <strand> %-10d "%(path.ljust(max([len(p) for p in G.graph['paths']])),0,go,sizes[G.graph['path2id'][path]]))
            sid=G.graph['path2id'][path]
            for node in nx.topological_sort(G):
                if type(node)!=str and sid in G.node[node]['offsets']:
                    while o<G.node[node]['graphoffset']:
                        maf.write("-")
                        o+=1
                    maf.write("%s"%G.node[node]['seq'])
                    o+=len(G.node[node]['seq'])
            
            maf.write("\n")




