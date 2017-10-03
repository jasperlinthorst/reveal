import networkx as nx
import utils
import sys
import logging
import os

def convert(args):
    for graph in args.graphs:
        #g=nx.DiGraph()
        g=nx.MultiDiGraph()
        g.graph['samples']=[]
        g.graph['sample2id']=dict()
        if graph.endswith(".gfa"): #gfa to gml/gfa
            utils.read_gfa(graph,None,None,g,minsamples=args.minsamples,
                                 maxsamples=args.maxsamples,
                                 targetsample=args.targetsample)
            if args.gfa:
                fn=graph.replace(".gfa",".rewrite.gfa")
                graph=utils.write_gfa(g,"", outputfile=fn)
                logging.info("gfa graph written to: %s"%fn)
            else:
                fn=utils.write_gml(g,"", hwm=args.hwm, outputfile=graph.replace(".gfa",""), partition=args.partition)
                logging.info("gml graph written to: %s"%fn)
        elif graph.endswith(".fa") or graph.endswith(".fasta") or graph.endswith(".fna"): #assume fasta to gfa
            i=0
            for name,seq in utils.fasta_reader(graph):
                g.graph['samples'].append(os.path.basename(graph))
                g.graph['sample2id'][os.path.basename(graph)]=0
                g.add_node(i,offsets={0:0},seq=seq)
                i+=1
            filename=graph[:graph.rfind(".")]+".gfa"
            utils.write_gfa(g,"", outputfile=filename)
            logging.info("gfa graph written to: %s"%filename)
        else:
            logging.fatal("Unknown filetype, need gfa or fasta extension.")
            return
