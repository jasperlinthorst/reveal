#!/usr/bin/env python

from intervaltree import Interval, IntervalTree
import networkx as nx
import argparse
import logging
import os
import sys
import time

#reveal imports
import schemes
import finish
import plot
import utils
import merge
import subgraph
import comp
import extract
import convert
import align
import realign
import bubbles
import matches

def main():
    desc="""
    Type 'reveal <positional_argument> --help' for help on a specific subcommand.\n
    Reveal constructs population reference graphs by aligning multiple whole 
    genomes using recursive exact matching.
    http://www.biorxiv.org/content/early/2015/07/17/022715.
    """
    
    parser = argparse.ArgumentParser(prog="reveal", usage="reveal -h for usage", description=desc)
    parser.add_argument("-l", "--log-level", type=int, dest="loglevel", default=20, help="Log level: 10=debug 20=info (default) 30=warn 40=error 50=fatal.")
    parser.add_argument("--64", dest="sa64", default=False, action="store_true", help="Use 64bit suffix array in the index.")
    
    subparsers = parser.add_subparsers()
    parser_aln = subparsers.add_parser('align',prog="reveal align", description="Construct a population graph from input genomes or other graphs.")
    parser_plot = subparsers.add_parser('plot', prog="reveal plot", description="Generate mumplot for two fasta files.")
    parser_convert = subparsers.add_parser('convert', prog="reveal convert", description="Convert gfa graph to gml.")
    parser_extract = subparsers.add_parser('extract', prog="reveal extract", description="Extract the input sequence from a graph.")
    parser_comp = subparsers.add_parser('comp', prog="reveal comp", description="Reverse complement the graph.")
    parser_finish = subparsers.add_parser('finish', prog="reveal finish", description="Finish a draft assembly by ordering and orienting contigs with respect to a finished reference assembly.")
    parser_matches = subparsers.add_parser('matches', prog="reveal matches", description="Outputs all (multi) m(u/e)ms.")
    parser_subgraph = subparsers.add_parser('subgraph', prog="reveal subgraph", description="Extract subgraph from gfa by specified node ids.")
    parser_bubbles = subparsers.add_parser('bubbles', prog="reveal bubbles", description="Extract all bubbles from the graph.")
    parser_realign = subparsers.add_parser('realign', prog="reveal realign", description="Realign between two nodes in the graph.")
    parser_merge = subparsers.add_parser('merge', prog="reveal merge", description="Combine multiple gfa graphs into a single gfa graph.")
    
    parser_aln.add_argument('inputfiles', nargs='*', help='Fasta or gfa files specifying either assembly/alignment graphs (.gfa) or sequences (.fasta). When only one gfa file is supplied, variants are called within the graph file.')
    parser_aln.add_argument("-o", "--output", dest="output", help="Prefix of the variant and alignment graph files to produce, default is \"sequence1_sequence2\"")
    #parser_aln.add_argument("-p", dest="pcutoff", type=float, default=1e-3, help="If, the probability of observing a MUM of the observed length by random change becomes larger than this cutoff the alignment is stopped (default 1e-3).")
    parser_aln.add_argument("-t", "--threads", dest="threads", type=int, default=0, help = "The number of threads to use for the alignment.")
    parser_aln.add_argument("-m", dest="minlength", type=int, default=15, help="Min length of an exact match (default 20).")
    parser_aln.add_argument("-c", dest="minscore", type=int, default=0, help="Min score of an exact match (default 0), exact maches are scored by their length and penalized by the indel they create with respect to previously accepted exact matches.")
    parser_aln.add_argument("-n", dest="minn", type=int, default=2, help="Only align graph on exact matches that occur in at least this many samples.")
    parser_aln.add_argument("-e", dest="exp", type=int, default=2, help="Increase \'e\' to prefer shorter matches observed in more genomes, over larger matches in less genomes (default=2).")
    parser_aln.add_argument("--wp", dest="wpen", type=int, default=1, help="Multiply penalty for a MUM by this number in scoring scheme.")
    parser_aln.add_argument("--ws", dest="wscore", type=int, default=3, help="Multiply length of MUM by this number in scoring scheme.")
    parser_aln.add_argument("--mumplot", dest="mumplot", action="store_true", default=False, help="Save a mumplot for the actual aligned chain of anchors (depends on matplotlib).")
    parser_aln.add_argument("-i", dest="interactive", action="store_true", default=False, help="Show an interactive visualisation of the mumplot (depends on matplotlib).")
    
    parser_aln.add_argument("--sa", dest="sa", default="", help="Specify a preconstructed suffix array to decouple suffix array construction.")
    parser_aln.add_argument("--lcp", dest="lcp", default="", help="Specify a preconstructed lcp array to decouple lcp array construction.")
    parser_aln.add_argument("--cache", dest="cache", default=False, action="store_true", help="When specified, it caches the suffix and lcp array to disk after construction.")
    
    parser_aln.add_argument("-g", dest="minsamples", type=int, default=1, help="Only index nodes that occur in this many samples or more (default 1).")
    parser_aln.add_argument("-x", dest="maxsamples", type=int, default=None, help="Only align nodes that have maximally this many samples (default None).")
    parser_aln.add_argument("-r", dest="reference", type=str, default=None, help="Name of the sequence that should be used as a coordinate system or reference.")
    parser_aln.add_argument("-s", dest="targetsample", type=str, default=None, help="Only align nodes in which this sample occurs.")
    parser_aln.add_argument("--gml", dest="gml", action="store_true", default=False, help="Produce a gml graph instead gfa.")
    parser_aln.add_argument("--gml-max", dest="hwm", default=4000, help="Max number of nodes per graph in gml output.")
    parser_aln.add_argument("--nometa", dest="nometa", action="store_true", default=False, help="Produce a gfa graph without node annotations, to ensure it's parseable by other programs.")
    parser_aln.add_argument("--paths", dest="paths", action="store_true", default=False, help="Output paths in GFA.")
    parser_aln.set_defaults(func=align.align_cmd)
    
    parser_extract.add_argument('graph', nargs=1, help='gfa file specifying the graph from which the genome should be extracted.')
    parser_extract.add_argument('samples', nargs='*', help='Name of the sample to be extracted from the graph.')
    parser_extract.add_argument("--width", dest="width", type=int, default=100 , help='Line width for fasta output.')
    parser_extract.set_defaults(func=extract.extract_cmd)
    
    parser_plot.add_argument('fastas', nargs='*', help='Two fasta files for which a mumplot should be generated.')
    parser_plot.add_argument("-m", dest="minlength", type=int, default=100, help="Minimum length of exact matches to vizualize (default=100).")
    parser_plot.add_argument("-i", dest="interactive", action="store_true", default=False, help="Wheter to produce interactive plots which allow zooming on the dotplot (default=False).")
    parser_plot.add_argument("-u", dest="uniq", action="store_true", default=False, help="Plot only maximal unique matches.")
    parser_plot.add_argument("-r", dest="region", default=None, help="Highlight interval (as \"<start>:<end>\") with respect to x-axis (first sequence).")
    parser_plot.set_defaults(func=plot.plot)
    
    parser_comp.add_argument('graph', nargs=1, help='The graph to be reverse complemented.')
    parser_comp.set_defaults(func=comp.comp_cmd)
    
    parser_finish.add_argument('reference', help='Graph or sequence to which query/contigs should be assigned.')
    parser_finish.add_argument('contigs', help='Graph or fasta that is to be reverse complemented with respect to the reference.')
    parser_finish.add_argument("-m", dest="minlength", type=int, default=100, help="Min length of maximal exact matches for considering (default 20).")
    parser_finish.add_argument("-i", dest="interactive", action="store_true", default=False, help="Output interactive plot.")
    parser_finish.add_argument("--plot", dest="plot", action="store_true", default=False, help="Output mumplots for the \'finished\' chromosomes (depends on matplotlib).")
    parser_finish.add_argument("--plotall", dest="plotall", action="store_true", default=False, help="Plot all matches, instead of only the chained matches.")
    parser_finish.add_argument("--split", dest="split", action="store_true", default=False, help="Split the \'finished\' genome by chromosome.")
    parser_finish.add_argument("--nofix", dest="fix", action="store_false", default=True, help="Don't attempt to fix misassembled contigs.")
    parser_finish.add_argument("--fixedgapsize", dest="fixedsize", action="store_true", default=False, help="Use fixed gapsizes of 100 N.")
    parser_finish.add_argument("--cache", dest="cache", default=False, action="store_true", help="When specified, it caches the text, suffix and lcp array to disk after construction.")
    parser_finish.add_argument("--sa1", dest="sa1", default="", help="Specify a preconstructed suffix array for extracting matches between the two genomes in their current orientation.")
    parser_finish.add_argument("--lcp1", dest="lcp1", default="", help="Specify a preconstructed lcp array for extracting matches between the two genomes in their current orientation.")
    parser_finish.add_argument("--sa2", dest="sa2", default="", help="Specify a preconstructed suffix array for extracting matches between the two genomes in which the sequence of the contigs was reverse complemented.")
    parser_finish.add_argument("--lcp2", dest="lcp2", default="", help="Specify a preconstructed lcp array for extracting matches between the two genomes in which the sequence of the contigs was reverse complemented.")
    parser_finish.set_defaults(func=finish.finish)

    parser_matches.add_argument('reference', help='Graph or sequence to which query/contigs should be assigned.')
    parser_matches.add_argument('contigs', help='Graph or fasta that is to be reverse complemented with respect to the reference.')    
    parser_matches.add_argument("-m", dest="minlength", type=int, default=100, help="Min length of maximal exact matches for considering (default 20).")
    parser_matches.add_argument("--cache", dest="cache", default=False, action="store_true", help="When specified, it caches the text, suffix and lcp array to disk after construction.")
    parser_matches.add_argument("--sa1", dest="sa1", default="", help="Specify a preconstructed suffix array for extracting matches between the two genomes in their current orientation.")
    parser_matches.add_argument("--lcp1", dest="lcp1", default="", help="Specify a preconstructed lcp array for extracting matches between the two genomes in their current orientation.")
    parser_matches.add_argument("--sa2", dest="sa2", default="", help="Specify a preconstructed suffix array for extracting matches between the two genomes in which the sequence of the contigs was reverse complemented.")
    parser_matches.add_argument("--lcp2", dest="lcp2", default="", help="Specify a preconstructed lcp array for extracting matches between the two genomes in which the sequence of the contigs was reverse complemented.")    
    parser_matches.set_defaults(func=matches.matches)
    
    parser_convert.add_argument('graphs', nargs='*', help='The gfa graph to convert to gml.')
    parser_convert.add_argument("-n", dest="minsamples", type=int, default=1, help="Only align nodes that occcur in this many samples (default 1).")
    parser_convert.add_argument("-x", dest="maxsamples", type=int, default=None, help="Only align nodes that have maximally this many samples (default None).")
    parser_convert.add_argument("-s", dest="targetsample", type=str, default=None, help="Only align nodes in which this sample occurs.")
    parser_convert.add_argument("--gml-max", dest="hwm", default=4000, help="Max number of nodes per graph in gml output.")
    parser_convert.add_argument("--gfa",  action="store_true", dest="gfa", default=False, help="Rewrite gfa file.")
    parser_convert.add_argument("--partition",  action="store_true", dest="partition", default=False, help="Output graph as multiple subgraphs if possible.")
    parser_convert.set_defaults(func=convert.convert)
    
    parser_subgraph.add_argument('inputfiles', nargs='*', help='The gfa graph followed by node ids that make up the subgraph.')
    parser_subgraph.add_argument("-o", dest="outfile", type=str, default="~tmp", help="Prefix of the file to which subgraph will be written.")
    parser_subgraph.add_argument("--gml", dest="gml", action="store_true", default=False, help="Produce a gml graph instead of gfa.")
    parser_subgraph.set_defaults(func=subgraph.subgraph)
    
    parser_bubbles.add_argument("graph", nargs=1, help='Graph in gfa format from which bubbles are to be extracted.')
    parser_bubbles.add_argument("-r", dest="reference", type=str, default=None, help="Name of the sequence that, if possible, should be used as a coordinate system or reference.")
    parser_bubbles.set_defaults(func=bubbles.bubbles_cmd)
    
    parser_realign.add_argument("graph", nargs=1, help='Graph in gfa format for which a bubble should be realigned.') 
    parser_realign.add_argument("source", nargs='?', type=int, help='Source node.')
    parser_realign.add_argument("sink", nargs='?', type=int, help='Sink node.')
    parser_realign.add_argument("-m", dest="minlength", type=int, default=20, help="Min length of an exact match (default 20).")
    parser_realign.add_argument("-c", dest="minscore", type=int, default=0, help="Min score of an exact match (default 0), exact maches are scored by their length and penalized by the indel they create with respect to previously accepted exact matches.")
    parser_realign.add_argument("-n", dest="minn", type=int, default=2, help="Only align graph on exact matches that occur in at least this many samples.")
    parser_realign.add_argument("-o", dest="outfile", type=str, default=None, help="File to which realigned graph is to be written.")
    parser_realign.add_argument("--all", action="store_true", dest="all", default=False, help="Trigger realignment for all complex bubbles.")
    parser_realign.add_argument("--maxlen", dest="maxlen", type=int, default=10000000, help="Maximum length of the cumulative sum of all paths that run through the complex bubble.")
    parser_realign.add_argument("--maxsize", dest="maxsize", type=int, default=500, help="Maximum allowed number of nodes that are contained in a complex bubble.")
    parser_realign.set_defaults(func=realign.realign_bubble_cmd)
    
    parser_merge.add_argument("graphs", nargs='*', help='Graphs in gfa format that should be merged.')
    parser_merge.add_argument("-o", dest="outprefix", type=str, default=None, help="Prefix of the file to which merged graph is written.")
    parser_merge.set_defaults(func=merge.merge_cmd)
    
    args = parser.parse_args()
    logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p', level=args.loglevel)
    logging.logThreads = 0
    
    args.func(args)
