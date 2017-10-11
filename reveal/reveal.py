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
import chain
import stats
import split

def main():
    desc="""
    Type 'reveal <positional_argument> --help' for help on a specific subcommand.\n
    Reveal constructs population reference graphs by aligning multiple whole 
    genomes using recursive exact matching.
    http://www.biorxiv.org/content/early/2015/07/17/022715.
    """
    
    parser = argparse.ArgumentParser(prog="reveal", usage="reveal -h for usage", description=desc, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-l", "--log-level", type=int, dest="loglevel", default=20, help="Log level: 10=debug 20=info 30=warn 40=error 50=fatal.")
    parser.add_argument("--64", dest="sa64", default=False, action="store_true", help="Use 64bit suffix array in the index.")
    
    subparsers = parser.add_subparsers()
    parser_aln = subparsers.add_parser('align',prog="reveal align", description="Construct a population graph from input genomes or other graphs.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_plot = subparsers.add_parser('plot', prog="reveal plot", description="Generate a mumplot that shows all mums between two fasta files.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_convert = subparsers.add_parser('convert', prog="reveal convert", description="Convert gfa graph to gml.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_extract = subparsers.add_parser('extract', prog="reveal extract", description="Extract the input sequence from a graph.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_comp = subparsers.add_parser('comp', prog="reveal comp", description="Reverse complement the graph.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_finish = subparsers.add_parser('finish', prog="reveal finish", description="Finish a draft assembly by ordering and orienting contigs with respect to a finished reference assembly.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    #parser_matches = subparsers.add_parser('matches', prog="reveal matches", description="Outputs all (multi) m(u/e)ms.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_subgraph = subparsers.add_parser('subgraph', prog="reveal subgraph", description="Extract subgraph from gfa by specified node ids.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_bubbles = subparsers.add_parser('bubbles', prog="reveal bubbles", description="Extract all bubbles from the graph.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_variants = subparsers.add_parser('variants', prog="reveal variants", description="Extract variant calls from the graph.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_realign = subparsers.add_parser('realign', prog="reveal realign", description="Realign between two nodes in the graph.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_merge = subparsers.add_parser('merge', prog="reveal merge", description="Combine multiple gfa graphs into a single gfa graph.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_chain = subparsers.add_parser('chain', prog="reveal chain", description="Use default chaining scheme to construct GFA graph based on a global multi-alignment of all input genomes.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_stats = subparsers.add_parser('stats', prog="reveal stats", description="Output statistics (number of node, edges, genomes etc.) for a graph.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_gplot = subparsers.add_parser('gplot', prog="reveal gplot", description="Generate a plot that represents the alignment of two samples in a graph.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_split = subparsers.add_parser('split', prog="reveal split", description="Split a graph file into a connected component per file.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser_aln.add_argument('inputfiles', nargs='*', help='Fasta or gfa files specifying either assembly/alignment graphs (.gfa) or sequences (.fasta).')
    parser_aln.add_argument("-o", "--output", dest="output", help="Prefix of the variant and alignment graph files to produce, default is \"sequence1_sequence2\"")
    parser_aln.add_argument("-p", dest="pcutoff", type=float, default=None, help="If, the probability of observing a MUM of the observed length by random change becomes larger than this cutoff the alignment is stopped.")
    parser_aln.add_argument("-t", "--threads", dest="threads", type=int, default=0, help = "The number of threads to use for the alignment.")
    parser_aln.add_argument("-m", dest="minlength", type=int, default=20, help="Min length of an exact match.")
    parser_aln.add_argument("-c", dest="minscore", type=int, default=None, help="Min score of an exact match, exact maches are scored by their length and penalized by the indel they create with respect to previously accepted exact matches.")
    parser_aln.add_argument("-n", dest="minn", type=int, default=None, help="Only align graph on exact matches that occur in at least this many samples (if not set, equal to total number of genomes in the (sub)index).")
    parser_aln.add_argument("-e", dest="exp", type=int, default=1, help="Increase \'e\' to prefer shorter matches observed in more genomes, over larger matches in less genomes.")
    parser_aln.add_argument("--gcmodel", dest="gcmodel", choices=["sumofpairs","star-avg","star-med"], default="sumofpairs", help="Which gap-cost model to use for multi-alignment.")
    parser_aln.add_argument("--wp", dest="wpen", type=int, default=1, help="Multiply penalty for a MUM by this number in scoring scheme.")
    parser_aln.add_argument("--ws", dest="wscore", type=int, default=3, help="Multiply length of MUM by this number in scoring scheme.")
    parser_aln.add_argument("--plot", dest="mumplot", action="store_true", default=False, help="Save a mumplot for the actual aligned chain of anchors (depends on matplotlib).")
    parser_aln.add_argument("-i", dest="interactive", action="store_true", default=False, help="Show an interactive visualisation of the mumplot (depends on matplotlib).")
    
    parser_aln.add_argument("--sa", dest="sa", default="", help="Specify a preconstructed suffix array to decouple suffix array construction.")
    parser_aln.add_argument("--lcp", dest="lcp", default="", help="Specify a preconstructed lcp array to decouple lcp array construction.")
    parser_aln.add_argument("--cache", dest="cache", default=False, action="store_true", help="When specified, it caches the suffix and lcp array to disk after construction.")
    
    parser_aln.add_argument("-g", dest="minsamples", type=int, default=1, help="Only index nodes that occur in this many samples or more.")
    parser_aln.add_argument("-x", dest="maxsamples", type=int, default=None, help="Only align nodes that have maximally this many samples.")
    parser_aln.add_argument("-r", dest="reference", type=str, default=None, help="Name of the sequence that should be used as a coordinate system or reference.")
    parser_aln.add_argument("-s", dest="targetsample", type=str, default=None, help="Only align nodes in which this sample occurs.")
    parser_aln.add_argument("--gml", dest="gml", action="store_true", default=False, help="Produce a gml graph instead gfa.")
    parser_aln.add_argument("--gml-max", dest="hwm", default=4000, help="Max number of nodes per graph in gml output.")
    parser_aln.add_argument("--nometa", dest="nometa", action="store_true", default=False, help="Produce a gfa graph without node annotations, to ensure it's parseable by other programs.")
    parser_aln.set_defaults(func=align.align_cmd)
    
    parser_extract.add_argument('graph', nargs=1, help='gfa file specifying the graph from which the genome should be extracted.')
    parser_extract.add_argument('input', nargs='*', help='Name of the sample or path to be extracted from the graph.')
    parser_extract.add_argument("-t", dest="type", default="pathname", choices=["pathname","path"], help="Type of input, either pathname or comma-seperated sequence of node-ids.")
    parser_extract.add_argument("--width", dest="width", type=int, default=100 , help='Line width for fasta output.')
    parser_extract.set_defaults(func=extract.extract_cmd)
    
    parser_plot.add_argument('fastas', nargs='*', help='Two fasta files for which a mumplot should be generated.')
    parser_plot.add_argument("-m", dest="minlength", type=int, default=100, help="Minimum length of exact matches to vizualize.")
    parser_plot.add_argument("-i", dest="interactive", action="store_true", default=False, help="Wheter to produce interactive plots which allow zooming on the dotplot.")
    parser_plot.add_argument("-u", dest="uniq", action="store_true", default=True, help="Plot only maximal unique matches.")
    parser_plot.add_argument("-e", dest="uniq", action="store_false", default=True, help="Plot all maximal exact matches.")
    parser_plot.add_argument("--norc", dest="rc", action="store_false", default=True, help="Don't draw reverse complement matches.")
    parser_plot.add_argument("--maxn", dest="maxn", type=int, default=10000, help="Cap the number of MUMs to draw.")
    parser_plot.add_argument("--nogaps", dest="showgaps", action="store_false", default=True, help="Don't mark gapped sequence.")
    parser_plot.add_argument("-r", dest="region", default=None, help="Highlight interval (as \"<start>:<end>\") with respect to x-axis (first sequence).")
    parser_plot.set_defaults(func=plot.plot)
    
    parser_gplot.add_argument('graph', help='A graph representing the genomes of two or more samples.')
    parser_gplot.add_argument("-m", dest="minlength", type=int, default=1, help="Minimum length of exact matches to vizualize.")
    parser_gplot.add_argument("-x", default=None, help='Name of sample 1 (x-axis), when graph contains more than two samples, assignment is random.')
    parser_gplot.add_argument("-y", default=None, help='Name of sample 2 (y-axis), when graph contains more than two samples, assignment is random.')
    parser_gplot.add_argument("-i", dest="interactive", action="store_true", default=False, help="Wheter to produce interactive plots which allow zooming on the dotplot.")
    parser_gplot.add_argument("-r", dest="region", default=None, help="Highlight interval (as \"<start>:<end>\") with respect to x-axis (first sequence).")
    parser_gplot.set_defaults(func=plot.gplot)
    
    parser_comp.add_argument('graph', nargs=1, help='The graph to be reverse complemented.')
    parser_comp.set_defaults(func=comp.comp_cmd)
    
    parser_finish.add_argument('reference', help='Multi-fasta reference sequence.')
    parser_finish.add_argument('contigs', help='Multi-fasta draft assembly that contains contigs that are to be oriented and ordered with respect to the reference.')
    parser_finish.add_argument("-o", "--output", dest="output", help="Prefix of fasta file for the \'finished\' genome.")
    parser_finish.add_argument("-m", dest="minlength", type=int, default=20, help="Min length of maximal exact matches for considering.")
    parser_finish.add_argument("-i", dest="interactive", action="store_true", default=False, help="Output interactive plot.")
    parser_finish.add_argument("--nproc", dest="nproc", default=1, type=int, help="Use multiprocessing to do MUM extraction (max: 2 proc) and mapping (max: number of contigs) in parallel (increases mem usage!).")
    parser_finish.add_argument("--gcmodel", dest="gcmodel", choices=["sumofpairs","star-avg","star-med"], default="sumofpairs", help="Which gap-cost model to use for multi-alignment.")
    parser_finish.add_argument("--plot", dest="plot", action="store_true", default=False, help="Output mumplots for the \'finished\' chromosomes (depends on matplotlib).")
    parser_finish.add_argument("--outputgraph", dest="outputgraph", action="store_true", default=False, help="Output a graph based representation, in which nodes are contigs/chains and edges represent the reference based order.")
    parser_finish.add_argument("--allcontigs", dest="allcontigs", action="store_true", default=False, help="Output all contigs as separate paths through the graph (if --outputgraph is specified), otherwise only report contigs that contain structural variants.")
    parser_finish.add_argument("--filter", dest="filtermums", action="store_true", default=False, help="Reduce search space by filtering exact matches.")
    parser_finish.add_argument("--plotall", dest="plotall", action="store_true", default=False, help="Plot all matches, instead of only the chained matches.")
    parser_finish.add_argument("--split", dest="split", action="store_true", default=False, help="Split the \'finished\' genome by chromosome.")
    parser_finish.add_argument("--order", dest="order", default="contigs", choices=["contigs","chains"], help="Determine the order for either contigs or chains. With \'chains\' large scale rearrangements within contigs can be undone in order to obtain colinear genomes.")
    parser_finish.add_argument("--maxgapsize", dest="maxgapsize", type=int, default=1500, help="Maxgapsize between MUMs before breaking chain.")
    parser_finish.add_argument("--maxn", dest="maxn", type=int, default=100000000, help="Max number of MUMs to consider for chaining.")
    parser_finish.add_argument("--cutn", dest="cutn", type=int, default=1000, help="Cut contigs at N-stretches longer than this value, to force re-estimation of gapsizes (set to 0, to switch off).")
    parser_finish.add_argument("--minchainsum", dest="minchainsum", type=int, default=10000,help="Minimal sum of the length of the MUMs in a chain before its reported.")
    parser_finish.add_argument("--fixedgapsize", dest="fixedsize", action="store_true", default=False, help="Do not estimate gapsize based on reference, instead use fixed gapsizes of length that can be set with \'gapsize\'.")
    parser_finish.add_argument("--gapsize", dest="gapsize", type=int, default=100, help="Use this number of N's between adjacent (only in case of fixedgapsizes) or  partially overlapping contigs.")
    parser_finish.set_defaults(func=finish.finish)
    
    #parser_matches.add_argument('reference', help='Graph or sequence to which query/contigs should be assigned.')
    #parser_matches.add_argument('contigs', help='Graph or fasta that is to be reverse complemented with respect to the reference.')    
    #parser_matches.add_argument("-m", dest="minlength", type=int, default=20, help="Min length of maximal exact matches for considering.")
    #parser_matches.add_argument("-u", dest="uniq", action="store_true", default=True, help="Output only unique exact matches.")
    #parser_matches.add_argument("-e", dest="uniq", action="store_false", default=True, help="Output all exact matches.")
    #parser_matches.add_argument("--norc", dest="rc", action="store_false", default=True, help="Also output reverse complement matches.")
    #parser_matches.add_argument("--cache", dest="cache", default=False, action="store_true", help="When specified, it caches the text, suffix and lcp array to disk after construction.")
    #parser_matches.add_argument("--sa1", dest="sa1", default="", help="Specify a preconstructed suffix array for extracting matches between the two genomes in their current orientation.")
    #parser_matches.add_argument("--lcp1", dest="lcp1", default="", help="Specify a preconstructed lcp array for extracting matches between the two genomes in their current orientation.")
    #parser_matches.add_argument("--sa2", dest="sa2", default="", help="Specify a preconstructed suffix array for extracting matches between the two genomes in which the sequence of the contigs was reverse complemented.")
    #parser_matches.add_argument("--lcp2", dest="lcp2", default="", help="Specify a preconstructed lcp array for extracting matches between the two genomes in which the sequence of the contigs was reverse complemented.")    
    #parser_matches.set_defaults(func=matches.matches)
    
    parser_convert.add_argument('graphs', nargs='*', help='The gfa graph to convert to gml.')
    parser_convert.add_argument("-n", dest="minsamples", type=int, default=1, help="Only align nodes that occcur in this many samples.")
    parser_convert.add_argument("-x", dest="maxsamples", type=int, default=None, help="Only align nodes that have maximally this many samples.")
    parser_convert.add_argument("-s", dest="targetsample", type=str, default=None, help="Only align nodes in which this sample occurs.")
    parser_convert.add_argument("--gml-max", dest="hwm", default=4000, type=int, help="Max number of nodes per graph in gml output.")
    parser_convert.add_argument("--gfa",  action="store_true", dest="gfa", default=False, help="Rewrite gfa file.")
    parser_convert.add_argument("--partition",  action="store_true", dest="partition", default=False, help="Output graph as multiple subgraphs if possible.")
    parser_convert.set_defaults(func=convert.convert)
    
    parser_subgraph.add_argument('inputfiles', nargs='*', help='The gfa graph followed by node ids that make up the subgraph.')
    parser_subgraph.add_argument("-o", dest="outfile", type=str, default="subgraph", help="Prefix of the file to which subgraph will be written.")
    parser_subgraph.add_argument("--gml", dest="gml", action="store_true", default=False, help="Produce a gml graph instead of gfa.")
    parser_subgraph.set_defaults(func=subgraph.subgraph)
    
    parser_bubbles.add_argument("graph", nargs=1, help='Graph in gfa format from which bubbles are to be extracted.')
    parser_bubbles.add_argument("-e", dest="exportcomplex", action="store_true", default=False, help="Output complex bubble structures in a separate gfa file.")
    parser_bubbles.add_argument("-s", dest="separate", action="store_true", default=False, help="Write a seperate gfa file for each complex bubble structure.")
    parser_bubbles.add_argument("--gml", dest="gml", action="store_true", default=False, help="Output gml instead of gfa.")
    parser_bubbles.set_defaults(func=bubbles.bubbles_cmd)
    
    parser_variants.add_argument("graph", nargs=1, help='Graph in gfa format from which bubbles are to be extracted.')
    parser_variants.add_argument("-r", dest="reference", type=str, default=None, help="Name of the sequence that, if possible, should be used as a coordinate system or reference.")
    parser_variants.add_argument("--fasta", dest="fastaout", action="store_true", default=False, help="Output variant sequence in a fasta format.")
    parser_variants.add_argument("--minsize", dest="minsize", default=0, type=int, help="Only sequence variants larger than this.")
    parser_variants.set_defaults(func=bubbles.variants_cmd)
    
    parser_realign.add_argument("graph", nargs=1, help='Graph in gfa format for which a bubble should be realigned.') 
    parser_realign.add_argument("source", nargs='?', type=int, help='Source node.')
    parser_realign.add_argument("sink", nargs='?', type=int, help='Sink node.')
    parser_realign.add_argument("-m", dest="minlength", type=int, default=20, help="Min length of an exact match.")
    parser_realign.add_argument("-c", dest="minscore", type=int, default=0, help="Min score of an exact match, exact maches are scored by their length and penalized by the indel they create with respect to previously accepted exact matches.")
    parser_realign.add_argument("-n", dest="minn", type=int, default=2, help="Only align graph on exact matches that occur in at least this many samples.")
    parser_realign.add_argument("-o", dest="outfile", type=str, default=None, help="File to which realigned graph is to be written.")
    parser_realign.add_argument("--all", action="store_true", dest="all", default=False, help="Trigger realignment for all complex bubbles.")
    parser_realign.add_argument("--maxlen", dest="maxlen", type=int, default=10000000, help="Maximum length of the cumulative sum of all paths that run through the complex bubble.")
    parser_realign.add_argument("--maxsize", dest="maxsize", type=int, default=500, help="Maximum allowed number of nodes that are contained in a complex bubble.")
    parser_realign.add_argument("-e", dest="exp", type=int, default=1, help="Increase \'e\' to prefer shorter matches observed in more genomes, over larger matches in less genomes.")
    parser_realign.add_argument("--wp", dest="wpen", type=int, default=1, help="Multiply penalty for a MUM by this number in scoring scheme.")
    parser_realign.add_argument("--ws", dest="wscore", type=int, default=3, help="Multiply length of MUM by this number in scoring scheme.")
    parser_realign.add_argument("-p", dest="pcutoff", type=float, default=None, help="If, the probability of observing a MUM of the observed length by random change becomes larger than this cutoff the alignment is stopped.")
    parser_realign.set_defaults(func=realign.realign_bubble_cmd)
    
    parser_merge.add_argument("graphs", nargs='*', help='Graphs in gfa format that should be merged.')
    parser_merge.add_argument("-o", dest="outprefix", type=str, default=None, help="Prefix of the file to which merged graph is written.")
    parser_merge.set_defaults(func=merge.merge_cmd)
    
    parser_chain.add_argument('fastas', nargs='*', help='Fasta files specifying sequences to be (multi-)aligned into a graph.')
    parser_chain.add_argument("-o", "--output", dest="output", help="Prefix of the variant and alignment graph files to produce, default is \"sequence1_sequence2\"")
    parser_chain.add_argument("-m", dest="minlength", type=int, default=20, help="Min length of an exact (multi-)match to consider for chaining.")
    parser_chain.add_argument("-n", dest="minn", type=int, default=2, help="Only align graph on exact matches that occur in at least this many samples.")
    parser_chain.add_argument("-a", dest="maxmums", type=int, default=5000, help="Number of largest mums to use for chaining.")
    
    parser_chain.add_argument("-e", dest="exp", type=int, default=1, help="Increase \'e\' to prefer shorter matches observed in more genomes, over larger matches in less genomes.")
    parser_chain.add_argument("--wp", dest="wpen", type=int, default=1, help="Multiply penalty for a MUM by this number in scoring scheme.")
    parser_chain.add_argument("--ws", dest="wscore", type=int, default=3, help="Multiply length of MUM by this number in scoring scheme.")
    parser_chain.add_argument("--gcmodel", dest="gcmodel", choices=["sumofpairs","star-avg","star-med"], default="sumofpairs", help="Which gap-cost model to use for multi-alignment.")
    parser_chain.add_argument("--recurse", dest="recurse", action="store_true", default=False, help="Use recursive approach to chain gaps.")
    parser_chain.add_argument("--plot", dest="mumplot", action="store_true", default=False, help="Save a mumplot for the actual aligned chain of anchors (depends on matplotlib).")
    parser_chain.add_argument("-i", dest="interactive", action="store_true", default=False, help="Show an interactive visualisation of the mumplot (depends on matplotlib).")
    parser_chain.add_argument("--nometa", dest="nometa", action="store_true", default=False, help="Produce a gfa graph without node annotations, to ensure it's parseable by other programs.")
    parser_chain.set_defaults(func=chain.chain_cmd)
    
    parser_stats.add_argument('gfa', nargs=1, help='GFA file for which statistics should be calculated.')
    parser_stats.set_defaults(func=stats.stats_cmd)
    
    parser_split.add_argument('gfa', nargs=1, help='GFA file which has to  be split into a GFA file per connected component.')
    parser_split.set_defaults(func=split.split_cmd)
    
    args = parser.parse_args()
    
    logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p', level=args.loglevel)
    logging.logThreads = 0
    
    args.func(args)
