#!/usr/bin/env python
from intervaltree import Interval, IntervalTree
import networkx as nx
import argparse
import os
import sys
import time

#reveal imports
import schemes
import transform
import transformold
import plot
import utils
import merge
import subgraph
import comp
import extract
import convert
import rem
import refine
import bubbles
import matches
import chain
import stats
import split
import align
import unzip
import chop
import annotate


import logging
#add custom loglevel TRACE
logging.TRACE = 1
logging.addLevelName(logging.TRACE, "TRACE")
logging.logThreads = 0
logging.Logger.trace = lambda inst, msg, *args, **kwargs: inst.log(logging.TRACE, msg, *args, **kwargs)
logging.trace = lambda msg, *args, **kwargs: logging.log(logging.TRACE, msg, *args, **kwargs)

def main():
    desc="""
    Type 'reveal <positional_argument> --help' for help on a specific subcommand.\n
    Reveal constructs population reference graphs by aligning multiple whole 
    genomes using recursive exact matching.
    http://www.biorxiv.org/content/early/2015/07/17/022715.
    """
    
    parser = argparse.ArgumentParser(prog="reveal", usage="reveal -h for usage", description=desc, formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    global_parser = argparse.ArgumentParser(add_help=False) #parser for arguments that apply to all subcommands
    global_parser.add_argument("-l", "--log-level", type=int, dest="loglevel", default=20, help="Log level: 1=trace 10=debug 20=info 30=warn 40=error 50=fatal.")
    global_parser.add_argument("--64", dest="sa64", default=False, action="store_true", help="Use 64bit suffix array in the index.")
    
    subparsers = parser.add_subparsers()
    
    parser_aln = subparsers.add_parser('align',prog="reveal align", description="Output a bash-script that decribes the various steps that need to be performed to construct a graph-based multi-genome for the input (draft-)genomes.", formatter_class=argparse.ArgumentDefaultsHelpFormatter, parents=[global_parser])
    parser_aln.add_argument('reference', nargs=1, help='Reference assembly to which draft-asemblies should be organized (this may be a draft assembly as well, just used for organizing structural rearrangements).')
    parser_aln.add_argument('inputfiles', nargs='+', help='(Multi-)Fasta files that specify the draft-assemblies that are to be aligned (.fasta).')
    parser_aln.add_argument("-o", "--output", dest="output", default="prg", help="Prefix for the filename of the resulting graph.")
    parser_aln.add_argument("-m", dest="m", default=20, type=int, help="Min length of an anchor to constrain the alignment.")
    parser_aln.add_argument("-n", dest="n", default=None, type=int, help="Number of genomes for anchor to constrain the alignment.")
    parser_aln.add_argument("--chunksize", dest="chunksize", type=int, default=2, help="If order is sequential, use this many genomes per alignment.")
    parser_aln.add_argument("--nproc", dest="nproc", default=1, type=int, help="Number of processes to use for the individual steps.")
    parser_aln.add_argument("--order", dest="order", default="simultaneous", choices=["simultaneous","sequential","tree"], help="Order in which graph is constructed.")
    parser_aln.add_argument("--norefine", dest="refine", action="store_false", default=True, help="Do not use consistency based msa, just produce the anchor graph.")
    parser_aln.add_argument("--minconf", dest="minconf", type=float, default=90, choices=range(0,100), metavar="[0-99]", help="Use cutoff on confidence values during refinement of the graph.")
    parser_aln.add_argument("--nounzip", dest="unzip", action="store_false", default=True, help="Do not unzip bubbles before refining.")
    parser_aln.add_argument("--notransform", dest="transform", action="store_false", default=True, help="Do not use transform to account for structural events or draft assemblies. Assume colinear alignment of complete genomes.")
    parser_aln.add_argument("--novariants", dest="variants", action="store_false", default=True, help="Do not output bubbles as variants.")
    parser_aln.set_defaults(func=align.align)

    parser_rem = subparsers.add_parser('rem',prog="reveal rem", description="Use recursive exact matching to obtain a graph from multiple input genomes or other graphs.", formatter_class=argparse.ArgumentDefaultsHelpFormatter, parents=[global_parser])
    parser_rem.add_argument('inputfiles', nargs='+', help='Fasta or gfa files specifying either assembly/alignment graphs (.gfa) or sequences (.fasta).')
    parser_rem.add_argument("-o", "--output", dest="output", help="Prefix of the alignment graph.")
    parser_rem.add_argument("-t", "--threads", dest="threads", type=int, default=0, help = "The number of threads to use for the alignment.")
    parser_rem.add_argument("-m", dest="minlength", type=int, default=20, help="Min length of an exact match.")
    parser_rem.add_argument("-p", dest="pcutoff", type=float, default=1e-08, help="Use this significance threshold for exact matches, when -m=0")
    parser_rem.add_argument("-n", dest="minn", type=int, default=2, help="Only align graph on exact matches that occur in at least this many samples (if not set, equal to total number of genomes in the (sub)index).")
    parser_rem.add_argument("--gcmodel", dest="gcmodel", choices=["sumofpairs","star-avg","star-med"], default="sumofpairs", help="Which gap-cost model to use.")
    parser_rem.add_argument("--wp", dest="wpen", type=int, default=1, help="Weight of penalty during chaining.")
    parser_rem.add_argument("--ws", dest="wscore", type=int, default=1, help="Weight of score during chaining.")
    parser_rem.add_argument("--seedsize", dest="seedsize", type=int, default=10000, help="Skip recursion for chained mums larger than this size (when 0 don't seed).")
    parser_rem.add_argument("--maxmums", dest="maxmums", type=int, default=1000, help="Number of largest MUMs to consider for chaining per iteration (when 0 use all).")
    parser_rem.add_argument("--plot", dest="mumplot", action="store_true", default=False, help="Save a mumplot for the actual aligned chain of anchors (depends on matplotlib).")
    parser_rem.add_argument("-i", dest="interactive", action="store_true", default=False, help="Show an interactive visualisation of the mumplot (depends on matplotlib).")
    parser_rem.add_argument("--sa", dest="sa", default="", help="Specify a preconstructed suffix array to decouple suffix array construction.")
    parser_rem.add_argument("--lcp", dest="lcp", default="", help="Specify a preconstructed lcp array to decouple lcp array construction.")
    parser_rem.add_argument("--cache", dest="cache", default=False, action="store_true", help="When specified, it caches the suffix and lcp array to disk after construction.")
    parser_rem.add_argument("-g", dest="minsamples", type=int, default=1, help="Only index nodes that occur in this many samples or more.")
    parser_rem.add_argument("-x", dest="maxsamples", type=int, default=None, help="Only align nodes that have maximally this many samples.")
    parser_rem.add_argument("-r", dest="reference", type=str, default=None, help="Name of the sequence that should be used as a coordinate system or reference.")
    parser_rem.add_argument("-s", dest="targetsample", type=str, default=None, help="Only align nodes in which this sample occurs.")    
    parser_rem.add_argument("--gml", dest="gml", action="store_true", default=False, help="Produce a gml graph instead gfa.")
    parser_rem.add_argument("--gml-max", dest="hwm", default=4000, help="Max number of nodes per graph in gml output.")
    parser_rem.add_argument("--noupper", dest="toupper", action="store_false", default=True, help="Do not force uppercase for sequence, in case we want to prevent matches in repeatmasked sequence.")
    parser_rem.add_argument("--maxbubblesize", dest="maxsize", type=int, default=None, help="Apply recursion until largest allele within a bubble is smaller than this size.")
    parser_rem.add_argument("--nocontigs", dest="contigs", default=True, action="store_false", help="Don't treat multi-fasta files as contigs, use every sequence as a target.")
    parser_rem.add_argument("--notrim", dest="trim", default=True, action="store_false", help="Don't trim overlap between MUMs, thus more greedy positioning of indels.")
    parser_rem.set_defaults(func=rem.align_cmd)

    parser_unzip = subparsers.add_parser('unzip',prog="reveal unzip", description="Opens up bubbles to account for uncertainty of indel placement and edge-wander. Specify --source and --sink to unzip a specific bubble.", formatter_class=argparse.ArgumentDefaultsHelpFormatter, parents=[global_parser])
    parser_unzip.add_argument("graph", nargs=1, help='Graph in gfa format for which bubbles should be unzipped.')
    parser_unzip.add_argument("-u", dest="minunzip", type=int, default=0, help="Try to unzip all bubbles at least this many bases.")
    parser_unzip.add_argument("-o", "--output", dest="output", default=None, help="Prefix for the filename of the resulting graph.")
    parser_unzip.add_argument("--mindiff", dest="mindiff", default=1, type=int, help="Only unzip bubbles where the difference between the min- and max-allele size is larger than this many bp, by default 1, so don't unzip SNPs.")
    parser_unzip.add_argument("--maxdiff", dest="maxdiff", default=10000, type=int, help="Only unzip bubbles where the difference between the min- and max-allele size is smaller than this many bp.")
    parser_unzip.add_argument("--source", dest="source", type=int, default=None, help="Source for specific bubble.")
    parser_unzip.add_argument("--sink", dest="sink", type=int, default=None, help="Sink for specific bubble.")
    parser_unzip.set_defaults(func=unzip.unzip)

    parser_chop = subparsers.add_parser('chop',prog="reveal chop", description="Uses the chop algorithm to introduce overlap (of length k-1) onto the edges of the graph such that reads of length k can be mapped on to the graph.", formatter_class=argparse.ArgumentDefaultsHelpFormatter, parents=[global_parser])
    parser_chop.add_argument("graph", nargs=1, help='Graph in gfa format which has to be chopped.')
    parser_chop.add_argument("-k", dest="k", type=int, default=100, help="Max length of the reads that need to be mapped to the graph.")
    parser_chop.add_argument("-o", "--output", dest="output", default=None, help="Prefix for the resulting overlap graph and fasta file.")
    parser_chop.add_argument("--noextend", dest="extend", default=True, action="store_false", help="Do not add prefix/suffix, just apply duplicate and contract.")
    parser_chop.add_argument("--fasta", dest="fasta", default=False, action="store_true", help="Write node sequence to a fasta file, for read mapping.")
    parser_chop.add_argument("--width", dest="lw", type=int, default=100, help="Line width for fasta output.")
    parser_chop.add_argument("--check", dest="check", default=False, action="store_true", help="Check if all k-length substrings of the input haplotype are covered in the flat representation.")
    parser_chop.set_defaults(func=chop.chop_cmd)

    parser_refine = subparsers.add_parser('refine', prog="reveal refine", description="Refine bubbles in the graph by multiple sequence alignment.", formatter_class=argparse.ArgumentDefaultsHelpFormatter, parents=[global_parser])
    parser_refine.add_argument("graph", nargs=1, help='Graph in gfa format for which bubbles should be realigned.')
    parser_refine.add_argument("source", nargs='?', default=None, type=int, help='Source node.')
    parser_refine.add_argument("sink", nargs='?', default=None, type=int, help='Sink node.')
    parser_refine.add_argument("--nproc", dest="nproc", default=1, type=int, help="Use multiprocessing to realign bubbles.")
    parser_refine.add_argument("--chunksize", dest="chunksize", type=int, default=10, help="Process in parallel refined bubbles in chunks of this size.")
    parser_refine.add_argument("--method", dest="method", choices=["reveal_probcons","reveal_rem","muscle","probcons","msaprobs","pecan"], default="reveal_probcons", help="Use external multiple sequence aligner for the alignment of bubbles (expects methods to be accessible through the $PATH variable, reveal_* methods use internal memory, other methods uese external memory")
    parser_refine.add_argument("--params", dest="parameters", default="", help="Add this value when calling the external methods.")
    parser_refine.add_argument("-o", dest="outfile", type=str, default=None, help="File to which realigned graph is to be written.")
    parser_refine.add_argument("--all", action="store_true", dest="all", default=False, help="Trigger realignment for all bubbles.")
    parser_refine.add_argument("--complex", action="store_true", dest="complex", default=False, help="Trigger realignment for all complex bubbles.")
    parser_refine.add_argument("--nogaps", action="store_true", dest="nogaps", default=False, help="Skip realignment for bubbles that span gaps.")
    parser_refine.add_argument("--simple", action="store_true", dest="simple", default=False, help="Trigger realignment for all simple bubbles.")
    parser_refine.add_argument("--minsize", dest="minsize", type=int, default=0, help="Only realign bubbles if the smallest allele contains at least this many bases.")
    parser_refine.add_argument("--maxsize", dest="maxsize", type=int, default=10000, help="Only realign bubbles if the largest allele contains less than this many bases.")
    # parser_refine.add_argument("--minmaxsize", dest="minmaxsize", type=int, default=2, help="Only realign bubbles if the largest allele contains more than this many bases.")
    parser_refine.add_argument("--mindiff", dest="mindiff", default=1, type=int, help="Only refine variants where the difference between the min- and max-allele size is larger than this many bp.")
    parser_refine.add_argument("--maxdiff", dest="maxdiff", default=None, type=int, help="Only refine variants where the difference between the min- and max-allele size is smaller than this many bp.")
    parser_refine.add_argument("--maxcumsize", dest="maxcumsize", type=int, default=None, help="Maximum length of the cumulative sum of all paths that run through the bubble.")
    parser_refine.add_argument("--mincumsize", dest="mincumsize", type=int, default=0, help="Minimum length of the cumulative sum of all paths that run through the bubble.")
    parser_refine.add_argument("--minconf", dest="minconf", type=float, default=0, choices=range(0,101), metavar="[0-100]", help="Use cutoff on confidence values from the MSA in graph construction ().")
    parser_refine.add_argument("--uniqueonly", dest="uniqueonly", default=False, action="store_true", help="Only consider unique haplotypes for multiple sequence alignment.")
    parser_refine.add_argument("-c","--consistency", dest="constrans", type=int, default=2, help="Number of consistency transformations to apply before alignment (only applies to reveal_probcons).")
    parser_refine.add_argument("-r","--iterative-refinement", dest="nrefinements", type=int, default=100, help="Number of iterative refinements to apply after alignment (only applies to reveal_probcons).")
    parser_refine.add_argument("--no-gap-consistency", dest="consgap", action="store_false", default=True, help="Don't consider gaps in consistency transform (only applies to reveal_probcons).")
    parser_refine.add_argument("-m", dest="minlength", type=int, default=20, help="Min length of an exact match (only applies when method is 'reveal_rem').")
    parser_refine.add_argument("-n", dest="minn", type=int, default=2, help="Only align graph on exact matches that occur in at least this many samples (only applies when method is 'reveal_rem').")
    parser_refine.add_argument("--gcmodel", dest="gcmodel", choices=["sumofpairs","star-avg","star-med"], default="sumofpairs", help="Which gap-cost model to use for multi-alignment (only applies when method is 'reveal_rem').")
    parser_refine.add_argument("--wp", dest="wpen", type=int, default=1, help="Multiply penalty for a MUM by this number in scoring scheme (only applies when method is 'reveal_rem').")
    parser_refine.add_argument("--ws", dest="wscore", type=int, default=1, help="Multiply length of MUM by this number in scoring scheme (only applies when method is 'reveal_rem').")
    parser_refine.add_argument("--seedsize", dest="seedsize", type=int, default=10000, help="Skip recursion for chained mums larger than this size (when 0 don't seed) (only applies when method is 'reveal_rem').")
    parser_refine.add_argument("--maxmums", dest="maxmums", type=int, default=1000, help="Number of largest MUMs to consider for chaining (when 0 use all) (only applies when method is 'reveal_rem').")
    parser_refine.set_defaults(func=refine.refine_bubble_cmd)
    parser_realign = subparsers.add_parser('realign', prog="reveal realign", parents=[parser_refine], add_help=False)

    parser_extract = subparsers.add_parser('extract', prog="reveal extract", description="Extract the input sequence from a graph.", formatter_class=argparse.ArgumentDefaultsHelpFormatter, parents=[global_parser])
    parser_extract.add_argument('graph', nargs=1, help='gfa file specifying the graph from which the genome should be extracted.')
    parser_extract.add_argument('input', nargs='*', help='Name of the sample or path to be extracted from the graph.')
    parser_extract.add_argument("-t", dest="type", default="pathname", choices=["pathname","path"], help="Type of input, either pathname or comma-seperated sequence of node-ids.")
    parser_extract.add_argument("--width", dest="width", type=int, default=100 , help='Line width for fasta output.')
    parser_extract.add_argument("--all", dest="all", default=False, action="store_true", help="Extract all paths from the graph and output as fasta.")
    parser_extract.add_argument("--nocycles",  action="store_true", dest="nocycles", default=False, help="Parse only the directed acyclic layout of the graph, so ignore strucural rearrangements (cycles) in the graph.")
    parser_extract.set_defaults(func=extract.extract_cmd)
    
    parser_plot = subparsers.add_parser('plot', prog="reveal plot", description="Generate a mumplot that shows all mums between two fasta files.", formatter_class=argparse.ArgumentDefaultsHelpFormatter, parents=[global_parser])
    parser_plot.add_argument('fastas', nargs='*', help='Two fasta files for which a mumplot should be generated.')
    parser_plot.add_argument("-m", dest="minlength", type=int, default=20, help="Minimum length of exact matches to vizualize.")
    parser_plot.add_argument("-i", dest="interactive", action="store_true", default=False, help="Wheter to produce interactive plots which allow zooming on the dotplot.")
    parser_plot.add_argument("--endpoints", dest="endpoints", action="store_true", default=False, help="Mark mum start/end points.")
    parser_plot.add_argument("--norc", dest="rc", action="store_false", default=True, help="Don't draw reverse complement matches.")
    parser_plot.add_argument("--maxmums", dest="maxmums", type=int, default=10000, help="Cap the number of MUMs to draw.")
    parser_plot.add_argument("--nogaps", dest="showgaps", action="store_false", default=True, help="Don't mark gapped sequence.")
    parser_plot.add_argument("--extension", dest="extension", choices=["png","pdf","ps","eps","svg"], default="png", help="How to save the plot.")
    parser_plot.add_argument("-r","--xr", dest="xregion", default=None, help="Highlight and zoom on intervals (encoded as \"<start1>-<end1>,<start2>-<end2>\" etc.) with respect to x-axis (first sequence).")
    parser_plot.add_argument("--yr", dest="yregion", default=None, help="Highlight and zoom on intervals (encoded as \"<start1>-<end1>,<start2>-<end2>\" etc.) with respect to y-axis (second sequence).")
    parser_plot.add_argument("--flanksize", dest="flanksize", default=None, help="In case of (a) specified region include this many bases of flanking sequence (encode as \"<flanksize_region1>,<flanksize_region2>\" etc.).")
    parser_plot.set_defaults(func=plot.plot)
    
    parser_gplot = subparsers.add_parser('gplot', prog="reveal gplot", description="Generate a plot that represents the alignment of two samples in a graph.", formatter_class=argparse.ArgumentDefaultsHelpFormatter, parents=[global_parser])
    parser_gplot.add_argument('graph', help='A graph representing the genomes of two or more samples.')
    parser_gplot.add_argument("-m", dest="minlength", type=int, default=1, help="Minimum length of exact matches to vizualize.")
    parser_gplot.add_argument("-x", default=None, help='Name of sample 1 (x-axis), when graph contains more than two samples, assignment is random.')
    parser_gplot.add_argument("-y", default=None, help='Name of sample 2 (y-axis), when graph contains more than two samples, assignment is random.')
    parser_gplot.add_argument("-i", dest="interactive", action="store_true", default=False, help="Wheter to produce interactive plots which allow zooming on the dotplot.")
    parser_gplot.add_argument("-r", dest="region", default=None, help="Highlight interval (as \"<start>:<end>\") with respect to x-axis (first sequence).")
    parser_gplot.set_defaults(func=plot.gplot)
    
    parser_comp = subparsers.add_parser('comp', prog="reveal comp", description="Reverse complement the graph.", formatter_class=argparse.ArgumentDefaultsHelpFormatter, parents=[global_parser])
    parser_comp.add_argument('graph', nargs=1, help='The graph to be reverse complemented.')
    parser_comp.set_defaults(func=comp.comp_cmd)
    

    # #the old transform method
    # parser_transform = subparsers.add_parser('transform', prog="reveal transform", description="Transform a draft assembly into a graph that encodes the structural order of assembled segments with respect to a finished reference assembly.", formatter_class=argparse.ArgumentDefaultsHelpFormatter, parents=[global_parser])
    # parser_transform.add_argument('reference', help='(Multi-)fasta reference sequence.')
    # parser_transform.add_argument('contigs', help='(Multi-)fasta draft assembly that contains contigs that are to be oriented and ordered with respect to the reference.')
    # parser_transform.add_argument("-o", "--output", dest="output", help="Prefix of fasta file for the \'finished\' genome.")
    # parser_transform.add_argument("-m", dest="minlength", type=int, default=20, help="Min length of maximal exact matches for considering (if set to 0, try to extract all MUMs).")
    # parser_transform.add_argument("-i", dest="interactive", action="store_true", default=False, help="Output interactive plot.")
    # parser_transform.add_argument("--nproc", dest="nproc", default=1, type=int, help="Use multiprocessing to do MUM extraction (max: 2 proc) and mapping (max: number of contigs) in parallel (increases mem usage!).")
    # parser_transform.add_argument("--gcmodel", dest="gcmodel", choices=["sumofpairs","star-avg","star-med"], default="sumofpairs", help="Which gap-cost model to use for multi-alignment.")
    # parser_transform.add_argument("--plot", dest="plot", action="store_true", default=False, help="Output mumplots for the \'finished\' chromosomes (depends on matplotlib).")
    # parser_transform.add_argument("--graph", dest="outputtype", choices=["graph","fasta"], default="graph", help="Output a graph or fasta representation of the transformed genome.")
    # parser_transform.add_argument("--structonly", dest="allcontigs", action="store_false", default=True, help="Only output paths for contigs that contain structural rearrangements.")
    # # parser_transform.add_argument("--filter", dest="filtermums", action="store_true", default=False, help="Reduce search space by filtering exact matches.")
    # parser_transform.add_argument("--plotall", dest="plotall", action="store_true", default=False, help="Plot all matches, instead of only the chained matches.")
    # parser_transform.add_argument("--split", dest="split", action="store_true", default=False, help="Split the \'finished\' genome by chromosome.")
    # parser_transform.add_argument("--order", dest="order", default="chains", choices=["contigs","chains"], help="Determine the order for either contigs or chains.")
    # parser_transform.add_argument("--mineventsize", dest="mineventsize", type=int, default=200, help="Maximal distance between clusters/mums for chaining.")
    # parser_transform.add_argument("--minchainsum", dest="minchainsum", type=int, default=10000, help="Minimal sum of the length of the MUMs in a chain before its considered.")
    # parser_transform.add_argument("--maxmums", dest="maxmums", type=int, default=0, help="Max number of MUMs to consider for chaining (when 0 use all).")
    # parser_transform.add_argument("--cutn", dest="cutn", type=int, default=1000, help="Cut contigs at N-stretches longer than this value, to force re-estimation of gapsizes (set to 0, to switch off).")
    # parser_transform.add_argument("--fixedgapsize", dest="fixedsize", action="store_true", default=False, help="Do not estimate gapsize based on reference, instead use fixed gapsizes of length that can be set with \'gapsize\'.")
    # parser_transform.add_argument("--gapsize", dest="gapsize", type=int, default=100, help="Use this number of N's between adjacent (only in case of fixedgapsizes) or  partially overlapping contigs.")
    # parser_transform.add_argument("--maxdist", dest="maxdist", type=int, default=90, help="Max space between adjacent MUMs in a cluster.")
    # parser_transform.add_argument("--mincluster", dest="mincluster", type=int, default=65, help="Max space between adjacent MUMs in a cluster.")
    # parser_transform.add_argument("--extiter", dest="extiter", type=int, default=3, help="Number of extension iterations using locally unique MUMs.")
    # parser_transform.add_argument("--maxextend", dest="maxextend", type=int, default=1000, help="Size of the region to try to inspect for locally unique MUMs.")
    # parser_transform.add_argument("--ml", dest="minlocallength", type=int, default=20, help="Min size of locally unique mums.")
    # parser_transform.add_argument("--unmapped", dest="outputunmapped", action="store_true", default=False, help="Output a unmappable sequence to a separate fasta file.")
    # parser_transform.set_defaults(func=transform.transform)

    
    parser_transform = subparsers.add_parser('transform', prog="reveal transform", description="Transform a draft assembly into a graph that encodes the structural order of assembled segments with respect to a finished reference assembly.", formatter_class=argparse.ArgumentDefaultsHelpFormatter, parents=[global_parser])

    parser_transform.add_argument('reference', nargs=1, help='(Multi-)fasta reference sequence.')
    parser_transform.add_argument('contigs', nargs='+', help='(Multi-)fasta draft assembly that contains contigs that are to be oriented and ordered with respect to the reference.')

    parser_transform.add_argument("--cutn", dest="cutn", type=int, default=0, help="Cut contigs at N-stretches longer than this value (default is 0, off).")
    parser_transform.add_argument("-o", "--output", dest="output", help="Prefix of gfa file for the \'transformed\' genome.")
    parser_transform.add_argument("-c", dest="minctglength", type=int, default=10000, help="Skip transform for contigs short than this length.")
    parser_transform.add_argument("-m", dest="minlength", type=int, default=20, help="Min length of maximal exact matches for considering (if set to 0, try to extract all MUMs).")
    parser_transform.add_argument("-i", dest="interactive", action="store_true", default=False, help="Output interactive plot.")
    parser_transform.add_argument("--plot", dest="plot", action="store_true", default=False, help="Output mumplots for the \'finished\' chromosomes (depends on matplotlib).")

    parser_transform.add_argument("--rc", dest="rearrangecost", default=10000, type=int, help="Cost for chaining translocated segments.")
    parser_transform.add_argument("--ic", dest="inversioncost", default=5, type=int, help="Cost for chaining inverted segments.")
    
    parser_transform.add_argument("--alfa", dest="alfa", default=2, type=float, help="Weight for match (+).")
    parser_transform.add_argument("--lambda", dest="_lambda", default=3, type=float, help="Weight for indel penalty (-).")
    parser_transform.add_argument("--eps", dest="eps", default=2, type=float, help="Weight for substitution penalty (-).")
    parser_transform.add_argument("--gapopen", dest="gapopen", default=1, type=float, help="Fixed penalty for adding an achor to the chain (-).")
    
    parser_transform.add_argument("--nocluster", dest="cluster", action="store_false", default=True, help="Don't cluster MUMs by diagonals.")
    parser_transform.add_argument("--maxdist", dest="maxdist", type=int, default=30, help="Max space between adjacent MUMs (on the same diagonal) in a cluster.")
    parser_transform.add_argument("--mincluster", dest="mincluster", type=int, default=50, help="Minimal size (sum of mums) of a cluster.")
    parser_transform.add_argument("--minchainsum", dest="minchainsum", type=int, default=50, help="Minimal sum of the length of the MUMs in a chain before its considered.")

    parser_transform.add_argument("--noopt", dest="optimise", action="store_false", default=True, help="Don't perform naive optimisation of the glocal chain.")
    parser_transform.add_argument("--heap", dest="useheap", action="store_true", default=False, help="Use a priority queue to compute an optimal chain.")
    parser_transform.add_argument("--lastn", dest="lastn", type=int, default=50, help="Backtrack at least this many anchors while chaining fragments.")
    parser_transform.add_argument("--lastbp", dest="lastbp", type=int, default=20000, help="Backtrack at least this many bp while chaining fragments.")

    parser_transform.add_argument("--greedy", dest="greedy", action="store_true", default=False, help="Assign overlap between anchors in a greedy manner. Large anchors become larger.")
    parser_transform.add_argument("--outputbed", dest="outputbed", action="store_true", default=True, help="Produce a bed file that stores the rearrangement breakpoints on the reference assembly.")

    parser_transform.set_defaults(func=transform.transform_cmd)







    parser_finish = subparsers.add_parser('finish', prog="reveal finish", description="Order and orient the contigs of a draft assembly with respect to a finished reference assembly (essentially the same as transform but with different default parameters).", formatter_class=argparse.ArgumentDefaultsHelpFormatter, parents=[global_parser])
    parser_finish.add_argument('reference', help='(Multi-)fasta reference sequence.')
    parser_finish.add_argument('contigs', help='(Multi-)fasta draft assembly that contains contigs that are to be oriented and ordered with respect to the reference.')
    parser_finish.add_argument("-o", "--output", dest="output", help="Prefix of fasta file for the \'finished\' genome.")
    parser_finish.add_argument("-m", dest="minlength", type=int, default=15, help="Min length of maximal exact matches for considering (if not set, use the set of largest MUMs for which the genome wide coverage is below 1).")
    parser_finish.add_argument("-i", dest="interactive", action="store_true", default=False, help="Output interactive plot.")
    
    parser_finish.add_argument("--nproc", dest="nproc", default=1, type=int, help="Use multiprocessing to do MUM extraction (max: 2 proc) and mapping (max: number of contigs) in parallel (increases mem usage!).")
    parser_finish.add_argument("--gcmodel", dest="gcmodel", choices=["sumofpairs","star-avg","star-med"], default="sumofpairs", help="Which gap-cost model to use for multi-alignment.")
    parser_finish.add_argument("--plot", dest="plot", action="store_true", default=False, help="Output mumplots for the \'finished\' chromosomes (depends on matplotlib).")
    parser_finish.add_argument("--graph", dest="outputtype", choices=["graph","fasta"], default="fasta", help="Output a graph or fasta representation of the transformed genome.")
    parser_finish.add_argument("--allcontigs", dest="allcontigs", action="store_true", default=True, help="Output all contigs as separate paths through the graph, otherwise only report contigs that contain structural variants.")
    parser_finish.add_argument("--filter", dest="filtermums", action="store_true", default=False, help="Reduce search space by filtering exact matches.")
    parser_finish.add_argument("--plotall", dest="plotall", action="store_true", default=False, help="Plot all matches, instead of only the chained matches.")
    parser_finish.add_argument("--split", dest="split", action="store_true", default=False, help="Split the \'finished\' genome by chromosome.")
    parser_finish.add_argument("--order", dest="order", default="contigs", choices=["contigs","chains"], help="Determine the order for either contigs or chains.")
    parser_finish.add_argument("--mineventsize", dest="mineventsize", type=int, default=1500, help="Maximal distance between clusters/mums for chaining.")
    parser_finish.add_argument("--minchainsum", dest="minchainsum", type=int, default=1000, help="Minimal sum of the length of the MUMs in a chain before its considered.")
    parser_finish.add_argument("--maxmums", dest="maxmums", type=int, default=0, help="Max number of MUMs to consider for chaining (when 0 use all).")
    parser_finish.add_argument("--cutn", dest="cutn", type=int, default=1000, help="Cut contigs at N-stretches longer than this value, to force re-estimation of gapsizes (set to 0, to switch off).")
    parser_finish.add_argument("--fixedgapsize", dest="fixedsize", action="store_true", default=False, help="Do not estimate gapsize based on reference, instead use fixed gapsizes of length that can be set with \'gapsize\'.")
    parser_finish.add_argument("--gapsize", dest="gapsize", type=int, default=100, help="Use this number of N's between adjacent (only in case of fixedgapsizes) or  partially overlapping contigs.")
    
    parser_finish.add_argument("--maxdist", dest="maxdist", type=int, default=90, help="Max space between adjacent MUMs in a cluster.")
    parser_finish.add_argument("--mincluster", dest="mincluster", type=int, default=20, help="Max space between adjacent MUMs in a cluster.")
    
    parser_finish.add_argument("--extiter", dest="extiter", type=int, default=3, help="Number of iterations of alignment extension.")
    parser_finish.add_argument("--maxextend", dest="maxextend", type=int, default=200, help="Size of the region to try to inspect for locally unique MUMs.")
    parser_finish.add_argument("--ml", dest="minlocallength", type=int, default=20, help="Min size of locally unique mums.")
    
    parser_finish.add_argument("--nounmapped", dest="outputunmapped", action="store_false", default=True, help="Do not output unmappable sequence to a separate fasta file.")
    parser_finish.set_defaults(func=transformold.transform)




    parser_convert = subparsers.add_parser('convert', prog="reveal convert", description="Convert gfa graph to gml.", formatter_class=argparse.ArgumentDefaultsHelpFormatter, parents=[global_parser])
    parser_convert.add_argument('graphs', nargs='*', help='The gfa graph to convert to gml.')
    parser_convert.add_argument("-n", dest="minsamples", type=int, default=1, help="Only align nodes that occcur in this many samples.")
    parser_convert.add_argument("-x", dest="maxsamples", type=int, default=None, help="Only align nodes that have maximally this many samples.")
    parser_convert.add_argument("-s", dest="targetsample", type=str, default=None, help="Only align nodes in which this sample occurs.")
    parser_convert.add_argument("--gml-max", dest="hwm", default=4000, type=int, help="Max number of nodes per graph in gml output.")
    # parser_convert.add_argument("--gfa",  action="store_true", dest="gfa", default=False, help="Rewrite gfa file.")
    parser_convert.add_argument("--partition",  action="store_true", dest="partition", default=False, help="Output graph as multiple subgraphs if possible.")
    parser_convert.add_argument("--nocycles",  action="store_true", dest="nocycles", default=False, help="Do not allow rearrangements (cycles) in graph.")
    parser_convert.add_argument("--to", dest="type", default="gml", choices=['gml','gfa','maf'], help="Filetype to convert to.")
    parser_convert.add_argument("--aligned", dest="aligned", default=False, action="store_true", help="Whether multi fasta file is aligned.")
    parser_convert.set_defaults(func=convert.convert)
    
    parser_subgraph = subparsers.add_parser('subgraph', prog="reveal subgraph", description="Extract subgraph from gfa by specified node ids.", formatter_class=argparse.ArgumentDefaultsHelpFormatter, parents=[global_parser])
    parser_subgraph.add_argument('inputfiles', nargs='*', help="The gfa graph followed by a subgraph definition. Either comma-separated node ids (e.g. \"node1,node2,node3,...\"), topological range (all nodes between e.g. \"<node1>-<node2>\"), interval-based (e.g. \"chr4:34000-35000\").")
    parser_subgraph.add_argument("-o", dest="outfile", type=str, default="subgraph", help="Prefix of the file to which subgraph will be written.")
    parser_subgraph.add_argument("--gml", dest="gml", action="store_true", default=False, help="Produce a gml graph instead of gfa.")
    parser_subgraph.set_defaults(func=subgraph.subgraph)
    
    parser_bubbles = subparsers.add_parser('bubbles', prog="reveal bubbles", description="Extract all bubbles from the graph.", formatter_class=argparse.ArgumentDefaultsHelpFormatter, parents=[global_parser])
    parser_bubbles.add_argument("graph", nargs=1, help='Graph in gfa format from which bubbles are to be extracted.')
    parser_bubbles.add_argument("-e", dest="exportcomplex", action="store_true", default=False, help="Output complex bubble structures in a separate gfa file.")
    parser_bubbles.add_argument("-s", dest="separate", action="store_true", default=False, help="Write a seperate gfa file for each complex bubble structure.")
    parser_bubbles.add_argument("--gml", dest="gml", action="store_true", default=False, help="Output gml instead of gfa.")
    parser_bubbles.set_defaults(func=bubbles.bubbles_cmd)
    
    parser_variants = subparsers.add_parser('variants', prog="reveal variants", description="Extract variant calls from the graph.", formatter_class=argparse.ArgumentDefaultsHelpFormatter, parents=[global_parser])
    parser_variants.add_argument("graph", nargs=1, help='Graph in gfa format from which bubbles are to be extracted.')
    parser_variants.add_argument("-r", dest="reference", type=str, default=None, help="Name of the sequence that, if possible, should be used as a coordinate system or reference.")
    parser_variants.add_argument("--fasta", dest="fastaout", action="store_true", default=False, help="Output variant sequence in a fasta format.")
    parser_variants.add_argument("--bed", dest="bedout", action="store_true", default=False, help="Output position of variants in bed format.")
    parser_variants.add_argument("--vcf", dest="vcfout", action="store_true", default=False, help="Output variants in vcf format.")
    parser_variants.add_argument("--split", dest="split", action="store_true", default=False, help="Output a multi-fasta file per variant.")
    parser_variants.add_argument("--minsize", dest="minsize", default=0, type=int, help="Only output variants where max-allele size is larger than this many bp.")
    parser_variants.add_argument("--mindiff", dest="mindiff", default=0, type=int, help="Only output variants where the difference between the min- and max-allele size is larger than this many bp.")
    parser_variants.add_argument("--maxdiff", dest="maxdiff", default=None, type=int, help="Only output variants where the difference between the min- and max-allele size is smaller than this many bp.")
    parser_variants.add_argument("--minflank", dest="minflank", default=0, type=int, help="Only output variants with an exact matching flanking sequence of at least this length.")
    parser_variants.add_argument("--type", dest="type", default="all", choices=["all","snv","indel","multi-allelic","region","complex","undefined"], help="Only output variants of this type.")
    parser_variants.add_argument("--nogaps", dest="nogaps", default=False, action="store_true", help="Don't output variants that are caused by gaps (contain the N character).")
    parser_variants.add_argument("--refonly", dest="refonly", default=False, action="store_true", help="Don't output variants that are not positionable on the specified reference.")
    parser_variants.set_defaults(func=bubbles.variants_cmd)
    
    parser_rearrangements = subparsers.add_parser('rearrangements', prog="reveal rearrangements", description="Report on edges in the graph that describe rearrangements.", formatter_class=argparse.ArgumentDefaultsHelpFormatter, parents=[global_parser])
    parser_rearrangements.add_argument("-r", dest="reference", type=str, default=None, help="Name of the sequence that, should be used as a coordinate system or reference.")
    parser_rearrangements.add_argument("graph", nargs=1, help='Graph in gfa format for rearrangement edges are reported.')
    parser_rearrangements.set_defaults(func=bubbles.rearrangements_cmd)

    parser_annotate = subparsers.add_parser('annotate', prog="reveal annotate", description="Add annotations to variants in a vcf file using trf and repeatmasker.", formatter_class=argparse.ArgumentDefaultsHelpFormatter, parents=[global_parser])
    parser_annotate.add_argument("--species", dest="species", type=str, default="human", help="Which \'species\' flag to pass on to repeatmasker.")
    parser_annotate.add_argument("--nproc", dest="repmproc", type=int, default=1, help="How many processes repeatmasker should use for annotation (-pa flag).")
    parser_annotate.add_argument("--mindiff", dest="mindiff", default=50, type=int, help="Only annotate variants where the difference between the min- and max-allele size is larger or equal to this many bp.")
    parser_annotate.add_argument("--maxdiff", dest="maxdiff", default=100000, type=int, help="Only annotate variants where the difference between the min- and max-allele size is smaller than this many bp.")
    parser_annotate.add_argument("--notrf", dest="trf", default=True, action="store_false", help="Skip the tandem repeat finder for annotation.")
    parser_annotate.add_argument("--norepm", dest="repm", default=True, action="store_false", help="Skip the RepeatMasker for annotation.")
    parser_annotate.add_argument("vcffile", help='Variants from a graph in the vcf file format')
    parser_annotate.set_defaults(func=annotate.annotate)

    parser_merge = subparsers.add_parser('merge', prog="reveal merge", description="Combine multiple gfa graphs into a single gfa graph.", formatter_class=argparse.ArgumentDefaultsHelpFormatter, parents=[global_parser])
    parser_merge.add_argument("graphs", nargs='*', help='Graphs in gfa format that should be merged.')
    parser_merge.add_argument("-o", dest="outprefix", type=str, default=None, help="Prefix of the file to which merged graph is written.")
    parser_merge.set_defaults(func=merge.merge_cmd)

    parser_chain = subparsers.add_parser('chain', prog="reveal chain", description="Use default chaining scheme to construct GFA graph based on a global multi-alignment of all input genomes.", formatter_class=argparse.ArgumentDefaultsHelpFormatter, parents=[global_parser])    
    parser_chain.add_argument('fastas', nargs='*', help='Fasta files specifying sequences to be (multi-)aligned into a graph.')
    parser_chain.add_argument("-o", "--output", dest="output", help="Prefix of the variant and alignment graph files to produce, default is \"sequence1_sequence2\"")
    parser_chain.add_argument("-m", dest="minlength", type=int, default=20, help="Min length of an exact (multi-)match to consider for chaining.")
    parser_chain.add_argument("-n", dest="minn", type=int, default=2, help="Only align graph on exact matches that occur in at least this many samples.")
    parser_chain.add_argument("-a", dest="maxmums", type=int, default=0, help="Number of largest mums to use for chaining (when 0 use all).")
    parser_chain.add_argument("--wp", dest="wpen", type=int, default=1, help="Multiply penalty for a MUM by this number in scoring scheme.")
    parser_chain.add_argument("--ws", dest="wscore", type=int, default=1, help="Multiply length of MUM by this number in scoring scheme.")
    parser_chain.add_argument("--gcmodel", dest="gcmodel", choices=["sumofpairs","star-avg","star-med"], default="sumofpairs", help="Which gap-cost model to use for multi-alignment.")
    parser_chain.add_argument("--recurse", dest="recurse", action="store_true", default=False, help="Use recursive approach to chain gaps.")
    parser_chain.add_argument("--plot", dest="mumplot", action="store_true", default=False, help="Save a mumplot for the actual aligned chain of anchors (depends on matplotlib).")
    parser_chain.add_argument("-i", dest="interactive", action="store_true", default=False, help="Show an interactive visualisation of the mumplot (depends on matplotlib).")
    parser_chain.add_argument("--nometa", dest="nometa", action="store_true", default=False, help="Produce a gfa graph without node annotations, to ensure it's parseable by other programs.")
    parser_chain.set_defaults(func=chain.chain_cmd)
    
    parser_stats = subparsers.add_parser('stats', prog="reveal stats", description="Output statistics (number of node, edges, genomes etc.) for a graph.", formatter_class=argparse.ArgumentDefaultsHelpFormatter, parents=[global_parser])
    parser_stats.add_argument('gfa', nargs=1, help='GFA file for which statistics should be calculated.')
    parser_stats.set_defaults(func=stats.stats_cmd)

    parser_split = subparsers.add_parser('split', prog="reveal split", description="Split a graph file into a connected component per file.", formatter_class=argparse.ArgumentDefaultsHelpFormatter, parents=[global_parser])
    parser_split.add_argument('gfa', nargs=1, help='GFA file which has to  be split into a GFA file per connected component.')
    parser_split.add_argument("--nocycles",  action="store_true", dest="nocycles", default=False, help="Parse only the directed acyclic layout of the graph, so ignore strucural rearrangements (cycles) in the graph.")
    parser_split.set_defaults(func=split.split_cmd)
    
    args = parser.parse_args()
    
    logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p', level=args.loglevel)

    args.func(args)
