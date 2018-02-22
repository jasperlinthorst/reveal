# REVEAL

REVEAL (REcursiVe Exact-matching ALigner) can be used to (multi) align whole genomes.

Preprint can be found here: [Scalable multi whole-genome alignment using recursive exact matching](http://www.biorxiv.org/content/early/2015/07/17/022715)

## LATEST CHANGES

- New subcommands: rem, transform, refine and unzip (to substitute align, finish and realign)
- Package ProbCons code for refining bubbles

## INSTALL

REVEAL is written in Python and C code. To build it, it needs Python version 2.7 and a GCC compiler.

It uses libdivsufsort for suffix array construction and the probcons code for refinement of the graph.

Furthermore it uses the Python packages networkx (version 2), intervaltree and matplotlib.

REVEAL can be installed by executing the following command:

**python setup.py test install**

To install without executing the unit tests:

**python setup.py install**

## RUN

To validate whether everything is correctly installed you can run a test alignment from the tests directory, e.g. by executing the following command:

**reveal align tests/1a.fa tests/1b.fa**

This will output a shell script that outlines the typical steps to generate some graphs. If you're not interested in changing any parameters or the intermediate steps, you can immediately execute the script by piping it into your shell:

**reveal align tests/1a.fa tests/1b.fa | sh**

If everything ran correctly, various gfa files should have been produced. Most likely you will be interested in 'prg.unzipped.realigned.gfa'. This file contains a reference graph in GFA format (see [GFA](http://lh3.github.io/2014/07/19/a-proposal-of-the-grapical-fragment-assembly-format/)).

By default reveal will try to anchor the alignment by simultaneously aligning all genomes, however, if this is unwanted (due to e.g. memory constraints), you can run the following command to generate a shell script that anchors the alignment in a hierarchically way in batches of for instance 5 genomes at a time:

**reveal align tests/1a.fa tests/1b.fa --order=sequential --chunksize=5 | sh**

All commands in the shell script that are in between the comment lines can be run in parallel in case you're running on a compute cluster.

There are other subcommands for reveal, for which some are used by the generated shell script:

To generate an anchor graph using the recursive exact matching approach for more than two sequences you can either call:

**reveal rem tests/1a.fa tests/1b.fa tests/1c.fa**

or by aligning a sequence to an existing gfa graph:

**reveal rem 1a_1b.gfa tests/1c.fa**

or align two graphs:

**reveal rem 1a_1b.gfa 1c_1d.gfa**

Important parameters for **reveal rem**  are -m and -n. See subcommand help.

With REVEAL a global alignment between chromosome length assemblies is assumed. To address the issues that follow from draft assemblies, a 'finish' subcommand is supplied that orders and orients contigs/scaffolds with respect to a reference genome and produces pseudo molecules for the draft assembly.

**reveal finish reference.fasta draft.fasta**

To address large events (like translocations, inversions, but also misassemblies) that prevent a colinear alignment between two genomes, the following command can be used to transform a structurally rearranged (draft) genome such that it conforms to the layout of the reference sequence.

**reveal finish --order=chains reference.fasta draft.fasta**

Have a look at the various parameters, especially: --mineventsize, --minchainsum and -m.

To obtain a graph-based representation that encodes the original as well as the 'transformed' genome as separate paths through a graph, use:

**reveal transform reference.fasta draft.fasta**

Note that the resulting graphs may contain cycles, but can still be used in subsequent alignments using REVEAL rem, as only the transformed paths that correspond to the reference layout will be used for segmenting the graph. Paths in the graph prefixed with an asterisk (\*) correspond to the original (non-transformed) input genomes, which are ignored by REVEAL during graph traversal and mainly function as a way to record the structural events that are present in the graph.

To extract bubbles (a list of source/sink pairs and nodes within the bubble) from a graph run:

**reveal bubbles 1a&#95;1b.gfa**

Similar to bubbles, but will print the actual varying sequence.

**reveal variants 1a&#95;1b.gfa**

To output statistics with respect to the number of nodes, bubbles, variants, aligned sequence etc.:

**reveal stats 1a&#95;1b.gfa**

To realign parts of the graph using a basepair resolution multiple sequence alignment method (instead of MUMs):

**reveal refine \<graph\> \<source-node\> \<sink-node\>**

To realign all bubbles: 

**reveal refine \<graph\> --all**

Note that when bubbles are larger than let's say 10000bp, this won't work, so have a look at different filtering options (e.g. --maxsize).

As the boundaries of Maximal Unique Matches are somewhat greedy, more accurate variant calls are obtained by first 'unzipping' bubbles before applying **reveal refine**. To unzip all bubbles 10bp, run:

**reveal unzip \<graph\> -u10**

To construct an interactive (-i, for zooming purposes) mumplot of two fasta files:

**reveal plot genome1.fasta genome2.fasta -i**

Or to visualise a graph in a mumplot

**reveal gplot 1a&#95;1b.gfa -i**

NOTE that you need to have matplotlib installed for these commands.

In case you want to inspect the graph with software like cytoscape or gephi, you can produce a graph in gml format by calling reveal as follows:

**reveal convert prg.gfa**

To extract a genome/path from the graph:

**reveal extract \<graph\> \<pathname\>**

To reverse complement a graph:

**reveal comp \<graph\>**

To extract a subgraph of the graph, to for instance inspect a complex bubble structure:

**reveal subgraph \<graph\> \<node1\> ...  \<nodeN\>**

To merge multiple gfa graphs into a single gfa graph, while maintaining node-id space:

**reveal merge \<graph1\> \<graph2\> ...  \<graphN\>**

Or to do the opposite, split a graph by its connected components:

**reveal split \<graph\>**

For the rest, most commands should print a help function, when you specify **reveal \<subcommand\> -h**

How it works:

![Graph alignment using Recursive Exact Matching](https://github.com/jasperlinthorst/reveal/blob/master/reveal.gif)

