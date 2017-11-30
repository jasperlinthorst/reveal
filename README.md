
#REVEAL
-----------------

REVEAL (REcursiVe Exact-matching ALigner) can be used to (multi) align whole genomes.

Preprint can be found here: [Scalable multi whole-genome alignment using recursive exact matching](http://www.biorxiv.org/content/early/2015/07/17/022715)

How it works:

![Graph alignment using Recursive Exact Matching](https://github.com/jasperlinthorst/reveal/blob/master/reveal.gif)

##INSTALL

REVEAL is written in Python and C code. To build it, it needs Python version 2.7 and a GCC compiler.

It uses libdivsufsort for suffix array construction and uses the Python packages networkx, intervaltree and matplotlib.

1 - Build and install REVEAL  
REVEAL can be installed by executing the following command:

**python setup.py test install**

To install without executing the unit tests:

**python setup.py install**


##RUN

To validate whether everything is correctly installed you can run a test alignment from the tests directory, e.g. by executing the following command:

**reveal align tests/1a.fa tests/1b.fa**

If everything is correctly installed, a file called 1a\_1b.gfa should have been produced. This file contains a reference graph in GFA format (see [GFA](http://lh3.github.io/2014/07/19/a-proposal-of-the-grapical-fragment-assembly-format/)).

Important parameters to consider when running a (multi) alignment are -m and -n. See subcommand help.

In case you want to inspect the graph with software like cytoscape or gephi, you can produce a graph in gml format by calling reveal as follows:

**reveal align tests/1a.fa tests/1b.fa --gml**

or run:

**reveal convert 1a_1b.gfa**

To generate a graph for more than two sequences you can either call:

**reveal align tests/1a.fa tests/1b.fa tests/1c.fa**

or by aligning a sequence to an existing gfa graph:

**reveal align 1a_1b.gfa tests/1c.fa**

or align two graphs:

**reveal align 1a_1b.gfa 1c_1d.gfa**

With REVEAL a global alignment between chromosome length assemblies are assumed. To address the issues that follow from draft assemblies, a 'finish' subcommand is supplied that orders and orients contigs/scaffolds with respect to a reference genome and produces pseudo molecules for the draft assembly.

**reveal finish reference.fasta draft.fasta**

To address large events (like translocations, inversions, but also misassemblies) that prevent a colinear alignment between two genomes, the following command, can be used to transform a structurally rearranged (draft) genome such that it conforms to the layout of the reference sequence.

**reveal finish --order=chains reference.fasta draft.fasta**

Have a look at the various parameters, especially: --mineventsize, --minchainsum and -m.

To obtain a graph-based representation that encodes the original as well as the 'transformed' genome as separate paths through a graph, use the --outputgraph option. Note that these graphs may contain cycles, but can still be used in subsequent alignments using REVEAL. Paths in the graph prefixed with an asterisk (\*) correspond to the original (non-transformed) input genomes, which are ignored by REVEAL during graph traversal and mainly function as a way to record structural events in a multi-genome alignment. By default, only the contigs in which a structural rearrangements was detected are output as a \*-path in order to save space.

**reveal finish --order=chains reference.fasta draft.fasta --outputgraph**

To extract bubbles (a list of bubbles source/sink pairs and nodes within the bubble) from a graph run:

**reveal bubbles 1a&#95;1b.gfa**

Similar to bubbles, but will print the actual varying sequence.

**reveal variants 1a&#95;1b.gfa**

To realign parts of a graph (e.g. with different settings or after progressive alignment):

**reveal realign \<graph\> \<source-node\> \<sink-node\>**

To construct an interactive (-i, for zooming purposes) mumplot of two fasta files (single contig for now...):

**reveal plot contig1.fasta contig2.fasta -i**

Or to visualise a graph in a mumplot

**reveal gplot 1a&#95;1b.gfa -i**

NOTE that you need to have matplotlib installed for these commands.

To extract a genome (name of the original fasta file) from the graph:

**reveal extract \<graph\> \<genome\>**

To reverse complement a graph:

**reveal comp \<graph\>**

To extract a subgraph of the graph, to for instance inspect a complex bubble structure:

**reveal subgraph \<graph\> \<node1\> ...  \<nodeN\>**

To merge multiple gfa graphs into a single gfa graph, while maintaining node-id space:

**reveal merge \<graph1\> \<graph2\> ...  \<graphN\>**

For the rest, most commands should print a help function, when you specify **reveal \<subcommand\> -h**

Under notebooks you can find some IPython notebooks that show some experiments for typical usecases (not up-to-date(!)).
