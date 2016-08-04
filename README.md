
#REVEAL
-----------------

REVEAL (REcursiVe Exact-matching ALigner) can be used to (multi) align whole genomes.

Preprint can be found here: [Scalable multi whole-genome alignment using recursive exact matching](http://www.biorxiv.org/content/early/2015/07/17/022715)

How it works:

![Graph alignment using Recursive Exact Matching](https://github.com/jasperlinthorst/reveal/reveal.gif)

##INSTALL

REVEAL is written in Python and C code. To build it, it needs Python version 2.7 and a GCC compiler.

It also depends on the libdivsufsort 2.0.1 package.


1 - Build and install libdivsufsort 2.0.1  
To build and install libdivsufsort you need to have CMAKE. After that's installed do the following:  

- git clone https://github.com/y-256/libdivsufsort.git (or wget https://github.com/y-256/libdivsufsort/archive/master.zip)
- cd libdivsufsort
- mkdir build
- cd build
- cmake -DBUILD\_DIVSUFSORT64:BOOL=ON -DCMAKE\_BUILD\_TYPE="Release" -DCMAKE\_INSTALL\_PREFIX="/usr/local" ..
- make install

Libdivsufsort should now be installed into your default installation directory (most likely /usr/local/lib).  

2 - Build and install REVEAL  
If libdivsufsort is installed, REVEAL can be installed by executing the following command:

**python setup.py install**

##RUN

To validate whether everything is correctly installed you can run a test alignment from the directory in which this readme is placed, by executing the following command:  

**reveal align tests/1a.fa tests/1b.fa**  

If everything is correctly installed, a file called 1a\_1b.gfa should have been produced. This file contains a reference graph in GFA format (see [GFA](http://lh3.github.io/2014/07/19/a-proposal-of-the-grapical-fragment-assembly-format/)).

Important parameters to consider when running a (multi) alignment are -m, -c and -n. See subcommand help.

In case you want to inspect the graph with software like cytoscape or gephi, you can produce a graph in gml format by calling reveal as follows:

**reveal align tests/1a.fa tests/1b.fa --gml**

or run:

**reveal convert 1a_1b.gfa**

To generate a graph for more than two sequences you can either call:

**reveal align tests/1a.fa tests/1b.fa tests/1c.fa**

or progressively align a sequence to an existing gfa graph:

**reveal align 1a_1b.gfa tests/1c.fa**

or align two graphs:

**reveal align 1a_1b.gfa 1c_1d.gfa**

To extract variants from a graph run:

**reveal bubbles 1a&#95;1b.gfa**

This will print a list of bubbles source/sink pairs and nodes within the bubble, the different variants (SNPs, indels, highly variable alleles etc.) and their positions (with respect to the genome specified as a reference) and the actual varying alleles.

To realign parts of a graph (e.g. with different settings or after progressive alignment):

**reveal realign \<graph\> \<source-node\> \<sink-node\>**

To construct an interactive (-i, for zooming purposes) mumplot of two fasta files (single contig for now...):

**reveal plot contig1.fasta contig2.fasta -i**

NOTE, you need to have matplotlib for this.

To extract a genome (name of the original fasta file) from the graph:

**reveal extract \<graph\> \<genome\>**

To reverse complement a graph:

**reveal comp \<graph\>**

To extract a subgraph of the graph, to for instance inspect a complex bubble structure:

**reveal subgraph \<graph\> \<node1\> ...  \<nodeN\>**

For the rest, most commands should print a help function, when you specify **reveal \<subcommand\> -h**

Under notebooks you can find some IPython notebooks that show some experiments for some typical usecases.
