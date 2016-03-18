
#REVEAL
-----------------

REVEAL (REcursiVe Exact-matching ALigner) can be used to (multi) align whole genomes.

Preprint can be found here: [Scalable multi whole-genome alignment using recursive exact matching](http://www.biorxiv.org/content/early/2015/07/17/022715)

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

If everything is correctly installed, a file called 1a\_1b.gfa should have been produced. This file contains a reference graph in GFA format (see [GFA](http://lh3.github.io/2014/07/19/a-proposal-of-the-grapical-fragment-assembly-format/)). In case you want to inspect the graph with software like cytoscape or gephi, you can produce a graph in gml format by calling reveal as follows:

**reveal align tests/1a.fa tests/1b.fa --gml**

To generate a graph for more than two sequences you can either call:

**reveal align tests/1a.fa tests/1b.fa tests/1c.fa**

or progressively align a sequence to an existing gfa graph:

**reveal align 1a_1b.gfa tests/1c.fa**

or align two graphs:

**reveal align 1a_1b.gfa 1c_1d.gfa**

To extract variants from a GFA graph, you need to install [seqal](https://github.com/mhulsman/seqal), and run:

**reveal call 1a&#95;1b.gfa --all**

This will print a list of different variants (SNPs, indels, inversions) and their positions (with respect to the genome used as a reference if specified during alignment) and the actual varying alleles. Consider the --minlen and --maxlen arguments to filter the output variants by size. For now, inverions are detected by needle-man wunsch alignment of the reverse complemented alleles so be aware of the memory requirements of detecting large inversions from the graph. 
