
#REVEAL
-----------------

REVEAL (REcursiVe Exact-matching ALigner) can be used to (multi) align whole genomes.

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

**reveal tests/1a.fa tests/1b.fa**  
  
If everything is correctly installed, the produced output should end with the following lines:  
  
...  
04/09/2015 09:59:37 AM Identity (tot_aligned_sequence/(len(s1)+len(s2))) between 1a.fa and 1b.fa is 76.63%.  
04/09/2015 09:59:37 AM Identity genome/graph 1a.fa: ((tot_aligned_sequence/2)/len(s1)) is 67.71%.  
04/09/2015 09:59:37 AM Identity genome/graph 1b.fa: ((tot_aligned_sequence/2)/len(s2)) is 88.26%.  
04/09/2015 09:59:37 AM 62.12% of the sequence in the graph is common to all input genomes  
04/09/2015 09:59:37 AM All graphs/sequences are aligned.  

The following files should now have been created:
- 1a\_1b.gfa - containing the alignment graph in GFA format (see [GFA](http://lh3.github.io/2014/07/19/a-proposal-of-the-grapical-fragment-assembly-format/))
- 1a\_1b.vcf.gz - containing variant calls in vcf format (experimental)

