
#REVEAL v0.1 alpha
-----------------

REVEAL (REcursiVe Exact-matching ALigner) can be used to (multi) align whole genomes.

##INSTALL

REVEAL is written in Python and C code. To build it, it needs Python version 2.7 and a GCC compiler.

It also depends on the libdivsufsort 2.0.1 package.


1 - Build and install libdivsufsort 2.0.1
To build and install libdivsufsort you need to have CMAKE. After that's installed do the following:

- wget https://libdivsufsort.googlecode.com/files/libdivsufsort-2.0.1.tar.gz
- tar xvf libdivsufsort-2.0.1.tar.gz
- cd libdivsufsort-2.0.1
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

**reveal -k1000 tests/1a.fa tests/1b.fa**  
  
If everything is correctly installed, the produced output should end with the following lines:  
  
...  
03/03/2015 01:15:01 PM Identity (tot\_aligned\_sequence/(len(s1)+len(s2))) between 1a.fa and 1b.fa is 75.70%.  
03/03/2015 01:15:01 PM 0 out of 318681 bases are not aligned.  
03/03/2015 01:15:01 PM 0 vertices are not aligned.  
03/03/2015 01:15:01 PM All graphs/sequences are aligned.  

The following files should now have been produced:
- 1a\_1b.gfasta.gz - containing the alignment graph
- 1a\_1b.vcf.gz - containing variant calls in the vcf format (experimental)

