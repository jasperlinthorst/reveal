
REVEAL v0.1 alpha
REVEAL can be used in order to (multi) align whole genomes.

INSTALL

REVEAL is written in Python and C code, so in order to build it it needs Python version 2.7 and a GCC compiler.

It also depends on the following 3d party packages:
	- libdivsufsort 2.0.1
	- seqal --> https://github.com/mhulsman/seqal



1) Build and install libdivsufsort 2.0.1
In order to build and install libdivsufsort you need to have CMAKE (versions should be available for most operating systems). After that's installed do the following:

- wget https://libdivsufsort.googlecode.com/files/libdivsufsort-2.0.1.tar.gz
- tar xvf libdivsufsort-2.0.1.tar.gz
- cd libdivsufsort-2.0.1
- mkdir build
- cd build
- cmake -DBUILD_DIVSUFSORT64:BOOL=ON -DCMAKE_BUILD_TYPE="Release" -DCMAKE_INSTALL_PREFIX="/usr/local" ..
- make install

Libdivsufsort should now be installed into your default installation directory (most likely /usr/local/lib).

2) Build and install seqal
Seqal is a python package that is used by REVEAL to perform Needleman-Wunsch alignments of larger bubbles to detect inversions.
To build and install it, do the following:

- git clone https://github.com/mhulsman/seqal.git
- cd seqal
- python setup.py install

3) Build and install REVEAL 
If libdivsufsort and seqal are installed, REVEAL can be installed by calling the following from the directory where this readme resides:
- python setup.py install

To validate whether everything is working the following can be run (again from the current directory) for a test alignment:
- python galn.py test/1a.fa tests/1b/fa

If everything is ok, the following output (or output alike):

03/03/2015 01:15:01 PM Identity (tot_aligned_sequence/(len(s1)+len(s2))) between 1a.fa and 1b.fa is 75.70%.
03/03/2015 01:15:01 PM 0 out of 318681 bases are not aligned.
03/03/2015 01:15:01 PM 0 vertices are not aligned.
03/03/2015 01:15:01 PM All graphs/sequences are aligned.

------------
