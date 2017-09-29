from unittest import TestCase
from nose.tools import with_setup
from StringIO import StringIO

from reveal import reveal,align,utils

import networkx as nx
import sys
import os

def teardown():
    os.remove("1a_1b.gfa")
    os.remove("1c_1d.gfa")
    os.remove("1a_1b_1c.gfa")
    os.remove("1c_1a_1b.gfa")
    os.remove("1a_1b_1c_1d.gfa")
    os.remove("123a_123b.gfa")

@with_setup(None, teardown)
class TestReveal(TestCase):
    #order is important here, nose uses alphabetical order of function names
    bakargv=None
    bakout=None
    
    def setup(self):
        bakargv=sys.argv
        bakout=sys.stdout

    def teardown(self):
        sys.argv=bakargv
        sys.stdout=bakout

    def test01_seqpair_align(self):
        G,idx=align.align([("1","ACTGATGTAGCTAGCTA"),("2","ACTAGCTAGCTAGTCAG")],minlength=1)
        self.assertTrue(isinstance(G, nx.DiGraph))
        self.assertTrue(G.number_of_nodes()>2)
        self.assertTrue(G.number_of_edges()>2)
    
    @with_setup(setup, teardown)
    def test02_fastapair_align_cmd(self):
        sys.argv=['reveal','align','tests/1a.fa','tests/1b.fa']
        reveal.main()
        self.assertTrue(os.path.exists("1a_1b.gfa"))
    
    @with_setup(setup, teardown)
    def test03_64_fastapair_align_cmd(self):
        sys.argv=['reveal','--64','align','tests/1c.fa','tests/1d.fa']
        reveal.main()
        self.assertTrue(os.path.exists("1c_1d.gfa"))
    
    @with_setup(setup, teardown)
    def test04_fastamulti_align_cmd(self):
        sys.argv=['reveal','align','tests/1a.fa','tests/1b.fa','tests/1c.fa']
        reveal.main()
        self.assertTrue(os.path.exists("1a_1b_1c.gfa"))
    
    @with_setup(setup, teardown)
    def test05_graph2graph_align_cmd(self):
        sys.argv=['reveal','align','1a_1b.gfa','1c_1d.gfa']
        reveal.main()
        self.assertTrue(os.path.exists("1a_1b_1c_1d.gfa"))
    
    @with_setup(setup, teardown)
    def test06_fasta2graph_align_cmd(self):
        sys.argv=['reveal','align','tests/1c.fa','1a_1b.gfa']
        reveal.main()
        self.assertTrue(os.path.exists("1c_1a_1b.gfa"))

    @with_setup(setup, teardown)
    def test07_multifastapair_align_cmd(self):
        sys.argv=['reveal','align','tests/123a.fa','tests/123b.fa','-m1000']
        reveal.main()
        self.assertTrue(os.path.exists("123a_123b.gfa"))
    
    @with_setup(setup, teardown)
    def test08_bubbles_cmd(self):
        sys.stdout=StringIO()
        sys.argv=['reveal','bubbles','1a_1b_1c.gfa']
        reveal.main()
        v=sys.stdout.getvalue()
        lines=v.split('\n')
        self.assertTrue(len(lines)>0)
        self.assertTrue(lines[0].startswith("#"))
    
    @with_setup(setup, teardown)
    def test09_variants_cmd(self):
        sys.stdout=StringIO()
        sys.argv=['reveal','variants','1a_1b_1c.gfa']
        reveal.main()
        v=sys.stdout.getvalue()
        lines=v.split('\n')
        self.assertTrue(len(lines)>0)
        self.assertTrue(lines[0].startswith("#"))
    
    @with_setup(setup, teardown)
    def test10_stats_cmd(self):
        sys.stdout=StringIO()
        sys.argv=['reveal','stats','1a_1b_1c.gfa']
        reveal.main()
        v=sys.stdout.getvalue()
        lines=v.split('\n')
        self.assertTrue(len(lines)>0)
        self.assertTrue(lines[0].find(":")!=-1)
    
    @with_setup(setup, teardown)
    def test11_comp_cmd(self):
        sys.argv=['reveal','comp','1a_1b_1c.gfa']
        reveal.main()
        self.assertTrue(os.path.exists('1a_1b_1c.rc.gfa'))
        os.remove("1a_1b_1c.rc.gfa")
    
    @with_setup(setup, teardown)
    def test12_split_cmd(self):
        sys.argv=['reveal','split','123a_123b.gfa']
        reveal.main()
        self.assertTrue(os.path.exists('0.gfa'))
        self.assertTrue(os.path.exists('1.gfa'))
        self.assertTrue(os.path.exists('2.gfa'))
        os.remove("0.gfa")
        os.remove("1.gfa")
        os.remove("2.gfa")

    @with_setup(setup, teardown)
    def test13_realign_cmd(self):
        sys.argv=['reveal','realign','1a_1b_1c.gfa','15','19','-n2'] #should cause a complex bubble
        reveal.main()
        self.assertTrue(os.path.exists('1a_1b_1c.realigned.gfa'))
    
    @with_setup(setup, teardown)
    def test14_complexbubble_cmd(self):
        sys.stdout=StringIO()
        sys.argv=['reveal','bubbles','1a_1b_1c.realigned.gfa']
        reveal.main()
        found=False
        v=sys.stdout.getvalue()
        self.assertTrue(v[0]=='#')
        for line in v.split('\n'):
            if line.split("\t")[3]=='complex':
                found=True
                break
        self.assertTrue(found)
        os.remove('1a_1b_1c.realigned.gfa')
    
    @with_setup(setup, teardown)
    def test15_extract_cmd(self):
        sys.stdout=StringIO()
        sys.argv=['reveal','extract','1a_1b.gfa','ACJE01000011_BB']
        reveal.main()
        for name,seq in utils.fasta_reader("tests/1a.fa"):
            pass
        extracted=sys.stdout.getvalue()
        extracted=extracted[extracted.find('\n')+1:].replace("\n","")
        self.assertTrue(seq==extracted)

    def test16_finish_cmd(self):
        sys.argv=['reveal','finish','tests/1a.fa','tests/1b.fa']
        reveal.main()
        self.assertTrue(os.path.exists("1a_1b.fasta"))
        self.assertTrue(os.path.exists("1a_1b.unplaced.fasta"))
        os.remove("1a_1b.fasta")
        os.remove("1a_1b.unplaced.fasta")
    
    def test17_pairchain_cmd(self):
        sys.argv=['reveal','chain','tests/1a.fa','tests/1b.fa','-o','1a_1b.chain']
        reveal.main()
        self.assertTrue(os.path.exists("1a_1b.chain.gfa"))
        os.remove("1a_1b.chain.gfa")
    
    def test18_multichain_cmd(self):
        sys.argv=['reveal','chain','tests/1a.fa','tests/1b.fa','tests/1c.fa','-o','1a_1b_1c.chain']
        reveal.main()
        self.assertTrue(os.path.exists("1a_1b_1c.chain.gfa"))
        os.remove("1a_1b_1c.chain.gfa")
    
    def test19_convert_cmd(self):
        sys.argv=['reveal','convert','1a_1b.gfa','123a_123b.gfa']
        reveal.main()
        self.assertTrue(os.path.exists("1a_1b.gml"))
        self.assertTrue(os.path.exists("123a_123b.gml"))
        os.remove("1a_1b.gml")
        os.remove("123a_123b.gml")
    
    def test20_merge_cmd(self):
        sys.argv=['reveal','merge','1a_1b.gfa','1c_1d.gfa','-o','merged']
        reveal.main()
        self.assertTrue(os.path.exists("merged.gfa"))
        os.remove("merged.gfa")
    
    def test21_plot_cmd(self):
        pass


