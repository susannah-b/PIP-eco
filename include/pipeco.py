#for local
import os, sys, glob, shutil
import os.path as path
import pandas as pd
import numpy as np
import subprocess
import multiprocessing
import toytree
import toyplot.svg
import toyplot.png
import seaborn as sns
from Bio import Phylo
import matplotlib.pyplot as plt
from IPython.display import clear_output
from prokka import prokka
import marker_alignment, phylo, pai

class pipeline : 
	def __init__( self ) :
		print( "PIPeco pipeline initiating..." )
		global path_data, path_program
		global dir_wgs, dir_genome, dir_blast, dir_pan, dir_pai, dir_result
		dir_user = path.abspath( path.join( os.getcwd() ) )
		path_main = "pipeline/"
		path_data = "data/"
		path_program = "program/" 
		dir_wgs = path_main + "genome/sequence/1.fna" 
		dir_genome = path_main + "genome/annot/2.anno" 
		dir_blast = path_main + "genome/alignment/3.align" 
		dir_pan = path_main + "genome/pan_phylo/4.pan"
		dir_pai = path_main + "genome/pai/5.pai_align"
		dir_result = path_main + "PIPeco_result"
		#print( dir_user )
		os.system( "mkdir -p %s" % path_main )
		os.system( "mkdir -p %s" % dir_wgs )
		os.system( "mkdir -p %s" % dir_genome )
		os.system( "mkdir -p %s" % dir_blast )
		os.system( "mkdir -p %s" % dir_pan )
		os.system( "mkdir -p %s" % dir_pai )
		os.system( "mkdir -p %s" % dir_result )

	class method() :
		def __init__( self ) :
			pass
		@staticmethod
		def prokka() :
			prokka.prokka_sc()
		@staticmethod
		def marker_alignment() :
			call_class1_1 = marker_alignment.local_alignment()
			call_class1_1.blastp_run()
		@staticmethod
		def draft_assignment() :
			call_class1_2 = marker_alignment.local_alignment()
			call_class1_2.pathotype_assign()
		@staticmethod
		def pan_phylo() :
			call_class2_1 = phylo.pan_phylo()
			call_class2_1.phylo_analysis()
		@staticmethod
		def pai_analysis() :
			call_class3_1 = pai.pai_analysis0()
			call_class3_1.pai_align()

