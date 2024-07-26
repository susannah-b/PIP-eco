#for local
import os, sys, glob, shutil, time
import os.path as path
import multiprocessing
from multiprocessing import Manager
from multiprocessing import pool
import pandas as pd
import numpy as np
import subprocess
import json
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Phylo
from functools import partial
from scipy.stats import sem
from scipy.interpolate import make_interp_spline
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns
from tqdm.notebook import tqdm
import toytree
import toyplot.svg
import toyplot.png
from genome_annot import gannot
import marker_align, vf_based_phylo, pai
from IPython.display import clear_output

class pathotype : 
	def __init__( self ) :
		print( "PIPeco pipeline initiating..." )
		global path_data, path_software
		global path_input, path_output
		global dir_wgs, dir_genome, dir_blast, dir_pan, dir_pai, dir_result
		dir_user = path.abspath( path.join( os.getcwd() ) )
		path_input = "input/fasta"; path_output = "pipeco_out/"
		path_data = "data/"
		path_software = "software/" 
		os.system( "mkdir -p %s" %path_input )
		os.system( "mkdir -p %s" %path_output )
		os.system( "mkdir -p %s" %path_software )

	class method() :
		def __init__( self ) :
			pass
		@staticmethod
		def genome_annotation() :
			gannot.prokka_sc()
		@staticmethod
		def marker_alignment() :
			call_class1 = marker_align.pathotype_alignment()
			call_class1.marker_pathotype()
		@staticmethod
		def vf_based_phylogenetic() :
			call_class2 = vf_based_phylo.vf_phylo()
			call_class2.vf_phylo_run()
		@staticmethod
		def pai_analysis() :
			call_class3 = pai.patho_pai()
			call_class3.pai_align()

		# dir_wgs = path_main + "genome/sequence/1.fna" 
		# dir_genome = path_main + "genome/annot/2.anno" 
		# dir_blast = path_main + "genome/alignment/3.align" 
		# dir_pan = path_main + "genome/pan_phylo/4.pan"
		# dir_pai = path_main + "genome/pai/5.pai_align"
		# dir_result = path_main + "PIPeco_result"
		# #print( dir_user )
		# os.system( "mkdir -p %s" % path_main )
		# os.system( "mkdir -p %s" % dir_wgs )
		# os.system( "mkdir -p %s" % dir_genome )
		# os.system( "mkdir -p %s" % dir_blast )
		# os.system( "mkdir -p %s" % dir_pan )
		# os.system( "mkdir -p %s" % dir_pai )
		# os.system( "mkdir -p %s" % dir_result )
