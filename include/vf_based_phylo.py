from pipeco import *
import pipeco
import re
from Bio import AlignIO
from Bio import SeqIO
import tempfile

class vf_phylo :
	def __init__( self ) :
		global path_data, path_software
		global dir_result, dir_blast, dir_out
		global path_genome_faa, path_genome_ffn
		path_data = os.path.join(os.getcwd(), "data") # Allows scripts to be run independently
		path_software = os.path.join(os.getcwd(), "software")
		path_output = os.path.join(os.getcwd(), "pipeco_out")
		path_genome0 = os.path.join(os.getcwd(), "input", "fasta")
		path_genome_faa = path_genome0.replace( "/fasta", "/faa" )
		path_genome_ffn = path_genome0.replace( "/fasta", "/ffn" )
		dir_out = path_output + "/02.vf_phylo_out"
		os.system( "mkdir -p %s" %dir_out )
		# dir_genome = pipeco.dir_genome
		# dir_blast = pipeco.dir_blast
		# dir_result = pipeco.dir_result
		# dir_out = dir_result + "/marker_out"	
		# os.system( "mkdir -p %s" %dir_out )
	@staticmethod

	def vf_phylo_run( ) :
		def read_fasta2( file_fasta ) :
			dic_fasta = dict()
			for record in SeqIO.parse( open( file_fasta, "r" ), "fasta" ) :
				seq_id = record.description
				seq = record.seq
				dic_fasta[ seq_id ] = seq
			return dic_fasta

		def vf_alignment_run( file_faa ) :
			os.system( "mkdir -p %s/alignment" %dir_out )
			file_udb = "%s/patho_vf/vf_hit.udb" %path_data
			file_usearch = "%s/usearch" %path_software
			blast_out0 = "%s/alignment/%s" %( dir_out, file_faa.split( "/" )[ -1 ].replace( ".faa", ".aln" ) )
			blast_out1 = "%s/alignment/%s" %( dir_out, file_faa.split( "/" )[ -1 ].replace( ".faa", ".b6" ) )
			blast_out2 = "%s/alignment/%s" %( dir_out, file_faa.split( "/" )[ -1 ].replace( ".faa", ".m8" ) )
			m8_userfields = "-userfields query+target+id+alnlen+mism+opens+qlo+qhi+tlo+thi+evalue+bits+ql+tl"
			command = "%s -ublast %s -db %s -evalue 1e-9 -alnout %s -blast6out %s -userout %s %s" %( file_usearch, file_faa, file_udb, blast_out0, blast_out1, blast_out2, m8_userfields )
			#print ( command )
			(exitstatus, outtext) = subprocess.getstatusoutput( "%s" % command )
			#print ( outtext )
		def align_parse1( grep_m8 ) :
			# path_temp = "/".join( grep_m8.split( "/" )[ :-1 ] )
			# w_cov_files = glob.glob( os.path.join( path_temp, "*_w_cov.m8" ) )
			# if not w_cov_files :
			grep_m8_file0 = grep_m8.split( "/" )[ -1 ].replace( ".m8", "_w_cov.m9" )
			grep_m8_file1 = grep_m8.split( "/" )[ -1 ].replace( ".m8", "_fil.tsv" )
			df_m8_edit0 = "%s/alignment/%s" %( dir_out, grep_m8_file0 )
			df_m8_edit1 = "%s/alignment/%s" %( dir_out, grep_m8_file1 )
			if os.path.exists( df_m8_edit0 ) and os.path.exists( df_m8_edit1 ) :
				print ( f"Files already exist: {df_m8_edit0} and {df_m8_edit1}" )
				return
			try :
				df_m8 = pd.read_csv( grep_m8, sep = "\t", header = None )
				#display( df_m8 )
				df_m8.columns = [ "query_locus", "target_hit", "pidentity", "length", "mismatch", "gapopen", "qstart", "qend", "tstart", "tend", "evalue", "bitscore", "qlen", "tlen" ]
				df_m8[ "coverage" ] = ( df_m8[ "length" ] / df_m8[ "tlen" ] ) * 100 
				df_m8 = df_m8.reset_index()
				df_m8[ "query_locus_drop" ] = df_m8.query_locus.str.split( " #" ).str[ 0 ]
				df_m8.loc[ df_m8.coverage > 100, "coverage" ] = 100.0
				df_m8_drop = df_m8[ [ "query_locus_drop", "target_hit", "bitscore", "evalue", "pidentity", "coverage" ] ]
				df_filter_group0 = df_m8_drop[ df_m8_drop.apply( lambda x: x.coverage >= 90, axis = 1 ) ]
				df_filter_group1 = df_m8_drop[ (df_m8_drop.coverage >= 90) & (df_m8_drop.pidentity >= 80) ]
				df_filter_group1 = df_filter_group1.reset_index( drop = True )
				df_fil_drop = df_filter_group1.drop_duplicates( [ "query_locus_drop" ], keep = "first" )
				df_filter_group0.to_csv( df_m8_edit0, index = False, sep = "\t" )
				df_fil_drop.to_csv( df_m8_edit1, index = False, sep = "\t", header = True )
			except :
				print ( "error" )
				shutil.copy( grep_m8, df_m8_edit0 )
				shutil.copy( grep_m8, df_m8_edit1 )

		def parse_hit_to_seq( file_hit, file_genome, out_file ) :
			read_genome = read_fasta2( file_genome )
			read_tsv = pd.read_csv( file_hit, sep = "\t", header = 0 )
			hit_locus = read_tsv[ "query_locus_drop" ]
			with open( out_file, "w" ) as f :
				for key in hit_locus :
					if key in read_genome : 
						f.write( ">%s\n%s\n" %( key, read_genome[ key ] ) )
					else :
						print(f"Key {key} not found.")

		def trans_fasta( lst_cds_fasta, out_path ):
			for genome in lst_cds_fasta :
				temp_genome_name = genome.split( "/")[ -1].replace( ".fasta", "" )
				dic_fasta = read_fasta2( genome ) 
				out_modi_fasta = "%s/%s_modi.fasta" %( out_path, temp_genome_name )
				out_modi_info = "%s/%s_info.txt" %( out_path, temp_genome_name )
				with open( out_modi_fasta, "w" ) as f0, open( out_modi_info, "w" ) as f1 :
					for idx, ( k, v ) in enumerate( dic_fasta.items() ) :
						f0.write( ">%s\n%s\n" %( idx, v ) )
						f1.write( ">%s|%s\n%s\n" %( idx, k, v ) )   

		def vf_hit_prokka( lst_fasta, path_out ) :
			for idx, genome in tqdm( enumerate( lst_fasta ), total = len( lst_fasta ), desc = "Hit-annot Processing" ) :
				genome_name = genome.split( "/" )[ -1 ].replace( ".fasta", "" )
				prokka_command = "prokka --outdir %s --prefix %s --force --cpus 32 %s" %( path_out, genome_name, genome )
				subprocess.run( prokka_command, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE )
				for file_name in os.listdir( path_out ) :
					if not file_name.endswith( ".gff" ) :
						file_path = os.path.join( path_out, file_name )
						os.remove( file_path )

		def panaroo_default( gff, dir_temp, num_cores ) :
			command = "panaroo -i %s -c 0.6 --len_dif_percent 0.9 -o %s -t %s --clean-mode moderate" %( gff, dir_temp, num_cores )
			#print ( command )
			log_file = os.path.join( dir_temp, "panaroo_log.txt" )
			with open( log_file, "w", buffering = 1 ) as log :  
				process = subprocess.Popen( command, shell = True, stdout = subprocess.PIPE, stderr = subprocess.STDOUT, universal_newlines = True )
				pbar = tqdm( total = 100, desc = "Panaroo Progress_default" )
				for line in iter(process.stdout.readline, '') :
					try :
						log.write( line )
						if "%" in line :
							try :
								progress = float( line.split( "%" )[ 0 ].split()[ -1 ] )
								pbar.update( progress - pbar.n )
							except ValueError :
								pass
					except IOError as e :
						print(f"Error writing to log file: {e}")
				pbar.close()
				process.stdout.close()
				return_code = process.wait()

		def panaroo_merged( draft_merged_graph, backbone, dir_out, num_cores ) :
			command = "panaroo-merge -d %s %s -t %s -o %s" %( draft_merged_graph, backbone, num_cores, dir_out )
			#print ( command )
			log_file = os.path.join( dir_out, "panaroo_log.txt" )
			with open( log_file, "w", buffering = 1 ) as log :  
				process = subprocess.Popen( command, shell = True, stdout = subprocess.PIPE, stderr = subprocess.STDOUT, universal_newlines = True )
				pbar = tqdm( total = 100, desc = "Panaroo Progress_merged" )
				for line in iter(process.stdout.readline, '') :
					try :
						log.write( line )
						if "%" in line :
							try :
								progress = float( line.split( "%" )[ 0 ].split()[ -1 ] )
								pbar.update( progress - pbar.n )
							except ValueError :
								pass
					except IOError as e :
						print(f"Error writing to log file: {e}")
				pbar.close()
				process.stdout.close()
				return_code = process.wait()

		def transform_matrix( df ) :
			df_transposed = df.drop( "Gene", axis = 1 ).transpose()
			df_transposed = df_transposed.replace( 1, "A" )
			df_transposed = df_transposed.replace( 0, "T" )
			return df_transposed

		def create_matrix( df, file_out ) :
			with open( out_file, "w" ) as f :
				for strain, sequence in df.iterrows() :
					out_seq = "".join( sequence.astype( str ) )
					f.write( f">{strain}\n{out_seq}\n")
			print ( "Wrote... %s" %file_out )
		
		def unwrap_fasta(input_file): # TODO removes newlines due to issues with wrapped files - likely could be fixed more elegantly. 
			input_dir = os.path.dirname(input_file)
			with tempfile.NamedTemporaryFile(mode="w", dir=input_dir, delete=False) as tmp:
				# Write unwrapped sequences to temp file
				for record in SeqIO.parse(input_file, "fasta"): # Change format to RAxML-compatible
					seq=str(record.seq)
					seq = seq.replace("-", "0")  # Replace gaps with absences
					seq = seq.replace("a", "0")  # Replace absences
					seq = seq.replace("t", "1")  # Replace presence
					seq = seq.replace("\n", "")  # Remove newlines
					record.seq = record.seq.__class__(seq)
					SeqIO.write(record, tmp, "fasta")
				tmp_name = tmp.name
    
			# Replace original file
			shutil.move(tmp_name, input_file)
		
		def convert_aln_to_phy(aln_path, phylip_path, mapping_file):
			unwrap_fasta(aln_path)
			alignment = AlignIO.read(aln_path, "fasta")
			seq_mapping = {}
			# Assign each sequences a <10 character ID for .phy conversion
			for i, record in enumerate(alignment):
				short_id = f"ID_{i:04d}"
				seq_mapping[short_id] = record.id
				record.id = short_id
				record.description = ""
			AlignIO.write(alignment, phylip_path, "phylip")
			# Save seq_mapping file
			with open(mapping_file, 'w') as mf:
				for short, full in seq_mapping.items():
					mf.write(f"{short}\t{full}\n")
			

		def msa_mafft( in_file, out_aln, num_cores ) :
			command = "mafft --auto --thread %s %s > %s" %( num_cores, in_file, out_aln )
			with tqdm(total=100, desc="Seq-align Progress") as pbar :
				pbar.update( 1 )
				process = subprocess.Popen(command, shell=True, stderr=subprocess.PIPE, universal_newlines=True)
				while process.poll() is None :
					time.sleep( 0.1 )
				pbar.update( 99 )
			return os.path.exists(out_aln) and os.path.getsize(out_aln) > 0
		
		def make_phylo_raxml(out_aln, out_nwk): # Generate phylogenetic tree using RAxML with binary data model
			phylip_file = out_aln.replace(".aln", ".phy")
			convert_aln_to_phy(out_aln, phylip_file, mapping_file)
			raxml_dir = os.path.dirname(out_nwk)
			base_name = os.path.basename(out_nwk).replace(".nwk", "")
			# RAxML command for gene presence/absence (binary) data
			raxml_cmd = f"raxmlHPC -m BINCAT -p 12345 -s {phylip_file} -n {base_name} -w {raxml_dir} --no-bfgs"
			raxml_log = f"{raxml_dir}/RAxML_log_{base_name}.txt"
			
			with tqdm(total=100, desc="Tree(nwk) construction Progress") as pbar:
				pbar.update(1)
				try:
					# Run RAxML
					with open(raxml_log, 'w') as log:
						process = subprocess.Popen(raxml_cmd, shell=True, stderr=subprocess.PIPE, 
                                      stdout=subprocess.PIPE, universal_newlines=True)
					
						stdout, stderr = process.communicate()
						log.write("STDOUT:\n")
						log.write(stdout)
						log.write("\nSTDERR:\n")
						log.write(stderr)
					pbar.update(89)
					
					# Find the best tree and copy to destination
					best_tree = os.path.join(raxml_dir, f"RAxML_bestTree.{base_name}")
					if os.path.exists(best_tree):
						shutil.copyfile(best_tree, out_nwk)
						pbar.update(10)
						return True
					else:
						print("RAxML failed to generate output tree")
						# Check if there are any error files
						info_file = os.path.join(raxml_dir, f"RAxML_info.{base_name}")
						if os.path.exists(info_file):
							with open(info_file, 'r') as f:
								print("RAxML info contents:")
								print(f.read())
						return False
                
				except Exception as e:
					print(f"Error running RAxML: {str(e)}")
					return False
		
		def rename_tree_labels(tree_file, mapping_file):
			# Read the mapping file
			seq_mapping = {}
			with open(mapping_file, 'r') as mf:
				for line in mf:
					short, full = line.strip().split("\t")
					seq_mapping[short] = full
			# Read the tree file
			with open(tree_file, 'r') as tf:
				tree_str = tf.read()
			# Replace each short id with its corresponding full name
			for short, full in seq_mapping.items():
				tree_str = tree_str.replace(short, full)
    
			# Write the updated tree file back
			with open(tree_file, 'w') as tf:
				tf.write(tree_str)
				print(f"Tree labels updated in {tree_file}")
			
		num_cores = multiprocessing.cpu_count()
		#print ( num_cores )
		lst_genome = glob.glob( path_genome_faa + "/*.faa" )
		for genome0 in tqdm( lst_genome, desc = "Processing VF-based alignment" ) :
			vf_alignment_run( genome0 )
		lst_m8 = glob.glob( dir_out + "/alignment/*.m8" )
		for m8_file in tqdm( lst_m8, desc = "Processing pre-parse_1" ) :
			align_parse1( m8_file )
		os.system( "mkdir -p %s/vf_align_hit/faa" %dir_out )
		os.system( "mkdir -p %s/vf_align_hit/gff" %dir_out )
		os.system( "mkdir -p %s/vf_align_hit/cds_fasta" %dir_out )
		os.system( "mkdir -p %s/vf_align_hit/cds_fasta/mod_cds_fasta" %dir_out )
		lst_tsv = glob.glob( dir_out + "/alignment/*_fil.tsv" )
		#print ( lst_tsv )
		for tsv_file in tqdm( lst_tsv, desc = "Processing pre-parse_2" ) :
			temp_name = tsv_file.split( "/" )[ -1 ].replace( "_fil.tsv", "" )
			out_faa = "%s/vf_align_hit/faa/%s_vf_hit.faa" %( dir_out, temp_name )
			out_cds_fasta = "%s/vf_align_hit/cds_fasta/%s_vf_hit.fasta" %( dir_out, temp_name )
			file_genome0 = "%s/%s.faa" %( path_genome_faa, temp_name )
			file_genome1 = "%s/%s.ffn" %( path_genome_ffn, temp_name )
			parse_hit_to_seq( tsv_file, file_genome0, out_faa )
			parse_hit_to_seq( tsv_file, file_genome1, out_cds_fasta )
			lst_cds_fasta = glob.glob( dir_out + "/vf_align_hit/cds_fasta/*.fasta" )
			trans_fasta( lst_cds_fasta, "%s/vf_align_hit/cds_fasta/mod_cds_fasta" %dir_out )
			lst_vf_fasta = glob.glob( dir_out + "/vf_align_hit/cds_fasta/mod_cds_fasta/*.fasta" )
			#print ( lst_vf_fasta )
		vf_hit_prokka( lst_vf_fasta, "%s/vf_align_hit/gff" %dir_out )
		backbone = "%s/phylo_backbone" %path_data
		gff = "%s/vf_align_hit/gff/*.gff" %dir_out
		#file_muscle = path_software + "/muscle-linux-x86.v5.3" # No longer utilised
		if len( os.listdir( path_genome_faa ) ) > 1 :
			dir_temp_out = "%s/tmp/inGenome_out" %dir_out 
			dir_merged_out = "%s/tmp/Merged_pan_out" %dir_out 
			os.system( "mkdir -p %s" %dir_temp_out )
			os.system( "mkdir -p %s" %dir_merged_out )
			
			Rtab_file0 = "%s/gene_presence_absence.Rtab" %dir_temp_out
			Rtab_file1 = "%s/gene_presence_absence.Rtab" %dir_merged_out
			if os.path.exists( Rtab_file0 ) == True and os.path.getsize( Rtab_file0 ) >= 10 * 1024 :
				if os.path.exists( Rtab_file1 ) == True and os.path.getsize( Rtab_file1 ) >= 10 * 1024 :
					pass
				else : 
					panaroo_merged( dir_temp_out, backbone, dir_merged_out, num_cores ) # Swapped args 1 and 2 - Based on the function call I think this is the intended order
			else :
				panaroo_default( gff, dir_temp_out, num_cores )
				panaroo_merged( dir_temp_out, backbone, dir_merged_out, num_cores ) # Swapped args 1 and 2 - Based on the function call I think this is the intended order
			out_file = "%s/pan_matrix_input.fasta" %dir_merged_out
			out_aln = out_file.replace( ".fasta", ".aln" )
			out_nwk = out_file.replace( ".fasta", ".nwk" )
			mapping_file = out_aln.replace(".aln", ".tsv")
			df = pd.read_csv( Rtab_file1, sep = "\t" )
			df_transposed = transform_matrix( df )
			create_matrix( df_transposed, out_file )
			if os.path.exists( out_aln ) == True and os.path.getsize( out_aln  ) >= 100 * 1024 : 
				pass
			else : 
				msa_mafft( out_file, out_aln, num_cores )
				if msa_mafft(out_file, out_aln, num_cores):
					make_phylo_raxml(out_aln, out_nwk)	
					rename_tree_labels(out_nwk, mapping_file)
				else:
					print("MAFFT alignment failed, skipping tree construction")
		print ( "PIP-eco: VF-based Phylogenetic process Done..." )

if __name__ == "__main__":
    pa = vf_phylo()
    pa.vf_phylo_run()