from pipeco import *
import pipeco
import re
import os
import matplotlib.pyplot as plt
import matplotlib

class patho_pai :
	def __init__( self ) :
		global path_data, path_software
		global dir_result, dir_blast, dir_out
		global path_genome0, path_genome_faa, path_genome_gff
		path_data = os.path.join(os.getcwd(), "data")
		path_software = os.path.join(os.getcwd(), "software")
		path_output = os.path.join(os.getcwd(), "pipeco_out")
		path_genome0 = os.path.join(os.getcwd(), "input", "fasta")
		path_genome_faa = path_genome0.replace( "/fasta", "/faa" )
		path_genome_gff = path_genome0.replace( "/fasta", "/gff" )		
		dir_out = path_output + "/03.pai_out"
		os.system( "mkdir -p %s" %dir_out )
	@staticmethod

	def pai_align( ) :
		def pai_alignment_run( file_faa ) :
			os.system( "mkdir -p %s/alignment" %dir_out )
			file_udb = "%s/pai_ecoli/pai_ecoli.udb" %path_data
			file_usearch = "%s/usearch" %path_software
			blast_out0 = "%s/alignment/%s" %( dir_out, file_faa.split( "/" )[ -1 ].replace( ".faa", ".aln" ) )
			blast_out1 = "%s/alignment/%s" %( dir_out, file_faa.split( "/" )[ -1 ].replace( ".faa", ".b6" ) )
			blast_out2 = "%s/alignment/%s" %( dir_out, file_faa.split( "/" )[ -1 ].replace( ".faa", ".m8" ) )
			m8_userfields = "-userfields query+target+id+alnlen+mism+opens+qlo+qhi+tlo+thi+evalue+bits+ql+tl"
			command = "%s -ublast %s -db %s -evalue 1e-9 -alnout %s -blast6out %s -userout %s %s" %( file_usearch, file_faa, file_udb, blast_out0, blast_out1, blast_out2, m8_userfields )
			(exitstatus, outtext) = subprocess.getstatusoutput( "%s" % command )
			
		def align_parse1( grep_m8 ) :
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

		def align_parse2( fil_tsv ) :
			alignment_tsv = fil_tsv.replace( "_fil.tsv", "_pai_alignment.tsv" )
			genome = fil_tsv.split( "/" )[ -1 ].strip( "_fil.tsv" )
			gff_file = "%s/%s.gff" %( path_genome_gff, genome )
			try :
				df_fil = pd.read_csv( fil_tsv, sep = "\t" )
				gff_data = dict()
				with open( gff_file, "r") as file :
					for line in file :
						if not line.startswith( "#" ) :
							field = line.strip().split( "\t" )
							if len(field) >= 9 :
								attributes = dict( attr.split( "=" ) for attr in field[ 8 ].split( ";" ) )
								locus_tag = attributes.get( "locus_tag", "" )
								gff_data[ locus_tag ] = {
								"start": int( field[ 3 ] ),
								"end": int( field[ 4 ] ),
								"strand": field[ 6 ] }
				df_fil[ "locus_tag" ] = df_fil[ "query_locus_drop" ].str.split().str[ 0 ]
				df_fil[ "position" ] = df_fil[ "locus_tag" ].apply( lambda x: f"{gff_data[ x ][ 'start' ]}-{gff_data[ x ][ 'end' ]}" if gff_data[ x ][ 'strand' ] == "+" else f"c({gff_data[ x ][ 'start' ]}-{gff_data[ x ][ 'end' ]})" if x in gff_data else "" )
				df_fil[ "PAIDB_idx" ] = df_fil[ "target_hit" ].str.split( "_PAIDB" ).str[ 0 ]
				df_fil[ "PAI" ] = df_fil[ "target_hit" ].str.split( "|" ).str[ 1 ]
				df_fil[ "PAI_gene" ] = df_fil[ "target_hit" ].str.split( "|" ).str[ 5 ].str.split( "_\[" ).str[ 0 ]
				df_fil = df_fil.drop(columns=["query_locus_drop"])
				df_fil = df_fil[["locus_tag", "position", "PAIDB_idx", "PAI", "PAI_gene", "target_hit", "bitscore", "evalue", "pidentity", "coverage"]]
				df_fil = df_fil.sort_values(by="locus_tag").reset_index( drop = True )
				df_fil.to_csv( alignment_tsv, sep="\t", header = True )
			except Exception as e :
				print(f"An error occurred: {e}")
		
		def parse_gene_group_pai( grep_alignment, threshold ) :
			pai_result_tsv = grep_alignment.replace( "_pai_alignment.tsv", "_pai_group_parse.tsv" )
			try :
				df = pd.read_csv( grep_alignment, sep = "\t", index_col = 0 )
				n_rows = list()
				for _, row in df.iterrows() :
					if "," in row[ "locus_tag"] :
						spl_locus = row[ "locus_tag" ].split( "," )
						spl_position = row[ "position" ].split( ")," )
						for i, ( tag, pos ) in enumerate( zip( spl_locus, spl_position ) ) :
							n_row = row.copy()
							n_row[ "locus_tag"] = tag; n_row[ "position" ] = pos
							if "join" in pos : 
								if "c(" in pos :
									pos = pos.replace( "join(", "" ).strip( ")" )
									start = pos.split( "-" )[ 0 ].strip( "c(" )
									end = pos.split( "-" )[ -1 ]
									n_row[ "position" ] = f"c({start}-{end})"
								else : 
									pos = pos.strip( "join(" ).strip( ")" )
									start = pos.split( "-" )[ 0 ]; end = pos.split( "-" )[ -1 ]
									n_row[ "position" ] = f"{start}-{end}"
							else : 
								n_row[ "position" ] = pos.strip()
							n_rows.append( n_row )
					else : 
						n_rows.append( row ) 
				df_n = pd.DataFrame( n_rows )
				df_n[ [ "start", "end" ] ] = df_n[ "position" ].str.extract( r'.*?(\d+).*?(\d+)' )
				df_n[ "strand" ] = df[ "position" ].apply( lambda x: "-" if x.startswith( "c" ) else "+" )
				df_n[ "virulence" ] = df_n[ "target_hit" ].apply( lambda x: x.split( "|" )[ 4 ] if len( x.split( "|" ) ) > 4 else "N/A" )
				df_n[ "start" ] = df_n[ "start" ].astype( int )
				df_n[ "end" ] = df_n[ "end" ].astype( int )
				df_n[ "pai_idx" ] = df_n[ "PAIDB_idx" ]
				df_n[ "pai_gene" ] = df_n[ "PAI_gene" ]
				df_n["locus_tag_num"] = df_n[ "locus_tag" ].str.extract(r'PROKKA_(\d+)').astype(int)
				df_n = df_n.sort_values( "locus_tag_num" ).reset_index( drop = True )
				df_n[ "p_end" ] = df_n[ "end" ].shift( 1 )
				df_n[ "n_start" ] = df_n[ "start" ].shift( -1 )
				df_n[ "group_id" ] = 0
				current_group = 1
				df_n.at[ 0, "group_id" ] = current_group
				for i in range( 1, len( df_n ) ) :
					prev_end = df_n.at[ i-1, "end" ]
					curr_start = df_n.at[ i, "start" ]
					if curr_start - prev_end >= 0 and curr_start - prev_end <= threshold :
						df_n.at[i, "group_id"] = current_group
					else:
						current_group += 1
						df_n.at[ i, "group_id" ] = current_group
				df_n[ "group_id" ] = df_n[ "group_id" ].apply( lambda x: f"gene_group_{x}" )
				df_parse = df_n[ [ "locus_tag", "start", "end", "strand", "group_id", "pai_idx", "PAI", "pai_gene", "virulence", "bitscore", "evalue", "pidentity", "coverage" ] ]
				df_parse.to_csv( pai_result_tsv, sep = "\t", header = True )
			except Exception as e :
				print( f"An error occurred: {e}" )

		def pai_final( file_tsv ) :
			genome = file_tsv.split( "/" )[ -1 ].replace( "_pai_group_parse.tsv", "" )
			path_nt = path_genome0
			nt_files = glob.glob( f"{path_nt}/*.fna" )
			nt_file = None 
			for file in nt_files : 
				if genome in file : 
					nt_file = file
					break
			out_tsv = file_tsv.replace( "_pai_group_parse.tsv", "_pai_final_results.tsv" )
			def extract_seq( file_genome ) :
				dic_genome = dict()
				for record in SeqIO.parse( file_genome, "fasta" ) :
					genome_seq = str( record.seq )
					dic_genome[ record.id ] = str( record.seq )
				return dic_genome
			def cal_gc( sequence ) :
				gc_count = sequence.count( "G" ) + sequence.count( "C" )
				total_count = len( sequence )
				gc_ratio = ( gc_count / total_count ) * 100 
				return gc_ratio
			def host_gc( sequence ) :
				total_gc_content = 0; total_length = 0 
				for seq in sequence :
					gc_content = cal_gc( seq )
					seq_length = len( seq )
					total_gc_content += gc_content * seq_length
					total_length += seq_length 
				avg_gc = total_gc_content/ total_length if total_length > 0 else 0
				return avg_gc
			df = pd.read_csv( file_tsv, sep = "\t", index_col = 0 )
			extract_seq( nt_file )
			dic_genome_seq = extract_seq( nt_file )
			avg_host_gc = host_gc( dic_genome_seq.values() )
			results = list(); significant_groups = list()
			for group_id, group_df in df.groupby( "group_id" ) :
				c_len = 0; found = False; group_seq = ""
				group_start = group_df[ "start" ].min()
				group_end = group_df[ "end" ].max()
				#print ( group_start, group_end )
				for seq_id, seq in dic_genome_seq.items() :
					seq_len = len( seq ) 
					n_c_len = c_len + seq_len
					if not found and c_len <= group_start < n_c_len :
						adj_start = group_start - c_len
						if group_end <= n_c_len :
							adj_end = group_end - c_len
							group_seq += seq[ adj_start : adj_end ]
							found = True
							break
						else : 
							group_seq += seq[ adj_start : ]
							group_end_remain = group_end - n_c_len
							found = True 
					elif found :
						if group_end_remain < seq_len :
							group_seq += seq[ : group_end_remain ]
							break
						else :
							group_seq += seq
							group_end_remain -= seq_len
					c_len = n_c_len
				if group_seq : 
					group_gc = cal_gc( group_seq )
					difference = abs( group_gc - avg_host_gc )
					#significant = difference > 1.5 and any( group_df[ "virulence" ] == "virulence" )
					significant = difference > 1.5
					if significant :
						pai_counts = group_df[ "PAI" ].value_counts().to_dict()
						significant_groups.append( group_id )
					else :
						pai_counts = group_df[ "PAI" ].value_counts().to_dict()
					group_size = group_end - group_start + 1
					results.append( [ group_id, group_start, group_end, group_size, group_gc, avg_host_gc, difference, significant, pai_counts ] )
				else :
					results.append( [ group_id, None, None, None, None, None, None, None, "Could not found seq" ] )
			columns = [ "group_ID", "start", "end", "size", "group_GC", "host_GC", "diff_GC", "significant", "pai_counts" ]
			df_results = pd.DataFrame( results, columns = columns ) 
			df_results = df_results.sort_values( by = "group_ID" )
			df_results.to_csv( out_tsv, sep = "\t", header = True, index = False )
			significant_df = df[ df[ "group_id" ].isin( significant_groups ) ]
			significant_out_tsv = file_tsv.replace( "_pai_group_parse.tsv", "_pai_final_significant_results.tsv" )
			significant_df.to_csv( significant_out_tsv, sep = "\t", header = True, index = True )

		num_cores = multiprocessing.cpu_count()
		lst_genome = glob.glob( path_genome_faa + "/*.faa" )
		for genome0 in tqdm( lst_genome, desc = "Processing alignment" ) :
			pai_alignment_run( genome0 )
		lst_m8 = glob.glob( dir_out + "/alignment/*.m8" )
		for m8_file in tqdm( lst_m8, desc = "Processing pre-parse_1" ) :
			align_parse1( m8_file )
		lst_tsv = glob.glob( dir_out + "/alignment/*_fil.tsv" )
		for tsv_file in tqdm( lst_tsv, desc = "Processing pre-parse_2" ) :
			align_parse2( tsv_file )
		lst_pai_alignment = glob.glob( dir_out + "/alignment/*_pai_alignment.tsv" )
		for file_pai_align in tqdm( lst_pai_alignment, desc = "Processing gene grouping" ) :
			parse_gene_group_pai( file_pai_align, 2000 )
		lst_pai_group = glob.glob( dir_out + "/alignment/*_pai_group_parse.tsv" )
		for group_file in tqdm( lst_pai_group, desc = "Processing PAIs parse" ) :
			pai_final( group_file )
		lst_pai_files = glob.glob( dir_out + "/alignment/*_pai_final_results.tsv" )
		count_df = pd.DataFrame( columns = [ "genome", "detected_gene_group", "significant_group" ] )
		details_df = pd.DataFrame( columns = [ "idx", "genome", "sig_group_pai_count" ] )
		for file in lst_pai_files : 
			genome = os.path.basename( file ).replace( "_pai_final_results.tsv", "" )
			df = pd.read_csv( file, sep = "\t" )
			detected_gene_group = len( df[ "group_ID" ].unique() )
			significant_group = len( df[ ( df[ "significant" ] == True ) & ( df[ "size" ] >= 500 ) ] )
			count_df = count_df.append( {"genome": genome, "detected_gene_group": detected_gene_group, "significant_group": significant_group}, ignore_index = True )
			sig_groups = df[ ( df[ "significant" ] == True ) & ( df[ "size" ] >= 500 ) ]
			for idx, row in sig_groups.iterrows() :
				details_df = details_df.append( {"idx": idx, "genome": genome, "sig_group_pai_count": row["pai_counts"]}, ignore_index = True )
		count_df = count_df.sort_values( by = "genome" )
		count_df.to_csv( dir_out + "/count.tsv", sep="\t", index=False )
		details_df = details_df.sort_values( by = "genome" )
		details_df.to_csv( dir_out + "/pai_significant_details.tsv", sep = "\t", index = False )

		# Make a summary of PAI counts based on how many times they match to significant gene groups
		summary_sig_data = []
		genome_pai_counts = {}
		for _, row in details_df.iterrows():
			genome = row["genome"]
			pai_counts_str = row["sig_group_pai_count"]
			if genome not in genome_pai_counts:
				genome_pai_counts[genome] = {}
			if isinstance(pai_counts_str, str) and pai_counts_str != "Could not found seq":
				try:
					if pai_counts_str.startswith("{") and pai_counts_str.endswith("}"):
						import ast
						pai_dict = ast.literal_eval(pai_counts_str)
					else:
						continue
				except:
					continue
			else:
				pai_dict = pai_counts_str
			
			if not isinstance(pai_dict, dict):
				continue
				
			# Count PAI occurrences for this genome (count each PAI once per significant group)
			for pai in pai_dict.keys():
				if pai in genome_pai_counts[genome]:
					genome_pai_counts[genome][pai] += 1
				else:
					genome_pai_counts[genome][pai] = 1

		# Convert the nested dictionary to a list of dictionaries for DataFrame creation
		for genome, pai_counts in genome_pai_counts.items():
			for pai, count in pai_counts.items():
				summary_sig_data.append({
					"genome": genome,
					"PAI": pai,
					"count": count
				})

		# Create and save the summary DataFrame
		summary_sig_df = pd.DataFrame(summary_sig_data)
		summary_sig_df = summary_sig_df.sort_values(by=["genome", "count", "PAI"])
		summary_sig_df.to_csv(dir_out + "/summary_sig.tsv", sep="\t", index=False)

		df_tsv = pd.read_csv( dir_out + "/count.tsv", sep = "\t" )
		genomes = df_tsv[ "genome" ]
		detected_gene_groups = df_tsv[ "detected_gene_group" ]
		significant_groups = df_tsv[ "significant_group" ]
		non_significant_groups = detected_gene_groups - significant_groups
		x = range( len( genomes ) )
		x_labels = genomes.tolist() # Used full names for labelling
		fig, ax = plt.subplots(figsize=(8, 5))
		ax.bar( x, non_significant_groups, width = 0.6, label = "Non-significant Groups" )
		ax.bar( x, significant_groups, width = 0.6, bottom = non_significant_groups, label = "Significant Groups" )
		ax.set_xticks( x )
		ax.set_xticklabels( x_labels, rotation = 45, ha = "right")
		ax.set_ylabel( "# of Gene Groups" )
		ax.legend()
		for i, v in enumerate( detected_gene_groups ) :
			ax.text( i, v + 0.1, str( v ), color = "black", ha = "center" )
		plt.tight_layout()
		plt.savefig( dir_out + "/Sig_nonSig_group.svg" )
		plt.savefig( dir_out + "/Sig_nonSig_group.png" )

if __name__ == "__main__":
    pa = patho_pai()
    pa.pai_align()	