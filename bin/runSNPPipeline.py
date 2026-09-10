#!/usr/bin/env python3

import sys
import os
import pandas as pd
import datetime
from pybedtools import BedTool,helpers
import concurrent.futures
import time
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO
from itertools import combinations
import numpy as np
import uuid
import traceback
import shutil
import argparse
from csp2_common import (safe_int, safe_float, fetchHeaders, processBED, processSNPs, swapHeader, parseSNPDiffs, calculate_total_length, filterSNPs)


    


    


    
def screenSNPDiffs(snpdiffs_file,trim_name, min_cov, min_len, min_iden, ref_edge, query_edge, density_windows, max_snps,reference_id,log_directory):
    
    screen_start_time = time.time()
    
    if temp_dir != "":
        helpers.set_tempdir(temp_dir)

    # Set CSP2 variables to NA
    csp2_screen_snps = purged_length = purged_identity = purged_invalid = purged_indel = purged_lengthIdentity = purged_duplicate = purged_het = purged_density = filtered_ref_edge = filtered_query_edge = filtered_both_edge = "NA"
    filtered_snp_df = pd.DataFrame()
    good_reference_bed_df = pd.DataFrame()
    
    # Ensure snpdiffs file exists
    if not os.path.exists(snpdiffs_file) or not snpdiffs_file.endswith('.snpdiffs'):
        run_failed = True       
        sys.exit(f"Invalid snpdiffs file provided: {snpdiffs_file}")
        
    # Ensure header can be read in
    try:
        header_data = fetchHeaders(snpdiffs_file)
        header_query = header_data['Query_ID'][0].replace(trim_name,'')
        header_ref = header_data['Reference_ID'][0].replace(trim_name,'')
    except:
        run_failed = True       
        sys.exit(f"Error reading headers from snpdiffs file: {snpdiffs_file}")
        
    # Check snpdiffs orientation
    if (header_ref == reference_id):
        snpdiffs_orientation = 1
        query_id = header_query
    elif (header_query == reference_id):
        snpdiffs_orientation = -1
        query_id = header_ref
        header_data = swapHeader(header_data)
    else:
        run_failed = True       
        sys.exit(f"Error: Reference ID not found in header of {snpdiffs_file}...")      

    # Establish log file
    log_file = f"{log_directory}/{query_id}__vs__{reference_id}.log"
    with open(log_file,"w+") as log:
        log.write("Reference Screening for SNP Pipeline Analysis\n")
        log.write(f"Query Isolate: {query_id}\n")
        log.write(f"Reference Isolate: {reference_id}\n")
        log.write(str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))+"\n")
        log.write("-------------------------------------------------------\n\n")
        if snpdiffs_orientation == 1:
            log.write("\t- SNPDiffs file is in the forward orientation\n")
            log.write("-------------------------------------------------------\n\n")
        else:
            log.write("\t- SNPDiffs file is in the reverse orientation\n")
            log.write("-------------------------------------------------------\n\n")

    
    # Set variables from header data
    raw_snps = safe_int(header_data['SNPs'][0])
    raw_indels = safe_int(header_data['Indels'][0])
    raw_invalid = safe_int(header_data['Invalid'][0])

    kmer_similarity = safe_float(header_data['Kmer_Similarity'][0])
    shared_kmers = safe_int(header_data['Shared_Kmers'][0])
    query_unique_kmers = safe_int(header_data['Query_Unique_Kmers'][0])
    reference_unique_kmers = safe_int(header_data['Reference_Unique_Kmers'][0])
    mummer_gsnps = safe_int(header_data['gSNPs'][0])
    mummer_gindels = safe_int(header_data['gIndels'][0])

    query_bases = safe_int(header_data['Query_Assembly_Bases'][0])
    reference_bases = safe_int(header_data['Reference_Assembly_Bases'][0])

    query_contigs = safe_int(header_data['Query_Contig_Count'][0])
    reference_contigs = safe_int(header_data['Reference_Contig_Count'][0])

    raw_query_percent_aligned = safe_float(header_data['Query_Percent_Aligned'][0])
    raw_ref_percent_aligned = safe_float(header_data['Reference_Percent_Aligned'][0])
        
    # If the reference is not covered by at least min_cov, STOP
    if raw_ref_percent_aligned < min_cov:
        query_percent_aligned = raw_query_percent_aligned
        reference_percent_aligned = raw_ref_percent_aligned
        screen_category = "Low_Coverage"
        with open(log_file,"a+") as log:
            log.write(f"\t- Reference genome coverage: {raw_ref_percent_aligned}% \n")
            log.write(f"\t- Query covers less than --min_cov ({min_cov}%)...Screen halted...\n")
            log.write("-------------------------------------------------------\n\n")

    elif raw_snps + raw_indels + raw_invalid > 10000:
        query_percent_aligned = raw_query_percent_aligned
        reference_percent_aligned = raw_ref_percent_aligned
        screen_category = "SNP_Cutoff"
        with open(log_file,"a+") as log:
            log.write(f"\t- {raw_snps} detected...\n")
            log.write("\t- > 10,000 SNPs, indels, or invalid sites detected by MUMmer...Screen halted...\n")
            log.write("-------------------------------------------------------\n\n")
   
    else:
    
        ##### 02: Read in BED/SNP data #####
        with open(log_file,"a+") as log:
                log.write("Step 1: Reading in snpdiffs BED/SNP data...")
        try:
            bed_df,snp_df = parseSNPDiffs(snpdiffs_file,snpdiffs_orientation)
            
            with open(log_file,"a+") as log:
                log.write("Done!\n")
                log.write("-------------------------------------------------------\n\n")
                
        except Exception as e:
            run_failed = True       
            with open(log_file,"a+") as log:
                log.write(f"\nError reading BED/SNP data from file: {snpdiffs_file}\n{str(e)}")
            sys.exit(f"Error reading BED/SNP data from file: {snpdiffs_file}\n{str(e)}")
            
        ##### 03: Filter genome overlaps #####
        with open(log_file,"a+") as log:
            log.write("Step 2: Filtering for short overlaps and low percent identity...")
        
        good_bed_df = bed_df[(bed_df['Ref_Aligned'] >= min_len) & (bed_df['Perc_Iden'] >= min_iden)].copy()
        
        if good_bed_df.shape[0] == 0:
            query_percent_aligned = raw_query_percent_aligned
            reference_percent_aligned = raw_ref_percent_aligned
            screen_category = "Low_Quality_Coverage"
            with open(log_file,"a+") as log:
                log.write(f"\n\t- After filtering based on --min_len ({min_len}) and --min_iden ({min_iden}) , no valid alignments remain...Screen halted...\n")
                log.write("-------------------------------------------------------\n\n")

                
        else:
            # Create a BED file for alignments that pass basic QC
            good_query_bed_df = good_bed_df[['Query_Contig','Query_Start','Query_End']].copy()
            good_reference_bed_df = good_bed_df[['Ref_Contig','Ref_Start','Ref_End']].copy()
            good_reference_bed_df.loc[:, 'Query_ID'] = query_id

            good_query_aligned = calculate_total_length(BedTool.from_dataframe(good_query_bed_df).sort().merge())
            good_reference_aligned = calculate_total_length(BedTool.from_dataframe(good_reference_bed_df[['Ref_Contig','Ref_Start','Ref_End']]).sort().merge())
            
            query_percent_aligned = (good_query_aligned / query_bases) * 100
            reference_percent_aligned = (good_reference_aligned / reference_bases) * 100
            
            if reference_percent_aligned < min_cov:
                screen_category = "Low_Quality_Coverage"
                with open(log_file,"a+") as log:
                    log.write(f"\n\t- Raw reference genome coverage was {raw_ref_percent_aligned}% \n")
                    log.write(f"\t- After filtering based on --min_len ({min_len}) and --min_iden ({min_iden}), reference genome coverage was {reference_percent_aligned:.2f}% \n")
                    log.write(f"\t- Query covers less than --min_cov ({min_cov}%) of reference after filtering...Screen halted...\n")
                    log.write("-------------------------------------------------------\n\n")

            else:
                screen_category = "Pass"
                with open(log_file,"a+") as log:
                    log.write("Done!\n")
                    log.write(f"\t- Raw reference genome coverage was {raw_ref_percent_aligned}% \n")
                    log.write(f"\t- After filtering based on --min_len ({min_len}) and --min_iden ({min_iden}), reference genome coverage was {reference_percent_aligned:.2f}% \n")
                    log.write("-------------------------------------------------------\n\n")

            
                # Filter SNPs
                with open(log_file,"a+") as log:
                    log.write("Step 3: Filtering SNPs to get final SNP distances...")
                
                if raw_snps == 0:
                    csp2_screen_snps = purged_length = purged_identity = purged_lengthIdentity = purged_indel = purged_invalid = purged_duplicate = purged_het = purged_density = filtered_ref_edge = filtered_query_edge = filtered_both_edge = 0
                    with open(log_file,"a+") as log:
                        log.write("Done!\n")
                        log.write("\t- No SNPs detected in MUMmer output, no filtering required\n")
                        log.write("-------------------------------------------------------\n\n")

                else:
                    filtered_snp_df = filterSNPs(snp_df,bed_df,log_file, min_len, min_iden, ref_edge, query_edge, density_windows, max_snps, temp_dir)
                    
                    csp2_screen_snps = filtered_snp_df[filtered_snp_df.Cat == "SNP"].shape[0]
                    purged_length = filtered_snp_df[filtered_snp_df.Cat == "Purged_Length"].shape[0]
                    purged_identity = filtered_snp_df[filtered_snp_df.Cat == "Purged_Identity"].shape[0]
                    purged_lengthIdentity = filtered_snp_df[filtered_snp_df.Cat == "Purged_LengthIdentity"].shape[0]
                    purged_invalid = filtered_snp_df[filtered_snp_df.Cat == "Purged_Invalid"].shape[0]
                    purged_indel = filtered_snp_df[filtered_snp_df.Cat == "Purged_Indel"].shape[0]
                    purged_het = filtered_snp_df[filtered_snp_df.Cat == "Purged_Heterozygous"].shape[0]                    
                    purged_duplicate = filtered_snp_df[filtered_snp_df.Cat == "Purged_Duplicate"].shape[0]
                    purged_density = filtered_snp_df[filtered_snp_df.Cat == "Purged_Density"].shape[0]
                    filtered_query_edge = filtered_snp_df[filtered_snp_df.Cat == "Filtered_Query_Edge"].shape[0]
                    filtered_ref_edge = filtered_snp_df[filtered_snp_df.Cat == "Filtered_Ref_Edge"].shape[0]
                    filtered_both_edge = filtered_snp_df[filtered_snp_df.Cat == "Filtered_Both_Edge"].shape[0]
                                             
                    # Write filtered SNP data to file
                    snp_file = log_file.replace(".log","_SNPs.tsv")
                    with open(snp_file,"w") as f:
                        filtered_snp_df.to_csv(f, sep="\t", index=False)
                    
                    filtered_snp_df.loc[:, 'Query_ID'] = query_id
                    
                    with open(log_file,"a+") as log:
                        log.write("Done!\n")
                        log.write(f"\t- {csp2_screen_snps} SNPs detected between {query_id} and {reference_id} after filtering\n")
                        log.write(f"\t- Full SNP data saved to {snp_file}\n")
                        log.write("-------------------------------------------------------\n\n")
    
    screen_end_time = time.time()
    with open(log_file,"a+") as log:
        log.write(f"Screening Time: {screen_end_time - screen_start_time:.2f} seconds\n")
    
    # Clean up pybedtools temp
    helpers.cleanup(verbose=False, remove_all=False)
    
    return ([str(item) for item in [query_id,reference_id,screen_category,csp2_screen_snps,
            f"{query_percent_aligned:.2f}",f"{reference_percent_aligned:.2f}",
            query_contigs,query_bases,reference_contigs,reference_bases,
            raw_snps,purged_length,purged_identity,purged_lengthIdentity,purged_invalid,purged_indel,purged_duplicate,purged_het,purged_density,
            filtered_query_edge,filtered_ref_edge,filtered_both_edge,
            kmer_similarity,shared_kmers,query_unique_kmers,reference_unique_kmers,
            mummer_gsnps,mummer_gindels]],good_reference_bed_df,filtered_snp_df)

def assessCoverage(query_id,site_list):
    
    if temp_dir != "":
        helpers.set_tempdir(temp_dir)
        
    if len(site_list) == 0:
        return pd.DataFrame(columns=['Ref_Loc','Query_ID','Cat'])
    else:
        coverage_df = pass_filter_coverage_df[pass_filter_coverage_df['Query_ID'] == query_id].copy()
        
        if coverage_df.shape[0] == 0:
            uncovered_loc_df = pd.DataFrame({
                'Ref_Loc': site_list,
                'Query_ID': [query_id] * len(site_list),
                'Cat': ["Uncovered"] * len(site_list)
            })
            return uncovered_loc_df
        else:
            coverage_bed = BedTool.from_dataframe(coverage_df[['Ref_Contig','Ref_Start','Ref_End']]).sort()
            snp_bed_df = pd.DataFrame([item.split('/') for item in site_list], columns=['Ref_Contig','Ref_End'])
            snp_bed_df['Ref_Start'] = snp_bed_df['Ref_End'].astype(float).astype(int) - 1
            snp_bed_df['Ref_Loc'] = site_list
            snp_bed = BedTool.from_dataframe(snp_bed_df[['Ref_Contig','Ref_Start','Ref_End','Ref_Loc']]).sort()
            
            # Ref_Locs from snp_bed that intersect with coverage_bed go into covered_locs, the rest go into uncovered_locs
            covered_locs = snp_bed.intersect(coverage_bed, wa=True)
            uncovered_locs = snp_bed.intersect(coverage_bed, v=True, wa=True)

            covered_loc_df = pd.DataFrame({
                'Ref_Loc': [snp.fields[3] for snp in covered_locs],
                'Query_ID': [query_id] * covered_locs.count(),
                'Cat': ["Ref_Base"] * covered_locs.count()
            }) if covered_locs.count() > 0 else pd.DataFrame(columns=['Ref_Loc','Query_ID','Cat'])

            uncovered_loc_df = pd.DataFrame({
                'Ref_Loc': [snp.fields[3] for snp in uncovered_locs],
                'Query_ID': [query_id] * uncovered_locs.count(),
                'Cat': ["Uncovered"] * uncovered_locs.count()
            }) if uncovered_locs.count() > 0 else pd.DataFrame(columns=['Ref_Loc','Query_ID','Cat'])
            
            # Clean up pybedtools temp
            helpers.cleanup(verbose=False, remove_all=False)
            
            return pd.concat([covered_loc_df.drop_duplicates(['Ref_Loc']),uncovered_loc_df])

def getPairwise(chunk, sequences, ids):
    results = []
    
    for i, j in chunk:
        seq1, seq2 = sequences[i], sequences[j]
        actg_mask1 = np.isin(seq1, list('ACTGactg'))
        actg_mask2 = np.isin(seq2, list('ACTGactg'))
        cocalled_mask = actg_mask1 & actg_mask2
        
        snps_cocalled = np.sum(cocalled_mask)
        snp_distance = np.sum((seq1 != seq2) & cocalled_mask)
        
        results.append([ids[i], ids[j], snp_distance, snps_cocalled])
    
    return results

def parallelAlignment(alignment, chunk_size=5000):
    sequences = [np.array(list(record.seq)) for record in alignment]
    ids = [record.id for record in alignment]
    pairwise_combinations = list(combinations(range(len(sequences)), 2))

    # Create chunks of pairwise combinations
    chunks = [pairwise_combinations[i:i + chunk_size] for i in range(0, len(pairwise_combinations), chunk_size)]

    results = []
    with concurrent.futures.ProcessPoolExecutor() as executor:
        future_to_chunk = {executor.submit(getPairwise, chunk, sequences, ids): chunk for chunk in chunks}
        for future in concurrent.futures.as_completed(future_to_chunk):
            chunk_results = future.result()
            results.extend(chunk_results)
    return results

def getFinalPurge(df):
    # Returns the 'farthest along' category for a given Ref_Loc
    if "Purged_Density" in df['Cat'].values:
        return "Purged_Density"
    elif "Purged_Heterozygous" in df['Cat'].values:
        return "Purged_Heterozygous"
    elif "Purged_Indel" in df['Cat'].values:
        return "Purged_Indel"
    elif "Purged_Invalid" in df['Cat'].values:
        return "Purged_Invalid"
    elif "Purged_Length" in df['Cat'].values:
        return "Purged_Length"
    else:
        return "Purged_Identity"

    
# Read in arguments
global run_failed
run_failed = False

start_time = time.time()

parser = argparse.ArgumentParser(description='CSP2 SNP Pipeline Analysis')
parser.add_argument('--reference_id', type=str, help='Reference Isolate')
parser.add_argument('--output_directory', type=str, help='Output Directory')
parser.add_argument('--log_directory', type=str, help='Log Directory')
parser.add_argument('--snpdiffs_file', type=str, help='Path to SNPdiffs file')
parser.add_argument('--min_cov', default=85,type=float, help='Minimum coverage')
parser.add_argument('--min_len', default=500,type=int, help='Minimum length')
parser.add_argument('--min_iden', default=99,type=float, help='Minimum identity')
parser.add_argument('--ref_edge', default=150,type=int, help='Reference edge')
parser.add_argument('--query_edge', default=150,type=int, help='Query edge')
parser.add_argument('--density_windows', default="1000,125,15",type=str, help='Density windows')
parser.add_argument('--max_snps', default="3,2,1", type=str, help='Maximum SNPs')
parser.add_argument('--trim_name', nargs='?', const="", default="", type=str, help='Trim name')
parser.add_argument('--max_missing',default=50, type=float, help='Maximum missing')
parser.add_argument('--tmp_dir',default="", type=str, help='Temporary directory')
parser.add_argument('--rescue', default="norescue",type=str, help='Rescue edge SNPs (rescue/norescue)')
args = parser.parse_args()

reference_id = args.reference_id
output_directory = os.path.abspath(args.output_directory)
log_directory = os.path.abspath(args.log_directory)
log_file = f"{output_directory}/CSP2_SNP_Pipeline.log"
snpdiffs_file = args.snpdiffs_file
min_cov = args.min_cov
min_len = args.min_len
min_iden = args.min_iden
ref_edge = args.ref_edge
query_edge = args.query_edge
density_windows = [int(x) for x in args.density_windows.split(",")]
max_snps = [int(x) for x in args.max_snps.split(",")]
trim_name = args.trim_name
max_missing = args.max_missing

# Establish log file
with open(log_file,"w+") as log:
    log.write("CSP2 SNP Pipeline Analysis\n")
    log.write(f"Reference Isolate: {reference_id}\n")
    log.write(str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))+"\n")
    log.write("-------------------------------------------------------\n\n")
    log.write("Reading in SNPDiffs files...")


# Read in all lines and ensure each file exists
snpdiffs_list = [line.strip() for line in open(snpdiffs_file, 'r')]
snpdiffs_list = [line for line in snpdiffs_list if line]
for snpdiffs_file in snpdiffs_list:
    if not os.path.exists(snpdiffs_file):
        run_failed = True
        sys.exit("Error: File does not exist: " + snpdiffs_file)

snpdiffs_list = list(set(snpdiffs_list))

if len(snpdiffs_list) == 0:
    run_failed = True
    sys.exit("No SNPdiffs files provided...")

with open(log_file, "a+") as log:
    log.write("Done!\n")
    log.write(f"\t- Read in {len(snpdiffs_list)} SNPdiffs files\n")
    log.write("-------------------------------------------------------\n\n")

global temp_dir
if args.tmp_dir != "":
    random_temp_id = str(uuid.uuid4())
    temp_dir = f"{os.path.normpath(os.path.abspath(args.tmp_dir))}/{random_temp_id}"
    try:
        os.mkdir(temp_dir)
        helpers.set_tempdir(temp_dir)
    except OSError as e:
        run_failed = True       
        print(f"Error: Failed to create directory '{temp_dir}': {e}")
else:
    temp_dir = ""

rescue_edge = str(args.rescue)
if rescue_edge not in ["rescue","norescue"]:
    with open(log_file,"a+") as log:
        log.write(f"\t- Unexpected rescue_edge variable ('{rescue_edge}'). Not performing SNP edge rescuing (any SNPs found within {query_edge}bp of a query contig edge will be purged)...\n")
        log.write("-------------------------------------------------------\n\n")
    rescue_edge = "norescue"
elif rescue_edge == "rescue":
    with open(log_file,"a+") as log:
        log.write(f"\t- Rescuing edge SNPs within {query_edge}bp of query contig edges if found more centrally in another query...\n")
        log.write("-------------------------------------------------------\n\n")
else:
    with open(log_file,"a+") as log:
        log.write(f"\t- Not performing SNP edge rescuing (any SNPs found within {query_edge}bp of a query contig edge will be purged)...\n")
        log.write("-------------------------------------------------------\n\n")

try:
    # Establish output files
    reference_screening_file = f"{output_directory}/Reference_Screening.tsv"
    locus_category_file = f"{output_directory}/Locus_Categories.tsv"
    query_coverage_file = f"{output_directory}/Query_Coverage.tsv"
    raw_loclist = f"{output_directory}/snplist.txt"
    raw_alignment = f"{output_directory}/snpma.fasta"
    preserved_loclist = f"{output_directory}/snplist_preserved.txt"
    preserved_alignment_file = f"{output_directory}/snpma_preserved.fasta"
    raw_pairwise = f"{output_directory}/snp_distance_pairwise.tsv"
    raw_matrix = f"{output_directory}/snp_distance_matrix.tsv"
    preserved_pairwise = f"{output_directory}/snp_distance_pairwise_preserved.tsv"
    preserved_matrix = f"{output_directory}/snp_distance_matrix_preserved.tsv"

    with open(log_file,"a+") as log:
        log.write("Screening all queries against reference...")
    with concurrent.futures.ProcessPoolExecutor() as executor:
        results = [executor.submit(screenSNPDiffs,snp_diff_file,trim_name, min_cov, min_len, min_iden, ref_edge, query_edge, density_windows, max_snps,reference_id,log_directory) for snp_diff_file in snpdiffs_list]

    # Combine results into a dataframe
    output_columns = ['Query_ID','Reference_ID','Screen_Category','CSP2_Screen_SNPs',
                'Query_Percent_Aligned','Reference_Percent_Aligned',
                'Query_Contigs','Query_Bases','Reference_Contigs','Reference_Bases',
                'Raw_SNPs','Purged_Length','Purged_Identity','Purged_LengthIdentity','Purged_Invalid','Purged_Indel','Purged_Duplicate','Purged_Het','Purged_Density',
                'Filtered_Query_Edge','Filtered_Ref_Edge','Filtered_Both_Edge',
                'Kmer_Similarity','Shared_Kmers','Query_Unique_Kmers','Reference_Unique_Kmers',
                'MUMmer_gSNPs','MUMmer_gIndels']

    # Save reference screening
    results_df = pd.DataFrame([item.result()[0] for item in results], columns = output_columns)
    with open(reference_screening_file,"w") as f:
        results_df.to_csv(f, sep="\t", index=False)

    # Get reference bed dfs
    covered_df = pd.concat([item.result()[1] for item in results])

    # Get snp dfs
    filtered_snp_df = pd.concat([item.result()[2] for item in results])

    # Separate isolates that pass QC
    pass_qc_isolates = list(set(results_df[results_df['Screen_Category'] == "Pass"]['Query_ID']))
    fail_qc_isolates = list(set(results_df[results_df['Screen_Category'] != "Pass"]['Query_ID']))

    if len(pass_qc_isolates) == 0:
        with open(log_file,"a+") as log:
            log.write("Done!\n")
            log.write(f"\t- Reference screening data saved to {reference_screening_file}\n")
            log.write(f"\t- Of {len(snpdiffs_list)} comparisons, no isolates passed QC. Pipeline cannot continue.\n")
            log.write(f"\t- {len(fail_qc_isolates)} comparisons failed QC\n")
            for isolate in fail_qc_isolates:
                isolate_category = results_df[results_df['Query_ID'] == isolate]['Screen_Category'].values[0]
                log.write(f"\t\t- {isolate}: {isolate_category}\n")
            log.write("-------------------------------------------------------\n\n")
        sys.exit("No sequences passed QC filters.")
    else:    
        with open(log_file,"a+") as log:
            log.write("Done!\n")
            log.write(f"\t- Reference screening data saved to {reference_screening_file}\n")
            log.write(f"\t- Of {len(snpdiffs_list)} comparisons, {len(pass_qc_isolates)} covered at least {min_cov}% of the reference genome after removing poor alignments\n")
            if len(fail_qc_isolates) > 0:
                log.write(f"\t- {len(fail_qc_isolates)} comparisons failed QC\n")
                for isolate in fail_qc_isolates:
                    isolate_category = results_df[results_df['Query_ID'] == isolate]['Screen_Category'].values[0]
                    log.write(f"\t\t- {isolate}: {isolate_category}\n")
            log.write("-------------------------------------------------------\n\n")

    with open(log_file,"a+") as log:
        log.write(f"Compiling SNPs across {len(pass_qc_isolates)} samples...\n")

    if filtered_snp_df.shape[0] == 0:
        snp_count = 0
        with open(log_file,"a+") as log:
            log.write("\t- No SNPs detected across samples...Skipping to output...\n")     
    else:

        # Remove samples that failed QC
        pass_filter_snps = filtered_snp_df[filtered_snp_df['Query_ID'].isin(pass_qc_isolates)].copy()
        pass_filter_snp_list = list(set(pass_filter_snps['Ref_Loc']))
        pass_filter_snp_count = len(pass_filter_snp_list)
    
        global pass_filter_coverage_df
        pass_filter_coverage_df = covered_df[covered_df['Query_ID'].isin(pass_qc_isolates)].copy()

        # Get SNP counts
        snp_df = pass_filter_snps[pass_filter_snps['Cat'] == "SNP"].copy()
        
        if snp_df.shape[0] == 0:
            snp_count = 0
            with open(log_file,"a+") as log:
                log.write(f"\t- {pass_filter_snp_count} total SNPs detected across all samples...\n")     
                log.write("\t- No SNPs passed QC filtering in any sample...Skipping to output...\n")
        
        else:
            snp_list = list(set(snp_df['Ref_Loc']))
            snp_count = len(snp_list)
            
            with open(log_file,"a+") as log:
                log.write(f"\t- {pass_filter_snp_count} total SNPs detected across all samples...\n")     
                log.write(f"\t- {snp_count} unique SNPs passed QC filtering in at least one sample...\n")
        
            # Note SNPs lost irrevocably to reference edge trimming
            ref_edge_df = pass_filter_snps[pass_filter_snps['Cat'].isin(["Filtered_Ref_Edge",'Filtered_Both_Edge'])].copy()    
            if ref_edge_df.shape[0] > 0:        
                with open(log_file,"a+") as log:
                    log.write(f"\t- {ref_edge_df['Ref_Loc'].nunique()} unique SNPs were within {ref_edge}bp of a reference contig end and were not considered in any query...\n")
        
            # Create Ref_Base df
            ref_base_df = snp_df[['Ref_Loc', 'Ref_Base']].copy().drop_duplicates().rename(columns = {'Ref_Base':'Query_Base'})
            ref_base_df.loc[:,'Query_ID'] = reference_id
            ref_base_df.loc[:,'Cat'] = "Reference_Isolate"
            ref_base_df = ref_base_df.loc[:,['Ref_Loc','Query_ID','Query_Base','Cat']]
            
            # Rescue SNPs that are near the edge if they are valid SNPs in other samples
            rescued_edge_df = pass_filter_snps[(pass_filter_snps['Cat'] == "Filtered_Query_Edge") & (pass_filter_snps['Ref_Loc'].isin(snp_list))].copy()
            
            if rescue_edge == "norescue":
                with open(log_file,"a+") as log:
                    log.write("\t- Skipping edge resucing...\n")
            
            elif rescued_edge_df.shape[0] > 0:
                
                # Remove rescued sites from pass_filter_snps
                rescue_merge = pass_filter_snps.merge(rescued_edge_df, indicator=True, how='outer')
                pass_filter_snps = rescue_merge[rescue_merge['_merge'] == 'left_only'].drop(columns=['_merge']).copy()            
                
                # Add rescued SNPs to snp_df
                rescued_edge_df.loc[:,'Cat'] = "Rescued_SNP"    
                snp_df = pd.concat([snp_df,rescued_edge_df]).reset_index(drop=True)
                with open(log_file,"a+") as log:
                    log.write(f"\t- Rescued {rescued_edge_df.shape[0]} query SNPs that fell within {query_edge}bp of the query contig edge...\n")
            
            else:
                with open(log_file,"a+") as log:
                    log.write(f"\t- No query SNPs that fell within {query_edge}bp of the query contig edge were rescued...\n")
                                
            # Gather base data for all valid SNPs
            snp_base_df = snp_df[['Ref_Loc','Query_ID','Query_Base','Cat']].copy()
            
            # Process purged sites
            purged_snp_df = pd.DataFrame(columns=['Ref_Loc','Query_ID','Query_Base','Cat'])
            purged_df = pass_filter_snps[pass_filter_snps['Cat'] != "SNP"].copy()
                    
            if purged_df.shape[0] > 0:

                # Remove rows from purged_df if the Ref_Loc/Query_ID pair is already in snp_base_df
                purged_df = purged_df[~purged_df[['Ref_Loc','Query_ID']].apply(tuple,1).isin(snp_base_df[['Ref_Loc','Query_ID']].apply(tuple,1))].copy()
                
                if purged_df.shape[0] > 0:

                    # Get purged SNPs where no query has a valid SNP
                    non_snp_df = purged_df[~purged_df['Ref_Loc'].isin(snp_list)].copy()
                    if non_snp_df.shape[0] > 0:
                        non_snp_merge = purged_df.merge(non_snp_df, indicator=True, how='outer')
                        purged_df = non_snp_merge[non_snp_merge['_merge'] == 'left_only'].drop(columns=['_merge']).copy()

                        with open(log_file,"a+") as log:
                            log.write(f"\t- {len(list(set(non_snp_df['Ref_Loc'])))} unique SNPs were purged in all queries they were found in, and were not considered in the final dataset...\n")
            
            if purged_df.shape[0] > 0:

                purged_snp_df = purged_df[['Ref_Loc','Query_ID']].copy().drop_duplicates()
                purged_snp_df.loc[:, 'Query_Base'] = "N"
                final_purge_df = purged_df.groupby(['Ref_Loc','Query_ID']).apply(getFinalPurge).reset_index().rename(columns={0:'Cat'})
                purged_snp_df = purged_snp_df.merge(final_purge_df, on=['Ref_Loc','Query_ID'], how='inner')

            # Genomic positions that do not occur- in the SNP data are either uncovered or match the reference base
            missing_df = pd.DataFrame(columns=['Ref_Loc','Query_ID','Query_Base','Cat'])

            covered_snps = pd.concat([snp_base_df,purged_snp_df]).copy()
            ref_loc_sets = covered_snps.groupby('Query_ID')['Ref_Loc'].agg(set).to_dict()
            isolates_with_missing = [isolate for isolate in pass_qc_isolates if len(set(snp_list) - ref_loc_sets.get(isolate, set())) > 0]
            
            uncovered_df = pd.DataFrame()
            
            if isolates_with_missing:
                isolate_data = [(isolate, list(set(snp_list) - ref_loc_sets.get(isolate, set()))) for isolate in isolates_with_missing]        
                
                with concurrent.futures.ProcessPoolExecutor() as executor:
                    results = [executor.submit(assessCoverage, query, sites) for query, sites in isolate_data]
                    coverage_dfs = [result.result() for result in concurrent.futures.as_completed(results)]

                coverage_df = pd.concat(coverage_dfs)
                covered_df = coverage_df[coverage_df['Cat'] == 'Ref_Base']
                uncovered_df = coverage_df[coverage_df['Cat'] == 'Uncovered']
                
                if not uncovered_df.empty:
                    uncovered_df.loc[:, 'Query_Base'] = "?"
                    missing_df = pd.concat([missing_df, uncovered_df[['Ref_Loc', 'Query_ID', 'Query_Base', 'Cat']]])
        
                if not covered_df.empty:
                    ref_base_snp_df = covered_df.merge(ref_base_df[['Ref_Loc', 'Query_Base']], on='Ref_Loc', how='left')
                    missing_df = pd.concat([missing_df, ref_base_snp_df[['Ref_Loc', 'Query_ID', 'Query_Base', 'Cat']]])        
            
            with open(log_file,"a+") as log:
                log.write("\t- Processed coverage information...\n")

            final_snp_df = pd.concat([snp_base_df,purged_snp_df,missing_df,ref_base_df]).sort_values(by=['Ref_Loc','Query_ID']).reset_index(drop=True)
            snp_counts = final_snp_df.groupby('Query_ID')['Ref_Loc'].count().reset_index().rename(columns={'Ref_Loc':'SNP_Count'})

            # Assert that all snp_counts == snp_count
            assert snp_counts['SNP_Count'].nunique() == 1
            assert snp_counts['SNP_Count'].values[0] == snp_count

            # Get locus coverage stats        
            snp_coverage_df = final_snp_df[final_snp_df['Cat'].isin(['SNP','Rescued_SNP'])].groupby('Ref_Loc')['Query_ID'].count().reset_index().rename(columns={'Query_ID':'SNP_Count'})
            ref_base_coverage_df = final_snp_df[final_snp_df['Cat'].isin(["Ref_Base","Reference_Isolate"])].groupby('Ref_Loc')['Query_ID'].count().reset_index().rename(columns={'Query_ID':'Ref_Base_Count'}) 
            
            if uncovered_df.shape[0] > 0:
                uncovered_count_df = final_snp_df[final_snp_df['Cat'] == "Uncovered"].groupby('Ref_Loc')['Query_ID'].count().reset_index().rename(columns={'Query_ID':'Uncovered_Count'}).copy()
            else:
                uncovered_count_df = pd.DataFrame(columns=['Ref_Loc','Uncovered_Count'])

            possible_purged_cols = ['Purged_Length','Purged_Identity','Purged_Invalid','Purged_Indel','Purged_Heterozygous','Purged_Density']
            if purged_snp_df.shape[0] > 0:
                purged_count_df = final_snp_df[final_snp_df['Cat'].isin(possible_purged_cols)].groupby('Ref_Loc')['Query_ID'].count().reset_index().rename(columns={'Query_ID':'Purged_Count'}).copy()
            else:
                purged_count_df = pd.DataFrame(columns=['Ref_Loc','Purged_Count'])

            locus_coverage_df = snp_coverage_df.merge(ref_base_coverage_df, how='outer', on='Ref_Loc').merge(uncovered_count_df, how='outer', on='Ref_Loc').merge(purged_count_df, how='outer', on='Ref_Loc').fillna(0)
            locus_coverage_df.loc[:, ['SNP_Count','Ref_Base_Count','Uncovered_Count','Purged_Count']] = locus_coverage_df.loc[:, ['SNP_Count','Ref_Base_Count','Uncovered_Count','Purged_Count']].astype(int)
            locus_coverage_df['Missing_Ratio'] = ((locus_coverage_df['Uncovered_Count'] + locus_coverage_df['Purged_Count']) / (1+len(pass_qc_isolates))) * 100
            with open(locus_category_file,"w") as f:
                locus_coverage_df.to_csv(f, sep="\t", index=False)
            
            # Get isolate coverage stats
            min_isolate_cols = ['Query_ID','SNP','Ref_Base','Percent_Missing','Purged','Uncovered','Rescued_SNP','Purged_Ref_Edge']
            isolate_coverage_df = final_snp_df.groupby('Query_ID')['Cat'].value_counts().unstack().fillna(0).astype(float).astype(int).reset_index().drop(columns=['Reference_Isolate'])
            isolate_coverage_df.loc[isolate_coverage_df['Query_ID'] == reference_id, 'Ref_Base'] = snp_count
            
            if "Rescued_SNP" not in isolate_coverage_df.columns.tolist():
                isolate_coverage_df.loc[:,'Rescued_SNP'] = 0
            isolate_coverage_df['SNP'] = isolate_coverage_df['SNP'] + isolate_coverage_df['Rescued_SNP']
            
            for col in ['Uncovered'] + possible_purged_cols:
                if col not in isolate_coverage_df.columns.tolist():
                    isolate_coverage_df.loc[:,col] = 0            
                
            isolate_coverage_df['Purged'] = isolate_coverage_df[possible_purged_cols].sum(axis=1)

            isolate_coverage_df['Percent_Missing'] = (isolate_coverage_df['Uncovered'] + isolate_coverage_df['Purged'])/(isolate_coverage_df['Uncovered'] + isolate_coverage_df['Purged'] + isolate_coverage_df['Ref_Base'] + isolate_coverage_df['SNP']) * 100
            
            isolate_coverage_df.loc[:,'Purged_Ref_Edge'] = 0
            if ref_edge_df.shape[0] > 0:
                isolate_coverage_df.loc[isolate_coverage_df['Query_ID'] == reference_id, 'Purged_Ref_Edge'] = ref_edge_df['Ref_Loc'].nunique()
            
            isolate_coverage_df = isolate_coverage_df[min_isolate_cols + possible_purged_cols].sort_values(by = 'Percent_Missing',ascending = False).reset_index(drop=True)
            with open(query_coverage_file,'w') as f:
                isolate_coverage_df.to_csv(f, sep="\t", index=False)
            
            with open(log_file,"a+") as log:
                log.write(f"\t- SNP coverage information: {locus_category_file}\n")
                log.write(f"\t- Query coverage information: {query_coverage_file}\n")
                log.write("-------------------------------------------------------\n\n")
            
            with open(log_file,"a+") as log:
                log.write("Processing alignment data...")

            alignment_df = final_snp_df[['Query_ID','Ref_Loc','Query_Base']].copy().rename(columns={'Query_Base':'Base'}).pivot(index='Query_ID', columns='Ref_Loc', values='Base')
            csp2_ordered = alignment_df.columns

            with open(raw_loclist,"w+") as loclist:
                loclist.write("\n".join(csp2_ordered)+"\n")
            
            seq_records = [SeqRecord(Seq(''.join(row)), id=query,description='') for query,row in alignment_df.iterrows()]
            
            alignment = MultipleSeqAlignment(seq_records)
            AlignIO.write(alignment,raw_alignment,"fasta")

            with open(log_file,"a+") as log:
                log.write("Done!\n")
                log.write(f"\t- Saved alignment of {snp_count} SNPs to {raw_alignment}\n")
                log.write(f"\t- Saved ordered loc list to {raw_loclist}\n")
            
            if max_missing == float(100):
                locs_pass_missing = csp2_ordered
                preserved_alignment = alignment
                
                AlignIO.write(preserved_alignment,preserved_alignment_file,"fasta")
                with open(preserved_loclist,"w+") as loclist:
                    loclist.write("\n".join(csp2_ordered)+"\n")
                
                with open(log_file,"a+") as log:
                    log.write("Skipping SNP preservation step...\n")
                    log.write(f"\t- Saved duplicate alignment to {preserved_alignment_file}\n")
                    log.write(f"\t- Saved duplicate ordered loc list to {preserved_loclist}\n")
            else:
                with open(log_file,"a+") as log:
                    log.write(f"\t- Preserving SNPs with at most {max_missing}% missing data...\n")
                
                # Parse missing data
                locs_pass_missing = list(set(locus_coverage_df[locus_coverage_df['Missing_Ratio'] <= max_missing]['Ref_Loc']))
                
                if len(locs_pass_missing) == 0:
                    with open(log_file,"a+") as log:
                        log.write(f"\t- Of {snp_count} SNPs, no SNPs pass the {max_missing}% missing data threshold...\n")
                        log.write("-------------------------------------------------------\n\n")
                else:
                    preserved_alignment_df = alignment_df[locs_pass_missing].copy()
                    
                    preserved_ordered = preserved_alignment_df.columns
                    with open(preserved_loclist,"w+") as loclist:
                        loclist.write("\n".join(preserved_ordered)+"\n")
                    
                    seq_records = [SeqRecord(Seq(''.join(row)), id=query,description='') for query,row in preserved_alignment_df.iterrows()]
                    preserved_alignment = MultipleSeqAlignment(seq_records)
                    AlignIO.write(preserved_alignment,preserved_alignment_file,"fasta")
                    with open(log_file,"a+") as log:
                        log.write(f"\t- Of {snp_count} SNPs, {len(locs_pass_missing)} SNPs pass the {max_missing}% missing data threshold...\n")
                        log.write(f"\t- Saved preserved alignment to {preserved_alignment_file}\n")
                        log.write(f"\t- Saved preserved ordered loc list to {preserved_loclist}\n")
                        log.write("-------------------------------------------------------\n\n")

    with open(log_file,"a+") as log:
        log.write("Processing pairwise comparisons files...")

    # Get pairwise comparisons between all pass_qc_isolates and reference_id
    pairwise_combinations = [sorted(x) for x in list(combinations([reference_id] + pass_qc_isolates, 2))]

    if snp_count == 0:
        pairwise_df = pd.DataFrame([(pairwise[0], pairwise[1], 0,np.nan) for pairwise in pairwise_combinations],columns = ['Query_1','Query_2','SNP_Distance','SNPs_Cocalled'])
        preserved_pairwise_df = pairwise_df.copy()
        
        with open(raw_pairwise,"w") as f:
            pairwise_df.to_csv(f, sep="\t", index=False)
        with open(preserved_pairwise,"w") as f:
            preserved_pairwise_df.to_csv(f, sep="\t", index=False)
        
        # Create matrix
        idx = sorted(set(pairwise_df['Query_1']).union(pairwise_df['Query_2']))
        mirrored_distance_df = pairwise_df.pivot(index='Query_1', columns='Query_2', values='SNP_Distance').reindex(index=idx, columns=idx).fillna(0, downcast='infer').pipe(lambda x: x+x.values.T).applymap(lambda x: format(x, '.0f'))
        mirrored_distance_df.index.name = ''
        with open(raw_matrix,"w") as f:
            mirrored_distance_df.to_csv(f,sep="\t")
        with open(preserved_matrix,"w") as f:
            mirrored_distance_df.to_csv(f,sep="\t")

    else:
        raw_distance_results = parallelAlignment(alignment)
        raw_pairwise_df = pd.DataFrame(raw_distance_results, columns=['Query_1', 'Query_2', 'SNP_Distance', 'SNPs_Cocalled'])
        with open(raw_pairwise,"w") as f:
            raw_pairwise_df.to_csv(f, sep="\t", index=False)

        if len(locs_pass_missing) == snp_count:
            preserved_pairwise_df = raw_pairwise_df.copy()
            with open(preserved_pairwise,"w") as f:
                preserved_pairwise_df.to_csv(f, sep="\t", index=False)
        elif len(locs_pass_missing) == 0:
            preserved_pairwise_df = pd.DataFrame([(pairwise[0], pairwise[1], 0,np.nan) for pairwise in pairwise_combinations],columns = ['Query_1','Query_2','SNP_Distance','SNPs_Cocalled'])
            with open(preserved_pairwise,"w") as f:
                preserved_pairwise_df.to_csv(f, sep="\t", index=False)        
        else:
            preserved_distance_results = parallelAlignment(preserved_alignment)
            preserved_pairwise_df = pd.DataFrame(preserved_distance_results, columns=['Query_1', 'Query_2', 'SNP_Distance', 'SNPs_Cocalled'])
            with open(preserved_pairwise,"w") as f:
                preserved_pairwise_df.to_csv(f, sep="\t", index=False)    
        
        # Create matrix
        idx = sorted(set(raw_pairwise_df['Query_1']).union(raw_pairwise_df['Query_2']))
        mirrored_distance_df = raw_pairwise_df.pivot(index='Query_1', columns='Query_2', values='SNP_Distance').reindex(index=idx, columns=idx).fillna(0, downcast='infer').pipe(lambda x: x+x.values.T).applymap(lambda x: format(x, '.0f'))
        mirrored_distance_df.index.name = ''
        with open(raw_matrix,"w") as f:
            mirrored_distance_df.to_csv(f,sep="\t")
        
        idx = sorted(set(preserved_pairwise_df['Query_1']).union(preserved_pairwise_df['Query_2']))
        mirrored_distance_df = preserved_pairwise_df.pivot(index='Query_1', columns='Query_2', values='SNP_Distance').reindex(index=idx, columns=idx).fillna(0, downcast='infer').pipe(lambda x: x+x.values.T).applymap(lambda x: format(x, '.0f'))
        mirrored_distance_df.index.name = ''
        with open(preserved_matrix,"w") as f:
            mirrored_distance_df.to_csv(f,sep="\t")

    # Clean up pybedtools temp
    helpers.cleanup(verbose=False,remove_all = False)

    end_time = time.time()
    with open(log_file,"a+") as log:
        log.write("Done!\n")
        if snp_count == 0:
            log.write(f"\t- No SNPs detected, zeroed pairwise distance files saved to {raw_pairwise}/{preserved_pairwise}/{raw_matrix}/{preserved_matrix}\n")
        else:
            log.write(f"\t- Saved raw pairwise distances to {raw_pairwise}\n")
            log.write(f"\t- Saved raw pairwise matrix to {raw_matrix}\n")

            if max_missing == float(100):
                log.write("Skipped SNP preservation step...\n")
                log.write(f"\t- Saved duplicated preserved pairwise distances to {preserved_pairwise}\n")
                log.write(f"\t- Saved duplicated preserved pairwise matrix to {preserved_matrix}\n")
            elif len(locs_pass_missing) == 0:
                log.write(f"\t- No SNPs passed the {max_missing}% missing data threshold, zeroed pairwise distance files saved to {preserved_pairwise}/{preserved_matrix}\n")
            else:
                log.write(f"\t- Saved preserved pairwise distances to {preserved_pairwise}\n")
                log.write(f"\t- Saved preserved pairwise matrix to {preserved_matrix}\n")
        log.write(f"Total Time: {end_time - start_time:.2f} seconds\n")
        log.write("-------------------------------------------------------\n\n")

except:
    run_failed = True
    print("Exception occurred:\n", traceback.format_exc())
finally:
    helpers.cleanup(verbose=False, remove_all=False)
    if temp_dir != "":
        shutil.rmtree(temp_dir)
    if run_failed:
        sys.exit(1)