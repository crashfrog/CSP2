#!/usr/bin/env python3

import sys
import os
import pandas as pd
import datetime
from pybedtools import BedTool,helpers
import concurrent.futures
import time
import uuid
import traceback
import shutil
import argparse
from csp2_common import (safe_int, safe_float, fetchHeaders, processBED, processSNPs, swapHeader, parseSNPDiffs, calculate_total_length, filterSNPs)


    


    


 

def screenSNPDiffs(snpdiffs_file,trim_name, min_cov, min_len, min_iden, ref_edge, query_edge, density_windows, max_snps,ref_ids):
    
    if temp_dir != "":
        helpers.set_tempdir(temp_dir)

    screen_start_time = time.time()

    # Set CSP2 variables to NA
    csp2_screen_snps = purged_length = purged_identity = purged_invalid = purged_indel = purged_lengthIdentity = purged_duplicate = purged_het = purged_density = filtered_ref_edge = filtered_query_edge = filtered_both_edge = "NA"

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
    if ref_ids == []:
        snpdiffs_orientation = 1
        query_id = header_query
        reference_id = header_ref
    elif (header_query not in ref_ids) and (header_ref in ref_ids):
        snpdiffs_orientation = 1
        query_id = header_query
        reference_id = header_ref
    elif (header_query in ref_ids) and (header_ref not in ref_ids):
        snpdiffs_orientation = -1
        query_id = header_ref
        reference_id = header_query
        header_data = swapHeader(header_data)
    else:
        snpdiffs_orientation = 2
        query_id = header_query
        reference_id = header_ref        

    # Establish log file
    log_file = f"{log_dir}/{query_id}__vs__{reference_id}.log"
    with open(log_file,"w+") as log:
        log.write("Screening Analysis\n")
        log.write(f"Query Isolate: {query_id}\n")
        log.write(f"Reference Isolate: {reference_id}\n")
        log.write(str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))+"\n")
        log.write("-------------------------------------------------------\n\n")
        if ref_ids == []:
            log.write("\t- No explicit references set, processing in forward orientation\n")
            log.write("-------------------------------------------------------\n\n")    
        elif snpdiffs_orientation == 1:
            log.write("\t- SNPDiffs file is in the forward orientation\n")
            log.write("-------------------------------------------------------\n\n")
        elif snpdiffs_orientation == -1:
            log.write("\t- SNPDiffs file is in the reverse orientation\n")
            log.write("-------------------------------------------------------\n\n")
        else:
            snpdiffs_orientation = 1
            log.write("\t- SNPDiffs file not contain a reference and non-reference sample, processing in forward orientation\n")
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
        except:
            with open(log_file,"a+") as log:
                log.write(f"Error reading BED/SNP data from file: {snpdiffs_file}")
            run_failed = True
            sys.exit(f"Error reading BED/SNP data from file: {snpdiffs_file}")
            
        ##### 03: Filter genome overlaps #####
        with open(log_file,"a+") as log:
            log.write("Step 2: Filtering for short overlaps and low percent identity...")
        
        good_bed_df = bed_df[(bed_df['Ref_Aligned'] >= min_len) & (bed_df['Perc_Iden'] >= min_iden)].copy()
        
        if good_bed_df.shape[0] == 0:
            screen_category = "Low_Quality_Coverage"
            query_percent_aligned = raw_query_percent_aligned
            reference_percent_aligned = raw_ref_percent_aligned
            with open(log_file,"a+") as log:
                log.write(f"\n\t- After filtering based on --min_len ({min_len}) and --min_iden ({min_iden}) , no valid alignments remain...Screen halted...\n")
                log.write("-------------------------------------------------------\n\n")

                
        else:
            # Create a BED file for alignments that pass basic QC
            good_query_bed_df = good_bed_df[['Query_Contig','Query_Start','Query_End']].copy()
            good_reference_bed_df = good_bed_df[['Ref_Contig','Ref_Start','Ref_End']].copy()
            
            good_query_aligned = calculate_total_length(BedTool.from_dataframe(good_query_bed_df).sort().merge())
            good_reference_aligned = calculate_total_length(BedTool.from_dataframe(good_reference_bed_df).sort().merge())
            
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
                    
                    # Write filtered SNP data to file
                    snp_file = log_file.replace(".log","_SNPs.tsv")
                    with open(snp_file,"w") as f:
                        filtered_snp_df.to_csv(f, sep="\t", index=False)
                    
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
                   
                    with open(log_file,"a+") as log:
                        log.write("Done!\n")
                        log.write(f"\t- {csp2_screen_snps} SNPs detected between {query_id} and {reference_id} after filtering\n")
                        log.write(f"\t- SNP data saved to {snp_file}\n")
                        log.write("-------------------------------------------------------\n\n")
    
    screen_end_time = time.time()
    helpers.cleanup(verbose=False, remove_all=False)

    with open(log_file,"a+") as log:
        log.write(f"Screening Time: {screen_end_time - screen_start_time:.2f} seconds\n")
    
    # Clean up pybedtools temp
    helpers.cleanup(verbose=False, remove_all=False)
    
    return [str(item) for item in [query_id,reference_id,screen_category,csp2_screen_snps,
            f"{query_percent_aligned:.2f}",f"{reference_percent_aligned:.2f}",
            query_contigs,query_bases,reference_contigs,reference_bases,
            raw_snps,purged_length,purged_identity,purged_lengthIdentity,purged_invalid,purged_indel,purged_duplicate,purged_het,purged_density,
            filtered_query_edge,filtered_ref_edge,filtered_both_edge,
            kmer_similarity,shared_kmers,query_unique_kmers,reference_unique_kmers,
            mummer_gsnps,mummer_gindels]]

# Read in arguments
global run_failed
run_failed = False

parser = argparse.ArgumentParser()
parser.add_argument("--snpdiffs_file", help="Path to the file containing SNP diffs")
parser.add_argument("--log_dir", help="Path to the log directory")
parser.add_argument("--min_cov", default=85, type=float, help="Minimum coverage")
parser.add_argument("--min_len", default=500,type=int, help="Minimum length")
parser.add_argument("--min_iden", default=99,type=float, help="Minimum identity")
parser.add_argument("--ref_edge", default=150,type=int, help="Reference edge")
parser.add_argument("--query_edge", default=150,type=int, help="Query edge")
parser.add_argument("--density_windows",default="1000,125,15", help="Density windows (comma-separated)")
parser.add_argument("--max_snps", default="3,2,1",help="Maximum SNPs (comma-separated)")
parser.add_argument('--trim_name', nargs='?', const="", default="", type=str, help='Trim name')
parser.add_argument("--output_file", help="Output file")
parser.add_argument("--ref_id", help="Reference IDs file")
parser.add_argument("--tmp_dir",default="", help="TMP dir")

args = parser.parse_args()

snpdiffs_list = [line.strip() for line in open(args.snpdiffs_file, 'r')]
snpdiffs_list = [line for line in snpdiffs_list if line]
for snpdiffs_file in snpdiffs_list:
    if not os.path.exists(snpdiffs_file):
        run_failed = True
        sys.exit("Error: File does not exist: " + snpdiffs_file)

snpdiffs_list = list(set(snpdiffs_list))

log_dir = os.path.normpath(os.path.abspath(args.log_dir))

min_cov = args.min_cov
min_len = args.min_len
min_iden = args.min_iden

ref_edge = args.ref_edge
query_edge = args.query_edge

input_density = args.density_windows
input_maxsnps = args.max_snps

if input_density == "0":
    density_windows = []
    max_snps = []
else:
    density_windows = [int(x) for x in args.density_windows.split(",")]
    max_snps = [int(x) for x in args.max_snps.split(",")]
assert len(density_windows) == len(max_snps)

trim_name = args.trim_name

output_file = os.path.abspath(args.output_file)

if os.stat(args.ref_id).st_size == 0:
    ref_ids = []
else:
    ref_ids = [line.strip() for line in open(args.ref_id, 'r')]

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

try:    
    with concurrent.futures.ProcessPoolExecutor() as executor:
        results = [executor.submit(screenSNPDiffs,snp_diff_file,trim_name, min_cov, min_len, min_iden, ref_edge, query_edge, density_windows, max_snps,ref_ids) for snp_diff_file in snpdiffs_list]

    # Clean up pybedtools temp
    helpers.cleanup(verbose=False,remove_all = False)

    # Combine results into a dataframe
    output_columns = ['Query_ID','Reference_ID','Screen_Category','CSP2_Screen_SNPs',
                'Query_Percent_Aligned','Reference_Percent_Aligned',
                'Query_Contigs','Query_Bases','Reference_Contigs','Reference_Bases',
                'Raw_SNPs','Purged_Length','Purged_Identity','Purged_LengthIdentity','Purged_Invalid','Purged_Indel','Purged_Duplicate','Purged_Het','Purged_Density',
                'Filtered_Query_Edge','Filtered_Ref_Edge','Filtered_Both_Edge',
                'Kmer_Similarity','Shared_Kmers','Query_Unique_Kmers','Reference_Unique_Kmers',
                'MUMmer_gSNPs','MUMmer_gIndels']

    results_df = pd.DataFrame([item.result() for item in results], columns = output_columns)
    with open(output_file,"w") as f:
        results_df.to_csv(f, sep="\t", index=False)
except:
    run_failed = True
    print("Exception occurred:\n", traceback.format_exc())
finally:
    helpers.cleanup(verbose=False, remove_all=False)
    if temp_dir != "":
        shutil.rmtree(temp_dir)
    if run_failed:
        sys.exit(1)


