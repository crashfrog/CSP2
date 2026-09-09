#!/usr/bin/env python3
"""Shared .snpdiffs parsing + SNP QC filtering for CSP2.

Extracted verbatim from screenSNPDiffs.py / runSNPPipeline.py, which each carried a
byte-for-byte (modulo whitespace) copy of these helpers. Behaviour is unchanged; the
only signature change is filterSNPs(), which now takes temp_dir explicitly instead of
reading a module global (so it works when imported from either caller).
"""
import pandas as pd
from pybedtools import BedTool, helpers


def safe_int(val):
    try:
        return int(val) if val != "NA" else 0
    except:
        return 0


def safe_float(val):
    try:
        return float(val) if val != "NA" else 0
    except:
        return 0


def fetchHeaders(snpdiffs_file):
    
    with open(snpdiffs_file, 'r') as file:
        top_line = file.readline().strip().split('\t')[1:]

    header_cols = [item.split(':')[0] for item in top_line]
    header_vals = [item.split(':')[1] for item in top_line]
    
    header_data = pd.DataFrame(columns = header_cols)
    header_data.loc[0] = header_vals
    header_data.loc[:, 'File_Path'] = snpdiffs_file
      
    return header_data


def processBED(bed_rows,snpdiffs_orientation):
    
    bed_columns = ['Ref_Contig','Ref_Start','Ref_End','Ref_Length','Ref_Aligned',
                   'Query_Contig','Query_Start','Query_End','Query_Length','Query_Aligned',
                   'Perc_Iden']
    
    reverse_columns = ['Query_Contig','Query_Start','Query_End','Query_Length','Query_Aligned',
                     'Ref_Contig','Ref_Start','Ref_End','Ref_Length','Ref_Aligned',
                     'Perc_Iden']
    
    int_columns = ['Ref_Start', 'Ref_End', 'Ref_Length', 'Ref_Aligned',
                'Query_Start', 'Query_End', 'Query_Length', 'Query_Aligned']

    float_columns = ['Perc_Iden']
    
    if len(bed_rows) > 0:
        
        bed_df = pd.DataFrame(bed_rows, columns=bed_columns)
        
        # Swap columns if reversed
        if snpdiffs_orientation == -1:
            bed_df = bed_df[reverse_columns].copy()
            bed_df.columns = bed_columns
            
        # Remove any rows where Query_Contig or Ref_Contig == "." (Unaligned)
        covered_bed_df = bed_df[(bed_df['Ref_Start'] != ".") & (bed_df['Query_Start'] != ".")].copy()
        
        if covered_bed_df.shape[0] > 0:
            for col in int_columns:
                covered_bed_df.loc[:, col] = covered_bed_df.loc[:, col].astype(float).astype(int)
            for col in float_columns:
                covered_bed_df.loc[:, col] = covered_bed_df.loc[:, col].astype(float)
            return covered_bed_df        
        else:
            return pd.DataFrame(columns=bed_columns)           
    else:
        return pd.DataFrame(columns=bed_columns)


def processSNPs(snp_rows,snpdiffs_orientation):
    
    snp_columns = ['Ref_Contig','Start_Ref','Ref_Pos',
                'Query_Contig','Start_Query','Query_Pos',
                'Ref_Loc','Query_Loc',
                'Ref_Start','Ref_End',
                'Query_Start','Query_End',
                'Ref_Base','Query_Base',
                'Dist_to_Ref_End','Dist_to_Query_End',
                'Ref_Aligned','Query_Aligned',
                'Query_Direction','Perc_Iden','Cat']
        
    reverse_columns = ['Query_Contig','Start_Query','Query_Pos',
                    'Ref_Contig','Start_Ref','Ref_Pos',
                    'Query_Loc','Ref_Loc',
                    'Query_Start','Query_End',
                    'Ref_Start','Ref_End',
                    'Query_Base','Ref_Base',
                    'Dist_to_Query_End','Dist_to_Ref_End',
                    'Query_Aligned','Ref_Aligned',
                    'Query_Direction','Perc_Iden','Cat']
    
    return_columns = ['Ref_Contig','Start_Ref','Ref_Pos',
            'Query_Contig','Start_Query','Query_Pos',
            'Ref_Loc','Query_Loc',
            'Ref_Start','Ref_End',
            'Query_Start','Query_End',
            'Ref_Base','Query_Base',
            'Dist_to_Ref_End','Dist_to_Query_End',
            'Ref_Aligned','Query_Aligned',
            'Perc_Iden','Cat']
    
    reverse_complement = {'A':'T','T':'A','G':'C','C':'G',
                          'a':'T','t':'A','c':'G','g':'C'}
    
    # Columns to convert to integer
    int_columns = ['Start_Ref', 'Ref_Pos', 'Start_Query', 'Query_Pos',
                'Dist_to_Ref_End', 'Dist_to_Query_End', 'Ref_Aligned', 'Query_Aligned']

    # Columns to convert to float
    float_columns = ['Perc_Iden']
    
    if len(snp_rows) > 0:
        snp_df = pd.DataFrame(snp_rows, columns= snp_columns).copy()
        
        if snpdiffs_orientation == -1:
            snp_df = snp_df[reverse_columns].copy()
            snp_df.columns = snp_columns
            
            # Replace Query_Base and Reference_Base with reverse complement if Query_Direction is -1 and base is in ['A','T','G','C','a','c','t','g']
            snp_df.loc[snp_df['Query_Direction'] == '-1','Query_Base'] = snp_df.loc[snp_df['Query_Direction'] == '-1','Query_Base'].apply(lambda x: reverse_complement[x] if x in reverse_complement else x)
            snp_df.loc[snp_df['Query_Direction'] == '-1','Ref_Base'] = snp_df.loc[snp_df['Query_Direction'] == '-1','Ref_Base'].apply(lambda x: reverse_complement[x] if x in reverse_complement else x)

            
        for col in int_columns:
            snp_df.loc[:, col] = snp_df.loc[:, col].astype(float).astype(int)
        for col in float_columns:
            snp_df.loc[:, col] = snp_df.loc[:, col].astype(float)
        
    else:
        snp_df = pd.DataFrame(columns = return_columns)

    return snp_df[return_columns]


def swapHeader(header_data):
            
    raw_header_cols = [x for x in header_data.columns]
    reverse_header_cols = [item.replace('Reference', 'temp').replace('Query', 'Reference').replace('temp', 'Query') for item in raw_header_cols]
    reversed_header_data = header_data[reverse_header_cols].copy()
    reversed_header_data.columns = raw_header_cols

    return reversed_header_data


def parseSNPDiffs(snpdiffs_file,snpdiffs_orientation):

    bed_rows = []
    snp_rows = []
    
    with open(snpdiffs_file, 'r') as file:
        lines = file.readlines()

    for line in lines:
        if line[0:2] == "#\t":
            pass
        elif line[0:3] == "##\t":
            bed_rows.append(line.strip().split("\t")[1:])
        else:
            snp_rows.append(line.strip().split("\t"))

    bed_df = processBED(bed_rows,snpdiffs_orientation)
    snp_df = processSNPs(snp_rows,snpdiffs_orientation)
    return (bed_df,snp_df)


def calculate_total_length(bedtool):
    return sum(len(interval) for interval in bedtool)


def filterSNPs(raw_snp_df,bed_df,log_file, min_len, min_iden, ref_edge, query_edge, density_windows, max_snps, temp_dir=""):
    
    if temp_dir != "":
        helpers.set_tempdir(temp_dir)
    
    # Grab raw data    
    total_snp_count = raw_snp_df.shape[0]
    
    # Get unique SNPs relative to the reference genome
    unique_ref_snps = raw_snp_df['Ref_Loc'].unique()
    unique_snp_count = len(unique_ref_snps)
    
    snp_tally_df = pd.DataFrame()
    
    with open(log_file,"a+") as log:
        log.write(f"\n\t- Raw SNP + indel count: {total_snp_count}\n")
        log.write(f"\n\t- Unique SNP positions in reference genome: {unique_snp_count}\n")
    
    # Set all sites to SNP
    raw_snp_df['Filter_Cat'] = "SNP"
    
    # Filter out SNPs based on --min_len and --min_iden
    reject_length = raw_snp_df.loc[(raw_snp_df['Ref_Aligned'] < min_len) & (raw_snp_df['Perc_Iden'] >= min_iden)].copy()
    if reject_length.shape[0] > 0:
        with open(log_file,"a+") as log:
            log.write(f"\t\t- Purged (Alignment Length): {reject_length.shape[0]}\n")
        reject_length['Filter_Cat'] = "Purged_Length"
        snp_tally_df = pd.concat([snp_tally_df,reject_length]).reset_index(drop=True)
        
    reject_iden = raw_snp_df.loc[(raw_snp_df['Ref_Aligned'] >= min_len) & (raw_snp_df['Perc_Iden'] < min_iden)].copy()
    if reject_iden.shape[0] > 0:
        with open(log_file,"a+") as log:
            log.write(f"\t\t- Purged (Alignment Identity): {reject_iden.shape[0]}\n")
        reject_iden['Filter_Cat'] = "Purged_Identity"
        snp_tally_df = pd.concat([snp_tally_df,reject_iden]).reset_index(drop=True)

    reject_lenIden = raw_snp_df.loc[(raw_snp_df['Ref_Aligned'] < min_len) & (raw_snp_df['Perc_Iden'] < min_iden)].copy()
    if reject_lenIden.shape[0] > 0:
        with open(log_file,"a+") as log:
            log.write(f"\t\t- Purged (Alignment Length + Identity): {reject_lenIden.shape[0]}\n")
        reject_lenIden['Filter_Cat'] = "Purged_LengthIdentity"
        snp_tally_df = pd.concat([snp_tally_df,reject_lenIden]).reset_index(drop=True)
    
    pass_filter = raw_snp_df.loc[(raw_snp_df['Ref_Aligned'] >= min_len) & (raw_snp_df['Perc_Iden'] >= min_iden)].copy().reset_index(drop=True)
    
    # Invalid processing
    reject_invalid = pass_filter[pass_filter['Cat'] == "Invalid"].copy()
    if reject_invalid.shape[0] > 0:
        with open(log_file,"a+") as log:
            log.write(f"\t\t- Purged (Invalid Base): {reject_invalid.shape[0]}\n")
        reject_invalid['Filter_Cat'] = "Purged_Invalid"
        snp_tally_df = pd.concat([snp_tally_df,reject_invalid]).reset_index(drop=True)
    pass_filter = pass_filter.loc[pass_filter['Cat'] != "Invalid"].copy()
    
    # Indel processing
    reject_indel = pass_filter[pass_filter['Cat'] == "Indel"].copy()
    if reject_indel.shape[0] > 0:
        with open(log_file,"a+") as log:
            log.write(f"\t\t- Purged (Indel): {reject_indel.shape[0]}\n")
        reject_indel['Filter_Cat'] = "Purged_Indel"
        snp_tally_df = pd.concat([snp_tally_df,reject_indel]).reset_index(drop=True)
    pass_filter = pass_filter.loc[pass_filter['Cat'] != "Indel"].copy()
    
    # Check for heterozygous SNPs
    check_heterozygous = pass_filter.groupby('Ref_Loc').filter(lambda x: x['Query_Base'].nunique() > 1)
    if check_heterozygous.shape[0] > 0:      
        reject_heterozygous = pass_filter.loc[pass_filter['Ref_Loc'].isin(check_heterozygous['Ref_Loc'])].copy()
        reject_heterozygous['Filter_Cat'] = "Purged_Heterozygous"
        with open(log_file,"a+") as log:
            log.write(f"\t\t- Purged (Heterozygotes): {reject_heterozygous.shape[0]}\n")  
        snp_tally_df = pd.concat([snp_tally_df,reject_heterozygous]).reset_index(drop=True)
        pass_filter = pass_filter.loc[~pass_filter['Ref_Loc'].isin(check_heterozygous['Ref_Loc'])].copy()
   
    # Check for duplicate SNPs and take the longest, best hit
    check_duplicates = pass_filter.groupby('Ref_Loc').filter(lambda x: x.shape[0] > 1)
    if check_duplicates.shape[0] > 0:
        reject_duplicate = pass_filter.loc[pass_filter['Ref_Loc'].isin(check_duplicates['Ref_Loc'])].copy()
        pass_filter = pass_filter.loc[~pass_filter['Ref_Loc'].isin(check_duplicates['Ref_Loc'])].copy()
                
        best_snp = reject_duplicate.groupby('Ref_Loc').apply(lambda x: x.sort_values(by=['Ref_Aligned', 'Perc_Iden'], ascending=[False, False]).head(1))
        pass_filter = pd.concat([pass_filter,best_snp]).reset_index(drop=True)

        dup_snps = reject_duplicate[~reject_duplicate.apply(lambda x: x in best_snp, axis=1)]        
        dup_snps['Filter_Cat'] = "Purged_Duplicate"
        
        snp_tally_df = pd.concat([snp_tally_df,dup_snps]).reset_index(drop=True)
                        
        with open(log_file,"a+") as log:
            log.write(f"\t\t- Purged (Duplicates): {dup_snps.shape[0]}\n")    
    
    # Assert that Ref_Loc and Query_Loc are unique in pass_filter
    helpers.cleanup(verbose=False,remove_all = False)
    assert pass_filter['Ref_Loc'].nunique() == pass_filter.shape[0]
    assert pass_filter['Query_Loc'].nunique() == pass_filter.shape[0]
    
    # Density filtering
    density_locs = []
    ref_locs = pass_filter['Ref_Loc'].tolist()
        
    if len(density_windows) == 0:
        with open(log_file,"a+") as log:
            log.write("\n\t- Density filtering disabled...\n")
    elif len(ref_locs) > 0:
        density_df = pd.DataFrame([item.split('/') for item in ref_locs], columns=['Ref_Contig','Ref_End'])
        density_df['Ref_Start'] = density_df['Ref_End'].astype(float).astype(int) - 1
        density_bed = BedTool.from_dataframe(density_df[['Ref_Contig','Ref_Start','Ref_End']])

        # For each density window, remove all SNPs that fall in a window with > max_snps
        for i in range(0,len(density_windows)):
            window_df = density_bed.window(density_bed,c=True, w=density_windows[i]).to_dataframe()
            problematic_windows = window_df[window_df['name'] > max_snps[i]].copy()
            if not problematic_windows.empty:
                temp_locs = []            
                for _, row in problematic_windows.iterrows():
                        purge_window_df = window_df[window_df['chrom'] == row['chrom']].copy()
                        purge_window_df['Dist'] = abs(purge_window_df['end'] - row['end'])
                        window_snps = purge_window_df.sort_values(by=['Dist'],ascending=True).head(row['name'])
                        temp_locs = temp_locs + ["/".join([str(x[0]),str(x[1])]) for x in list(zip(window_snps.chrom, window_snps.end))]
                density_locs.extend(list(set(temp_locs)))
    
    density_locs = list(set(density_locs))
    reject_density = pass_filter[pass_filter['Ref_Loc'].isin(density_locs)].copy() 
           
    if reject_density.shape[0] > 0:
        with open(log_file,"a+") as log:
            log.write(f"\t\t- Purged (Density): {reject_density.shape[0]}\n")
        reject_density['Filter_Cat'] = "Purged_Density"
        snp_tally_df = pd.concat([snp_tally_df,reject_density]).reset_index(drop=True)
        pass_filter = pass_filter[~pass_filter['Ref_Loc'].isin(density_locs)].copy()

    reject_query_edge = pass_filter[(pass_filter['Dist_to_Query_End'] < query_edge) & (pass_filter['Dist_to_Ref_End'] >= ref_edge)].copy()
    reject_ref_edge = pass_filter[(pass_filter['Dist_to_Ref_End'] < ref_edge) & (pass_filter['Dist_to_Query_End'] >= query_edge)].copy()
    reject_both_edge = pass_filter[(pass_filter['Dist_to_Query_End'] < query_edge) & (pass_filter['Dist_to_Ref_End'] < ref_edge)].copy()
    
    if reject_query_edge.shape[0] > 0:
        with open(log_file,"a+") as log:
            log.write(f"\t\t- Purged (Query Edge): {reject_query_edge.shape[0]}\n")
        reject_query_edge['Filter_Cat'] = "Filtered_Query_Edge"
        snp_tally_df = pd.concat([snp_tally_df,reject_query_edge]).reset_index(drop=True)
            
    if reject_ref_edge.shape[0] > 0:
        with open(log_file,"a+") as log:
            log.write(f"\t\t- Purged (Ref Edge): {reject_ref_edge.shape[0]}\n")
        reject_ref_edge['Filter_Cat'] = "Filtered_Ref_Edge"
        snp_tally_df = pd.concat([snp_tally_df,reject_ref_edge]).reset_index(drop=True)
            
    if reject_both_edge.shape[0] > 0:
        with open(log_file,"a+") as log:
            log.write(f"\t\t- Purged (Both Edge): {reject_both_edge.shape[0]}\n")
        reject_both_edge['Filter_Cat'] = "Filtered_Both_Edge"
        snp_tally_df = pd.concat([snp_tally_df,reject_both_edge]).reset_index(drop=True)
        
    pass_filter = pass_filter[(pass_filter['Dist_to_Query_End'] >= query_edge) & (pass_filter['Dist_to_Ref_End'] >= ref_edge)].copy()
    
    helpers.cleanup(verbose=False,remove_all = False)

    assert snp_tally_df.shape[0] + pass_filter.shape[0] == total_snp_count
    return_df = pd.concat([pass_filter,snp_tally_df]).reset_index(drop=True).sort_values(by=['Ref_Loc'])
        
    return return_df.drop(columns=['Cat']).rename({'Filter_Cat':'Cat'}, axis=1)
