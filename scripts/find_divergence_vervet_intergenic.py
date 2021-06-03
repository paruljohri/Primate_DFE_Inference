


#This is a Python script that calculates divergence rate values for vervet INTERGENIC REGIONS specifically. 
#Whereas most scripts that have separate "find x for chimps" and "find x for vervets" versions
#have none or extremely minor differences, these chimp and vervet divergence calculation 
#scripts are substantially different due to how they go about calculating whether or 
#not a site is diverged in addition to the vervet portion being split into separate scripts for 
#intergenic and exonic regions.

#This script is currently designed to be run locally and will require path changes on the user's part to the variables:
#'divergence_file' with input files "pull_out_intergenic_aligned_sites.txt"
#'region_data_file' with input file "248_homo_exons_reference_with_intergenic_lengths.txt"
#'result_file' with output file "vervet_all_exons_divergence.txt"




#This section reads in information from the "pull_out_intergenic_aligned_sites.txt" file, a file 
#created by scanning through all of the post alignment files with the ancestral allele already 
#calculated and the callability of that site already calculated. Any sites that fall within any 
#the 5' or 3' intergenic regions flanking the 248 identified homologous vervet exons are output. 
#With this information, we can know and an account for what percentage of an intergenic region is 
#lacking alignment. Reading in this information, it puts the vervet allele, calculated ancestral 
#allele, and whether or not the site is masked and puts that information in a dictionary for future 
#calculation of divergence.

divergence_file = 'C:/Users/Kellen/Downloads/pull_out_intergenic_aligned_sites.txt'
d_div = {}

with open(divergence_file, 'r') as f:
    skip_first_line = 0
    for line in f:
        if skip_first_line == 0:
            skip_first_line = 1
            continue
        
        split = line.strip('\n').split()
        total_position = split[0]
        v_nuc = split[1]
        v_h_nuc = split[2]
        v_h_m_nuc = split[3]
        ancestral_nuc = split[4]
        is_diverged = split[5]
        is_callable = split[6]
        
        d_div[total_position] = [v_nuc, ancestral_nuc, is_callable]
                




#This function takes in the start and end position of the intergenic regions one wishes to find the divergence 
#for and scans through the previously created dictionary. All site no present in the dictionary was not present 
#in the alignment, and is thus not counted towards the denominator when calculating divergence. For sites that 
#are present, it tallies up the number of diverged sites and then obtains a final divergence value by dividing
#that tally by the total number of callable sites in the region. Only sites that are callable are considered, 
#so positions that may be diverged but are not callable are not considered in the final calculation. Regions 
#with no callable sites are given a divergence value of "NA" in the output file.

def get_divergence(chr_num, start_pos, end_pos):
    possible_nuc = ['A', 'a', 'T', 't', 'G', 'g', 'C', 'c']
    
    div_numerator = 0
    div_denominator = 0
    
    i = int(start_pos)
    while i <= int(end_pos):
        s = chr_num + '_' + str(i)
        at_site = d_div.get(s, 'NA')
        if at_site != 'NA':
            v_nuc = at_site[0]
            ancestral_nuc = at_site[1]
            is_callable = at_site[2]
            
            if is_callable == 'Y':
                div_denominator += 1
                if ancestral_nuc in possible_nuc:
                    if v_nuc.upper() != ancestral_nuc.upper():
                        div_numerator += 1
        i += 1
    
    if div_denominator == 0:
        divergence = 'NA'
    else:
        divergence = div_numerator / div_denominator
        divergence = round(divergence * 10**9) / 10**9
    
    return (divergence, div_denominator)




#This section goes through all the 5' and 3' intergenic regions flanking the established homologous 
#vervet exons and calculates the hard divergence value, soft divergence value, and number of callable '
#sites for a particular region. 

"""For HOMO regions"""
region_data_file = 'C:/Users/Kellen/Downloads/primate_dfe_project/248_homo_exons_reference_with_intergenic_lengths.txt'
with open(region_data_file, 'r') as g:
    region_data = g.readlines()[1 : ]

    result_file = 'C:/Users/Kellen/Downloads/primate_dfe_project/chimp_and_vervet_data/vervet_homo_intergenic_divergence.txt'
    results = open(result_file, 'w+')
    results.write('name' + '\t' + 'exon' + '\t' + 'chr' + '\t' + 'strand' + '\t' + 'intergenic_5p_start' + '\t' + 'intergenic_5p_end' + '\t' + 'intergenic_3p_start' + '\t' + 'intergenic_3p_end' + '\t' + 'div_5p' + '\t' + '5p_callability' + '\t' + 'div_3p' + '\t' + '3p_callability' + '\n')
    
    total_num_lines = len(region_data)
    for w in range(0, total_num_lines):
        line = region_data[w].strip('\n')
        split = line.split()
        name = split[32]
        exon = split[33]
        chr_num = split[34]
        strand = split[35]
        start_5p = split[41]
        end_5p = split[42]
        len_5p = split[43]
        start_3p = split[44]
        end_3p = split[45]
        len_3p = split[46]
        
        if start_5p == 'NA' or end_5p == 'NA' or len_5p == 'NA':
            div_5p = 'NA'
            div_denominator_5p = 'NA'
        else:
            div_5p, div_denominator_5p = get_divergence(chr_num, int(start_5p), int(end_5p))

        if start_3p == 'NA' or end_3p == 'NA' or len_3p == 'NA':
            div_3p = 'NA'
            div_denominator_3p = 'NA'
        else:
            div_3p, div_denominator_3p = get_divergence(chr_num, int(start_3p), int(end_3p))   

        results.write(name + '\t' + exon + '\t' + chr_num + '\t' + strand + '\t' + start_5p + '\t' + end_5p + '\t' + start_3p + '\t' + end_3p + '\t' + str(div_5p) + '\t' + str(div_denominator_5p) + '\t' + str(div_3p) + '\t' + str(div_denominator_3p) + '\n') 
    results.close()

