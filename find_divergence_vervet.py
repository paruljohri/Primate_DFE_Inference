


#This is a Python script that calculates divergence rate values for vervets. 
#Whereas most scripts that have separate "find x for chimps" and "find x for vervets" versions
#have none or extremely minor differences, these chimp and vervet divergence calculation 
#scripts are substantially different due to how they go about calculating whether or 
#not a site is diverged.

#This script is currently designed to be run locally and will require path changes on the user's part to the variables:
#'mask_file' with input files 'chr' + chr_num + '.mask.fa'
#'divergence_file' with input files "all_sites_that_differ.txt"
#'region_data_file' with input file "248_homo_exons_reference_with_intergenic_lengths.txt"
#first instance of 'result_file' with output file "vervet_homo_exons_divergence.txt"
#second instance of 'result_file' with output file "vervet_homo_intergenic_divergence.txt"
#'exon_file' with input file "chlSab2.ensGenes.txt"
#third instance of 'result_file' with output file "vervet_all_exons_divergence.txt"
#'intergenic_file' with input file "vervet_all_intergenic_lengths.txt"
#fourth instance of 'result_file' with output file "vervet_all_intergenic_divergence.txt"




#This section reads in all of the data from the vervet mask files, allowing us to determine if 
#a site is callable. The number of callable sites is then used in the denominator of the divergence calculation.

all_chr_nums = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29']
all_mask_data = {}
for chr_num in all_chr_nums:
    mask_file = 'C:/Users/Kellen/Downloads/primate_dfe_project/chimp_and_vervet_data/vervet_mask/chr' + chr_num + '_mask.fasta'
    with open(mask_file) as f:
        mask_data = f.readlines()[1 : ]
    all_mask_data[chr_num] = mask_data[0]
mask_data = []




#This section reads in information from the "all_sites_that_differ.txt" file, a file 
#created by scanning through all of the post alignment files to find any vervet sites 
#where one or multiple of the alleles needed to determine the final ancestral allele 
#(vervet allele, vervet-human ancestral allele, or vervet-human-marmoset allele) 
#differs from the rest. With this information, this section then goes through that file,
#determines the ancestral allele, and then determines if the vervet allele is diverged at
#that site or not. Finally, it adds any of the diverged positions to a dictionary for 
#further calculation.

divergence_file = 'C:/Users/Kellen/Downloads/all_sites_that_differ_NEW.txt/all_sites_that_differ.txt'
d_div = {}
hard = ['A', 'T', 'G', 'C']
soft = ['a', 't', 'g', 'c']

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
        h_nuc = split[4]
        
        seq_a = v_h_nuc.upper()
        seq_b = v_nuc.upper()
        seq_c = v_h_m_nuc.upper()
        
#         Rules used to determine ancestral allele:
#         (*) uppercase when all 3 --(a), (b) and (c)-- sequences agree
#         (*) lowercase when:
#             - there is no ancestral sequence for the ancestral sequence, i.e. there are only two extant sequences in the alignment, but (a) and (b) agree.
#             - there is a gap in the sister sequence, but (a) and (c) agree.
#             - either (b) or (c) disagree with (a), but not both.
#         (*) N when both (b) and (c) disagree with (a)
#         (*) - (dash) when no there is no ancestral allele, this is a lineage-specific insertion
#         (*) . (dot) when there is no alignment, i.e. no data
            
        
        if seq_a == 'X' and seq_c == 'X':
            ancestral_nuc = '.'
        elif seq_a == '-' and seq_c == '-':
            ancestral_nuc = '-'
        elif seq_a == seq_b == seq_c:
            ancestral_nuc = seq_a.upper()
        elif seq_c == 'X' and seq_a == seq_b:
            ancestral_nuc = seq_a.lower()
        elif seq_b == '-' and seq_a == seq_c:
            ancestral_nuc = seq_a.lower()
        elif (seq_a == seq_b and seq_a != seq_c) or (seq_a == seq_c and seq_a != seq_b):
            ancestral_nuc = seq_a.lower()
        elif seq_a != seq_b and seq_a != seq_c:
            ancestral_nuc = 'N'

        if ancestral_nuc in hard:
            if v_nuc != ancestral_nuc:
                d_div[total_position] = 'H'
        elif ancestral_nuc in soft:
            if v_nuc.lower() != ancestral_nuc.lower():
                d_div[total_position] = 'S'
                




#This function takes in the start and end position of a region one wishes to find the divergence 
#for, tallies up the number of diverged sites, and then obtains a final divergence value by dividing
#that tally by the total number of callable sites in the region. Please note that this is really 
#just the raw number of callable sites, as positions that are diverged but said to be not callable
#based on the information from the mask file are still counted in the tally. 
#While I thought this was odd to not eliminate sites from being counted if they're not callable, 
#Susanne said it was acceptable to do it this way. Regions with no callable sites though are 
#given a divergence value of "NA".

def get_divergence(chr_num, start_pos, end_pos, true_denominator):
    hard_div_num = 0
    hard_div_denom = 0
    soft_div_num = 0 
    soft_div_denom = 0
    i = int(start_pos)
    while i <= int(end_pos):
        s = chr_num + '_' + str(i)
        #print(s)
        if d_div.get(s, 'NA') == 'H':
            hard_div_num += 1
            hard_div_denom += 1
            soft_div_num += 1
            soft_div_denom += 1
        elif d_div.get(s, 'NA') == 'S':
            soft_div_num += 1
            soft_div_denom += 1
        elif d_div.get(s, 'NA') == 'NA':
            hard_div_denom += 1
            soft_div_denom += 1
        i += 1
    
    
    if true_denominator == 0:
        hard_div = 'NA'
        soft_div = 'NA'
    else:
        hard_div = float(hard_div_num) / float(true_denominator)
        soft_div = float(soft_div_num) / float(true_denominator)
        
    return (hard_div, hard_div_denom, soft_div, soft_div_denom)




#This section goes through all the established homologous vervet exons and calculates the 
#hard divergence value, soft divergence value, and number of callable sites for a 
#particular region. This can be performed for either exons on intergenic regions 
#by changing the "find_for" variable to either 'exons' or 'intergenic'. 

"""For HOMO regions"""
find_for = 'exons'
# find_for = 'intergenic'
print('ran for: ' + find_for)

region_data_file = 'C:/Users/Kellen/Downloads/primate_dfe_project/248_homo_exons_reference_with_intergenic_lengths.txt'
# region_data_file = 'C:/Users/Kellen/Downloads/primate_dfe_project/FILTERED_homo_exons_reference_all_info.txt'
with open(region_data_file, 'r') as g:
    region_data = g.readlines()[1:]

if find_for == 'exons':
    result_file = 'C:/Users/Kellen/Downloads/primate_dfe_project/chimp_and_vervet_data/vervet_homo_exons_divergence.txt'
    results = open(result_file, 'w+')
    results.write('name' + '\t' + 'exon' + '\t' + 'chr' + '\t' + 'strand' + '\t' + 'start' + '\t' + 'end' + '\t' + 'hard_div' + '\t' + 'soft_div' + '\t' + 'callability' + '\n')

    total_num_lines = len(region_data)
    for w in range(0, total_num_lines):
        line = region_data[w].strip('\n')
        split = line.split()
        name = split[32]
        exon = split[33]
        chr_num = split[34]
        strand = split[35]
        start_pos = int(split[38])
        end_pos = int(split[39])
        
        if chr_num == 'chrX' or chr_num == 'chrY' or chr_num == 'chrUn' or chr_num == 'chrM' or 'random' in chr_num:
            continue
        else:
            chr_num_stripped = chr_num.strip('chr')
            cut_out = all_mask_data[chr_num_stripped][start_pos - 1 : end_pos]
            true_denominator = cut_out.count('0')
            
            t_div = get_divergence(chr_num, start_pos, end_pos, true_denominator)
            hard_div = t_div[0]
            hard_div_denom = t_div[1]
            soft_div = t_div[2]
            soft_div_denom = t_div[3]
            
            results.write(name + '\t' + exon + '\t' + chr_num + '\t' + strand + '\t' + str(start_pos) + '\t' + str(end_pos) + '\t' + str(hard_div) + '\t' + str(soft_div) + '\t' + str(true_denominator) + '\n')
    results.close()
    
if find_for == 'intergenic':
    result_file = 'C:/Users/Kellen/Downloads/primate_dfe_project/chimp_and_vervet_data/vervet_homo_intergenic_divergence.txt'
    results = open(result_file, 'w+')
    results.write('name' + '\t' + 'exon' + '\t' + 'chr' + '\t' + 'strand' + '\t' + 'intergenic_5p_start' + '\t' + 'intergenic_5p_end' + '\t' + 'intergenic_3p_start' + '\t' + 'intergenic_3p_end' + '\t' + 'hard_div_5p' + '\t' + 'soft_div_5p' + '\t' + '5p_callability' + '\t' + 'hard_div_3p' + '\t' + 'soft_div_3p' + '\t' + '3p_callability' + '\n')

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
            t_div_5p = ['NA', 'NA', 'NA', 'NA']
        else:
            chr_num_stripped = chr_num.strip('chr')
            cut_out = all_mask_data[chr_num_stripped][int(start_5p) - 1 : int(end_5p)]
            true_denominator_5p = cut_out.count('0')
            t_div_5p = get_divergence(chr_num, int(start_5p), int(end_5p), true_denominator_5p)
            
        hard_div_5p = t_div_5p[0]
        hard_div_denom_5p = t_div_5p[1]
        soft_div_5p = t_div_5p[2]
        soft_div_denom_5p = t_div_5p[3]

        if start_3p == 'NA' or end_3p == 'NA' or len_3p == 'NA':
            t_div_3p = ['NA', 'NA', 'NA', 'NA']
        else:
            chr_num_stripped = chr_num.strip('chr')
            cut_out = all_mask_data[chr_num_stripped][int(start_3p) - 1 : int(end_3p)]
            true_denominator_3p = cut_out.count('0')
            t_div_3p = get_divergence(chr_num, int(start_3p), int(end_3p), true_denominator_3p)
        hard_div_3p = t_div_3p[0]
        hard_div_denom_3p = t_div_3p[1]
        soft_div_3p = t_div_3p[2]
        soft_div_denom_3p = t_div_3p[3]
            
        results.write(name + '\t' + exon + '\t' + chr_num + '\t' + strand + '\t' + start_5p + '\t' + end_5p + '\t' + start_3p + '\t' + end_3p + '\t' + str(hard_div_5p) + '\t' + str(soft_div_5p) + '\t' + str(true_denominator_5p) + '\t' + str(hard_div_3p) + '\t' + str(soft_div_3p) + '\t' + str(true_denominator_3p) + '\n') 
    results.close()




#This section goes through all exons in the vervet UCSC exon or all intergenic regions 
#in the "vervet_all_intergenic_lengths.txt" file and calculates the 
#hard divergence value, soft divergence value, and number of callable sites for a 
#particular region. This can be performed for either exons on intergenic regions 
#by changing the "find_for" variable to either 'exons' or 'intergenic'. 

"""For ALL exons"""

find_for = 'exons'
# find_for = 'intergenic'
print('ran for: ' + find_for)

if find_for == 'exons':
    exon_file = 'C:/Users/Kellen/Downloads/primate_dfe_project/chimp_and_vervet_data/chlSab2.ensGenes/chlSab2.ensGenes.txt'
    with open(exon_file, 'r') as g:
        exon_data = g.readlines()[1:]

    result_file = 'C:/Users/Kellen/Downloads/primate_dfe_project/chimp_and_vervet_data/vervet_all_exons_divergence.txt'
    results = open(result_file, 'w+')
    results.write('name' + '\t' + 'exon' + '\t' + 'chr' + '\t' + 'strand' + '\t' + 'intergenic_5p_start' + '\t' + 'intergenic_5p_end' + '\t' + 'intergenic_3p_start' + '\t' + 'intergenic_3p_end' + '\t' + 'hard_div_5p' + '\t' + 'soft_div_5p' + '\t' + '5p_callability' + '\t' + 'hard_div_3p' + '\t' + 'soft_div_3p' + '\t' + '3p_callability' + '\n')
    total_num_lines = len(exon_data)
    for w in range(0, total_num_lines):
        line = exon_data[w].strip('\n')
        split = line.split()
        name = split[1]
        chr_num = split[2]
        strand = split[3]
        starts = split[9].split(',')[0 : -1]
        starts_int = [int(j) for j in starts]
        ends = split[10].split(',')[0 : -1]
        ends_int = [int(j) for j in ends]
        if chr_num == 'chrX' or chr_num == 'chrY' or chr_num == 'chrUn' or chr_num == 'chrM' or 'random' in chr_num:
            continue
        else:
            for j in range(0, len(ends_int)):
                start_pos = starts_int[j]
                end_pos = ends_int[j]

                chr_num_stripped = chr_num.strip('chr')
                cut_out = all_mask_data[chr_num_stripped][start_pos - 1 : end_pos]
                true_denominator = cut_out.count('0')

                t_div = get_divergence(chr_num, int(start_pos), int(end_pos), true_denominator)
                hard_div = t_div[0]
                hard_div_denom = t_div[1]
                soft_div = t_div[2]
                soft_div_denom = t_div[3]

                results.write(name + '_ex' + str(j + 1) + '\t' + chr_num + '\t' + strand + '\t' + str(start_pos) + '\t' + str(end_pos) + '\t' + str(hard_div) + '\t' + str(soft_div) + '\t' + str(true_denominator) + '\n')
    results.close()

    

if find_for == 'intergenic':
    intergenic_file = 'C:/Users/Kellen/Downloads/primate_dfe_project/chimp_and_vervet_data/vervet_all_intergenic_lengths.txt'
    with open(intergenic_file, 'r') as g:
        intergenic_data = g.readlines()[1:]

    result_file = 'C:/Users/Kellen/Downloads/primate_dfe_project/chimp_and_vervet_data/vervet_all_intergenic_divergence.txt'
    results = open(result_file, 'w+')
    results.write('name' + '\t' + 'exon' + '\t' + 'chr' + '\t' + 'strand' + '\t' + 'intergenic_5p_start' + '\t' + 'intergenic_5p_end' + '\t' + 'intergenic_3p_start' + '\t' + 'intergenic_3p_end' + '\t' + 'soft_div_5p' + '\t' + '5p_callability' + '\t' + 'soft_div_3p' + '\t' + '3p_callability' + '\n')

    total_num_lines = len(intergenic_data)
    for w in range(0, total_num_lines):
        line = intergenic_data[w].strip('\n')
        split = line.split()
        name = split[1]
        exon = split[2]
        chr_num = split[3]
        strand = split[4]
        start_5p = split[10]
        end_5p = split[11]
        len_5p = split[12]
        start_3p = split[13]
        end_3p = split[14]
        len_3p = split[15]
        
        if start_5p == 'NA' or end_5p == 'NA' or len_5p == 'NA':
            t_div_5p = ['NA', 'NA', 'NA', 'NA']
        else:
            chr_num_stripped = chr_num.strip('chr')
            cut_out = all_mask_data[chr_num_stripped][int(start_5p) - 1 : int(end_5p)]
            true_denominator_5p = cut_out.count('0')
            t_div_5p = get_divergence(chr_num, int(start_5p), int(end_5p), true_denominator_5p)
        hard_div_5p = t_div_5p[0]
        hard_div_denom_5p = t_div_5p[1]
        soft_div_5p = t_div_5p[2]
        soft_div_denom_5p = t_div_5p[3]

        if start_3p == 'NA' or end_3p == 'NA' or len_3p == 'NA':
            t_div_3p = ['NA', 'NA', 'NA', 'NA']
        else:
            chr_num_stripped = chr_num.strip('chr')
            cut_out = all_mask_data[chr_num_stripped][int(start_3p) - 1 : int(end_3p)]
            true_denominator_3p = cut_out.count('0')
            t_div_3p = get_divergence(chr_num, int(start_3p), int(end_3p), true_denominator_3p)
        hard_div_3p = t_div_3p[0]
        hard_div_denom_3p = t_div_3p[1]
        soft_div_3p = t_div_3p[2]
        soft_div_denom_3p = t_div_3p[3]
            
        results.write(name + '\t' + exon + '\t' + chr_num + '\t' + strand + '\t' + start_5p + '\t' + end_5p + '\t' + start_3p + '\t' + end_3p + '\t' + str(hard_div_5p) + '\t' + str(soft_div_5p) + '\t' + str(true_denominator_5p) + '\t' + str(hard_div_3p) + '\t' + str(soft_div_3p) + '\t' + str(true_denominator_3p) + '\n') 
    results.close()
    

