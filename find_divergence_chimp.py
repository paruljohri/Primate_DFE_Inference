


#This is a Python script that calculates divergence rate values for chimps. 
#Whereas most scripts that have separate "find x for chimps" and "find x for vervets" versions
#have none or extremely minor differences, these chimp and vervet divergence calculation 
#scripts are substantially different due to how they go about calculating whether or 
#not a site is diverged.

#This script is currently designed to be run locally and will require path changes on the user's part to the variables:
#'mask_file' with input files 'chr' + chr_num + '.mask.fa'
#'divergence_file' with input files "divergent_final_pt2.txt"
#'region_data_file' with input file "248_homo_exons_reference_with_intergenic_lengths.txt"
#first instance of 'result_file' with output file "chimp_homo_exons_divergence.txt"
#second instance of 'result_file' with output file "chimp_homo_intergenic_divergence.txt"
#'exon_file' with input file "panTro2.ensGenes.txt"
#third instance of 'result_file' with output file "chimp_all_exons_divergence.txt"
#'intergenic_file' with input file "chimp_all_intergenic_lengths.txt"
#fourth instance of 'result_file' with output file "chimp_all_intergenic_divergence.txt"




#This section reads in all of the data from the chimp mask files, allowing us to determine if 
#a site is callable. The number of callable sites is then used in the denominator of the divergence calculation.

all_chr_nums = ['1', '2a', '2b', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22']
all_mask_data = {}
for chr_num in all_chr_nums:
    mask_file = 'C:/Users/Kellen/Downloads/primate_dfe_project/chimp_and_vervet_data/chimp_mask/chr' + chr_num + '.mask.fa'
    with open(mask_file) as f:
        mask_data = f.readlines()[1 :]
    mask_data = [j.strip('\n') for j in mask_data]
    join_mask_data = ''.join(mask_data)
    all_mask_data[chr_num] = join_mask_data
mask_data = []




#This section reads in information from the file "divergent_final_pt2.txt" (based on the PanTro2 assembly)
#containing chimp sites where the chimp allele and aligned human differ, as well as the true ancestral 
#allele. It then determines if the chimp allele is diverged at that site or not. Finally, it adds any
#of the diverged positions to a dictionary for further calculation.

divergence_file = 'C:/Users/Kellen/Downloads/primate_dfe_project/chimp_and_vervet_data/chimp_divergence/divergent_final_pt2.txt'
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
        chr_num = split[0]
        position = split[1]
        total_position = chr_num + '_' + position
        human_nuc = split[1]
        chimp_nuc = split[2]
        ancestral_nuc = split[3]

        if ancestral_nuc in hard:
            if chimp_nuc != ancestral_nuc:
                d_div[total_position] = 'H'
        elif ancestral_nuc in soft:
            if chimp_nuc.lower() != ancestral_nuc.lower():
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




#This section goes through all the established homologous chimp exons and calculates the 
#hard divergence value, soft divergence value, and number of callable sites for a 
#particular region. This can be performed for either exons on intergenic regions 
#by changing the "find_for" variable to either 'exons' or 'intergenic'. 

"""For HOMO regions"""

find_for = 'exons'
# find_for = 'intergenic'
print('ran for: ' + find_for)

region_data_file = 'C:/Users/Kellen/Downloads/primate_dfe_project/chimp_and_vervet_data/248_homo_exons_reference_with_intergenic_lengths.txt'
with open(region_data_file, 'r') as g:
    region_data = g.readlines()[1:]

if find_for == 'exons':
    result_file = 'C:/Users/Kellen/Downloads/primate_dfe_project/chimp_and_vervet_data/chimp_homo_exons_divergence.txt'
    results = open(result_file, 'w+')
    results.write('name' + '\t' + 'exon' + '\t' + 'chr' + '\t' + 'strand' + '\t' + 'start' + '\t' + 'end' + '\t' + 'hard_div' + '\t' + 'soft_div' + '\t' + 'callability' + '\n')

    total_num_lines = len(region_data)
    for w in range(0, total_num_lines):
        line = region_data[w].strip('\n')
        split = line.split()
        name = split[16]
        exon = split[17]
        chr_num = split[18]
        strand = split[19]
        start_pos = int(split[22])
        end_pos = int(split[23])
        
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
    result_file = 'C:/Users/Kellen/Downloads/primate_dfe_project/chimp_and_vervet_data/chimp_homo_intergenic_divergence.txt'
    results = open(result_file, 'w+')
    results.write('name' + '\t' + 'exon' + '\t' + 'chr' + '\t' + 'strand' + '\t' + 'intergenic_5p_start' + '\t' + 'intergenic_5p_end' + '\t' + 'intergenic_3p_start' + '\t' + 'intergenic_3p_end' + '\t' + 'hard_div_5p' + '\t' + '5p_callability' + '\t' + 'hard_div_3p' + '\t' + '3p_callability' + '\n')

    total_num_lines = len(region_data)
    for w in range(0, total_num_lines):
        line = region_data[w].strip('\n')
        split = line.split()
        name = split[16]
        exon = split[17]
        chr_num = split[18]
        strand = split[19]
        start_5p = split[25]
        end_5p = split[26]
        len_5p = split[27]
        start_3p = split[28]
        end_3p = split[29]
        len_3p = split[30]
        
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
            
        results.write(name + '\t' + exon + '\t' + chr_num + '\t' + strand + '\t' + start_5p + '\t' + end_5p + '\t' + start_3p + '\t' + end_3p + '\t' + str(hard_div_5p) + '\t' + str(true_denominator_5p) + '\t' + str(hard_div_3p) + '\t' + str(true_denominator_3p) + '\n') 
    results.close()
    




#This section goes through all exons in the chimp UCSC exon or all intergenic regions 
#in the "vervet_all_intergenic_lengths.txt" file and calculates the 
#hard divergence value, soft divergence value, and number of callable sites for a 
#particular region. This can be performed for either exons on intergenic regions 
#by changing the "find_for" variable to either 'exons' or 'intergenic'.

"""For ALL exons"""

find_for = 'exons'
# find_for = 'intergenic'
print('ran for: ' + find_for)

if find_for == 'exons':
    exon_file = 'C:/Users/Kellen/Downloads/primate_dfe_project/chimp_and_vervet_data/panTro2.ensGenes/panTro2.ensGenes.txt'
    with open(exon_file, 'r') as g:
        exon_data = g.readlines()[1:]

    result_file = 'C:/Users/Kellen/Downloads/primate_dfe_project/chimp_and_vervet_data/chimp_all_exons_divergence.txt'
    results = open(result_file, 'w+')
    results.write('name' + '\t' + 'chr' + '\t' + 'strand' + '\t' + 'start' + '\t' + 'end' + '\t' + 'soft_div' + '\t' + 'callability' + '\n')

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

                results.write(name + '_ex' + str(j + 1) + '\t' + chr_num + '\t' + strand + '\t' + str(start_pos) + '\t' + str(end_pos) + '\t' + str(soft_div) + '\t' + str(true_denominator) + '\n')
    results.close()

    

if find_for == 'intergenic':
    intergenic_file = 'C:/Users/Kellen/Downloads/primate_dfe_project/chimp_and_vervet_data/chimp_all_intergenic_lengths.txt'
    with open(intergenic_file, 'r') as g:
        intergenic_data = g.readlines()[1:]

    result_file = 'C:/Users/Kellen/Downloads/primate_dfe_project/chimp_and_vervet_data/chimp_all_intergenic_divergence.txt'
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
            
        results.write(name + '\t' + exon + '\t' + chr_num + '\t' + strand + '\t' + start_5p + '\t' + end_5p + '\t' + start_3p + '\t' + end_3p + '\t' + str(soft_div_5p) + '\t' + str(true_denominator_5p) + '\t' + str(soft_div_3p) + '\t' + str(true_denominator_3p) + '\n') 
    results.close()
    

