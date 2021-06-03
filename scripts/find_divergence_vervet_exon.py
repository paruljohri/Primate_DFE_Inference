


#This is a Python script that calculates divergence rate values for vervet EXONS specifically. 
#Whereas most scripts that have separate "find x for chimps" and "find x for vervets" versions
#have none or extremely minor differences, these chimp and vervet divergence calculation 
#scripts are substantially different due to how they go about calculating whether or 
#not a site is diverged in addition to the vervet portion being split into separate scripts for 
#intergenic and exonic regions.

#This script is currently designed to be run locally and will require path changes on the user's part to the variables:
#'mask_file' with input files 'chr' + chr_num + '.mask.fa'
#'divergence_file' with input files "all_sites_that_differ.txt"
#'region_data_file' with input file "248_homo_exons_reference_with_intergenic_lengths.txt"
#first instance of 'result_file' with output file "vervet_homo_exons_divergence.txt"
#'exon_file' with input file "chlSab2.ensGenes.txt"
#second instance of 'result_file' with output file "vervet_all_exons_divergence.txt"




#This section reads in all of the data from the vervet mask files, allowing us to determine if 
#a site is callable. The number of callable sites is then used in the denominator of the divergence calculation.

all_chr_nums = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29']
all_mask_data = {}
for chr_num in all_chr_nums:
    mask_file = 'C:/Users/Kellen/Downloads/primate_dfe_project/chimp_and_vervet_data/vervet_mask/chr' + chr_num + '_mask.fasta'
    with open(mask_file) as f:
        mask_data = f.readlines()[1 : ]
    key = 'chr' + str(chr_num)
    all_mask_data[key] = mask_data[0]
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
possible_nuc = ['A', 'a', 'T', 't', 'G', 'g', 'C', 'c']
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

        if ancestral_nuc in possible_nuc:
            if v_nuc.upper() != ancestral_nuc.upper():
                d_div[total_position] = 'Y'
                




#This function takes in the start and end position of a region one wishes to find the divergence 
#for, tallies up the number of diverged sites, and then obtains a final divergence value by dividing
#that tally by the total number of callable sites in the region. Only sites that are callable are considered, 
#so positions that may be diverged but are not callable are not considered in the final calculation.
#Regions with no callable sites are given a divergence value of "NA" in the output file.

def get_divergence(chr_num, start_pos, end_pos):
    chr_sequence = all_mask_data[chr_num]

    div_numerator = 0
    div_denominator = 0
    
    i = int(start_pos)
    while i <= int(end_pos):
        is_diverged = 'N'
        s = chr_num + '_' + str(i)
        at_site = d_div.get(s, 'NA')
        if at_site != 'NA':
            is_diverged = 'Y'
            
        is_callable = 'N'    
        site_mask_value = chr_sequence[i - 1]
        if site_mask_value == '0':
            is_callable = 'Y'
            
        if is_callable == 'Y':
            div_denominator += 1
            if is_diverged == 'Y':
                div_numerator += 1
            
        i += 1
    
    
    if div_denominator == 0:
        divergence = 'NA'
    else:
        divergence = div_numerator / div_denominator
        divergence = round(divergence * 10**9) / 10**9
    
    return (divergence, div_denominator)




#This section goes through all the established homologous vervet exons and calculates the 
#hard divergence value, soft divergence value, and number of callable sites.

"""For HOMO regions"""
region_data_file = 'C:/Users/Kellen/Downloads/primate_dfe_project/248_homo_exons_reference_with_intergenic_lengths.txt'
with open(region_data_file, 'r') as g:
    region_data = g.readlines()[1 : ]

result_file = 'C:/Users/Kellen/Downloads/primate_dfe_project/chimp_and_vervet_data/vervet_homo_exons_divergence.txt'
results = open(result_file, 'w+')
results.write('name' + '\t' + 'exon' + '\t' + 'chr' + '\t' + 'strand' + '\t' + 'start' + '\t' + 'end' + '\t' + 'div' + '\t' + 'callability' + '\n')

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
        div, div_denominator = get_divergence(chr_num, start_pos, end_pos)    

    results.write(name + '\t' + exon + '\t' + chr_num + '\t' + strand + '\t' + str(start_pos) + '\t' + str(end_pos) + '\t' + str(div) + '\t' + str(div_denominator) + '\n')
results.close()




#This section goes through all exons in the vervet UCSC exon or all intergenic regions 
#in the "vervet_all_intergenic_lengths.txt" file and calculates the 
#hard divergence value, soft divergence value, and number of callable sites.

"""For ALL exons"""
exon_file = 'C:/Users/Kellen/Downloads/primate_dfe_project/chimp_and_vervet_data/chlSab2.ensGenes/chlSab2.ensGenes.txt'
with open(exon_file, 'r') as g:
    exon_data = g.readlines()[1 : ]

result_file = 'C:/Users/Kellen/Downloads/primate_dfe_project/chimp_and_vervet_data/vervet_all_exons_divergence.txt'
results = open(result_file, 'w+')
results.write('name' + '\t' + 'chr' + '\t' + 'strand' + '\t' + 'start' + '\t' + 'end' + '\t' + 'div' + '\t' + 'callability' + '\n')

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

            div, div_denominator = get_divergence(chr_num, int(start_pos), int(end_pos)) 

            results.write(name + '_ex' + str(j + 1) + '\t' + chr_num + '\t' + strand + '\t' + str(start_pos) + '\t' + str(end_pos) + '\t' + str(div) + '\t' + str(div_denominator) + '\n')
results.close()

