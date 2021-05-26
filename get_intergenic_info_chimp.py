


#This Python script calculates the lengths of flanking 5’ and 3’ intergenic regions around both the identified 
#homologous exons, as well as for all exons, specifically for chimps.

#This script is currently designed to be run locally, and will require path changes on the user's part to the variables:
#'results_file' with output file "chimp_all_intergenic_lengths.txt"
#'ucsc_file' with input files 'panTro2_' + compare_chr +'.ensGenes.txt'
#'homo_exons_file' with input file "all_homo_exons_reference_exons_only.txt"
#'chimp_intergenic_file' with with input file "chimp_all_intergenic_lengths.txt"
#'results_file' with output file "chimp_homo_intergenic_lengths.txt"




import sys
import math




#This is a fairly complicated section that essentially reads in all of the exons in the CHIMP UCSC file and builds a 
#theoretical positive and negative strand for each chromosome consisting of all of the exons in the list. We then 
#assume that all regions between these exons are intergenic regions. A problem we ran into is that some exons clearly 
#show up multiple times, but each listing may be different by a few base pairs. This creates a problem in that, by 
#taking the regions in between, if the same exon has one start given as 114, and another given as 115, then the 
#intergenic region between them is taken as zero. To compensate for this, we TRY OUR BEST (this is still imperfect) 
#to find any duplicated exons and take them to span from the earliest given start position to the furthest given end 
#position. For example, if an exon is given once as 10 – 21, again as 10 – 20, and again as 9 – 20, we would take the
#exon to span from 9 – 21. The issue is that we can only identify exons that are listed multiple times if two or more
#of those listings share either a common start or end. It’s not a perfect solution, but it’s better than doing nothing.
#Once the positive and negative strands have been assembled, it then calculates and outputs the intergenic lengths for all exons.
#Additionally, any exons found to have duplicates are marked with a "N" to let you know not to use them in future analysis. 
#Otherwise, they're fine and marked with a "Y".

results_file = 'C:/Users/Kellen/Downloads/primate_dfe_project/chimp_and_vervet_data/chimp_all_intergenic_lengths.txt'
results = open(results_file,'w+')
results.write('use?' + '\t' + 'name' + '\t' + 'exon' + '\t' + 'chr' + '\t' + 'strand' + '\t' + 'UCSC_start' + '\t' + 'UCSC_end' + '\t' + 'exon_start' + '\t' + 'exon_end' + '\t' + 'exon_len' + '\t' + 'intergenic_5p_start' + '\t' + 'intergenic_5p_end' + '\t' + 'intergenic_5p_len' + '\t' + 'intergenic_3p_start' + '\t' + 'intergenic_3p_end' + '\t' + 'intergenic_3p_len'+'\n')

acceptable_chr = ['1', '2a', '2b', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22']
for chr_num in acceptable_chr:
    compare_chr = 'chr' + chr_num
    ucsc_file = 'C:/Users/Kellen/Downloads/primate_dfe_project/chimp_and_vervet_data/panTro2.ensGenes/panTro2_' + compare_chr +'.ensGenes.txt'
    with open(ucsc_file, 'r') as f:
        ucsc_data = f.readlines()
    
    total_num_lines = len(ucsc_data)
    
    dict_starts_check = {}
    dict_ends_check = {}
    positive_true_starts_reference = {}
    negative_true_starts_reference = {}
    positive_true_ends_reference = {}
    negative_true_ends_reference = {}
    positive_strand_assembly_starts_check = []
    negative_strand_assembly_starts_check = []
    positive_strand_assembly_ends_check = []
    negative_strand_assembly_ends_check = []
    
    
    
    #This sub-section makes the dictionary keys with a list initially containing the strand info
    for w in range(0, total_num_lines):
        line = ucsc_data[w].strip('\n')
        split = line.split()
        strand = split[3]
        starts = split[9].split(',')[0 : -1]
        starts_int = [int(q) for q in starts]
        ends = split[10].split(',')[0 : -1]
        ends_int = [int(q) for q in ends]

        for j in range(0, len(ends)):
            start_key = str(starts_int[j]) + strand
            end_key = str(ends_int[j]) + strand
            dict_starts_check[start_key] = [strand]
            dict_ends_check[end_key] = [strand]

            
        
    #This sub-section adds all of the start and end positions to the respective lists
    for y in range(0, total_num_lines):
        line = ucsc_data[y].strip('\n')
        split = line.split()
        strand = split[3]
        starts = split[9].split(',')[0 : -1]
        starts_int = [int(q) for q in starts]
        ends = split[10].split(',')[0: -1]
        ends_int = [int(q) for q in ends]

        for j in range(0,len(ends)):
            start = starts_int[j]
            start_key = str(starts_int[j]) + strand
            end = ends_int[j]
            end_key = str(ends_int[j]) + strand
            dict_starts_check[start_key].append(end)
            dict_ends_check[end_key].append(start)

            
    #This sub-section builds a positive and negative strand as a list consisting of exon start and end positions with the expectation of all gaps being intergenic regions, trying its best to settle on a singluar size for any duplicated exons with varying listed start and end positions 
    for start in dict_starts_check:
        key_info = dict_starts_check.get(start)
        strand = key_info[0]
        start = int(start.strip(strand))
        end_positions = key_info[1 :]
        end_positions_duplicate_check = list(set(end_positions))
                
        if len(end_positions_duplicate_check) == 1:
            end = end_positions_duplicate_check[0]
            if strand == '+':
                positive_strand_assembly_starts_check.append(start)
                positive_strand_assembly_starts_check.append(end)
                positive_true_ends_reference[end] = [end, 'Y']
            elif strand == '-':
                negative_strand_assembly_starts_check.append(start)
                negative_strand_assembly_starts_check.append(end)
                negative_true_ends_reference[end] = [end, 'Y']
        else:
            max_end_position = max(end_positions)
            if strand == '+':
                positive_strand_assembly_starts_check.append(start)
                positive_strand_assembly_starts_check.append(max_end_position)
            elif strand == '-':
                negative_strand_assembly_starts_check.append(start)
                negative_strand_assembly_starts_check.append(max_end_position)
            for bad_end in end_positions_duplicate_check:
                if strand == '+':
                    positive_true_ends_reference[bad_end] = [max_end_position, 'N']
                elif strand == '-':
                    negative_true_ends_reference[bad_end] = [max_end_position, 'N']
      


    for end in dict_ends_check:
        key_info = dict_ends_check.get(end)
        strand = key_info[0]
        end = int(end[0 : -1])
        start_positions = key_info[1 :]
        start_positions_duplicate_check = list(set(start_positions))
                
        if len(start_positions_duplicate_check) == 1:
            start = start_positions_duplicate_check[0]
            if strand == '+':
                positive_strand_assembly_ends_check.append(start)
                positive_strand_assembly_ends_check.append(end)
                positive_true_starts_reference[start] = [start, 'Y']
            elif strand == '-':
                negative_strand_assembly_ends_check.append(start)
                negative_strand_assembly_ends_check.append(end)
                negative_true_starts_reference[start] = [start, 'Y']
        else:
            min_start_position = min(start_positions)
            if strand == '+':
                positive_strand_assembly_ends_check.append(min_start_position)
                positive_strand_assembly_ends_check.append(end)
            elif strand == '-':
                negative_strand_assembly_ends_check.append(min_start_position)
                negative_strand_assembly_ends_check.append(end)
            for start in start_positions_duplicate_check:
                if strand == '+':
                    positive_true_starts_reference[start] = [min_start_position, 'N']
                elif strand == '-':
                    negative_true_starts_reference[start] = [min_start_position, 'N']
        
    

    positive_strand_assembly = list(set(positive_strand_assembly_starts_check).intersection(set(positive_strand_assembly_ends_check)))
    positive_strand_assembly = list(set(positive_strand_assembly))
    positive_strand_assembly.sort()
    positive_strand_assembly.insert(0, 0)
    positive_strand_assembly.append('NA')

    negative_strand_assembly = list(set(negative_strand_assembly_starts_check).intersection(set(negative_strand_assembly_ends_check)))
    negative_strand_assembly = list(set(negative_strand_assembly))
    negative_strand_assembly.sort()
    negative_strand_assembly.insert(0, 0)
    negative_strand_assembly.append('NA')
    
#     print(negative_true_ends_reference)
#     print(dict_ends_check)
#     print(dict_starts_check)
#     print(positive_strand_assembly_starts_check)
#     print(positive_strand_assembly_ends_check)
#     print(positive_strand_assembly)
#     sys.exit()

    #This sub-section finds an exon's position in the constructed positive or negative strand, finds its flanking 5' and 3' intergenic regions, and outputs that information
    for i in range(0, total_num_lines):
        line = ucsc_data[i].strip('\n')
        split = line.split()
        name1 = split[1]
        name2 = split[12]
        chr_num = split[2]
        strand = split[3]

        starts = split[9].split(',')[0 : -1]
        starts_int = [int(j) for j in starts]
        ends = split[10].split(',')[0 : -1]
        ends_int = [int(j) for j in ends]

        all_pairs = []
        for j in range(0, len(ends_int)):
            start = starts_int[j]
            end = ends_int[j]
            if strand == '+':
                true_start = positive_true_starts_reference[start][0]
                true_end = positive_true_ends_reference[end][0]
                start_useability = positive_true_starts_reference[start][1]
                end_useability = positive_true_ends_reference[end][1]
            elif strand == '-':
                true_start = negative_true_starts_reference[start][0]
                true_end = negative_true_ends_reference[end][0]
                start_useability = negative_true_starts_reference[start][1]
                end_useability = negative_true_ends_reference[end][1]
        
            if start_useability == 'N' or end_useability == 'N':
                exon_useability = 'N'
            else:
                exon_useability = 'Y'
                
#             print(str(start)+' - '+str(true_start)+'   /   '+str(end)+' - '+str(true_end))
            pair = [true_start, true_end, exon_useability, start, end]
            all_pairs.append(pair)
        
        num_exons_in_transcript = len(all_pairs)
        for k in range(0, len(all_pairs)):
            this_exon_start = all_pairs[k][0]
            this_exon_end = all_pairs[k][1]
            this_exon_length = this_exon_end - this_exon_start + 1
            exon_useability = all_pairs[k][2]
            UCSC_start = all_pairs[k][3]
            UCSC_end = all_pairs[k][4]

            if strand == '+':
                previous_exon_end = positive_strand_assembly[positive_strand_assembly.index(this_exon_start) - 1]
                next_exon_start = positive_strand_assembly[positive_strand_assembly.index(this_exon_end) + 1]

                if previous_exon_end == 'NA':
                    intergenic_5p_start = 'NA'
                    intergenic_5p_end = this_exon_start - 1
                    intergenic_5p_length = 'NA'
                else:
                    doubled_intergenic_5p_start = previous_exon_end + 1
                    intergenic_5p_end = this_exon_start - 1
                    intergenic_5p_start = int(math.ceil((doubled_intergenic_5p_start + intergenic_5p_end) / 2))
                    intergenic_5p_length = intergenic_5p_end - intergenic_5p_start + 1

                if next_exon_start == 'NA':
                    intergenic_3p_start = this_exon_end + 1
                    intergenic_3p_end = 'NA'
                    intergenic_3p_length = 'NA'
                else:
                    intergenic_3p_start = this_exon_end + 1
                    doubled_intergenic_3p_end = next_exon_start - 1
                    intergenic_3p_end = int(math.ceil((intergenic_3p_start + doubled_intergenic_3p_end) / 2))
                    intergenic_3p_length = intergenic_3p_end - intergenic_3p_start + 1

            elif strand == '-':
                previous_exon_end = negative_strand_assembly[negative_strand_assembly.index(this_exon_start) - 1]
                next_exon_start = negative_strand_assembly[negative_strand_assembly.index(this_exon_end) + 1]

                if previous_exon_end == 'NA':
                    intergenic_3p_start = 'NA'
                    intergenic_3p_end = this_exon_start - 1
                    intergenic_3p_length = 'NA'
                else:
                    doubled_intergenic_3p_start = previous_exon_end + 1
                    intergenic_3p_end = this_exon_start - 1
                    intergenic_3p_start = int(math.ceil((doubled_intergenic_3p_start + intergenic_3p_end) / 2))
                    intergenic_3p_length = intergenic_3p_end - intergenic_3p_start + 1

                if next_exon_start == 'NA':
                    intergenic_5p_start = this_exon_end + 1
                    intergenic_5p_end = 'NA'
                    intergenic_5p_length = 'NA'
                else:
                    intergenic_5p_start = this_exon_end + 1
                    doubled_intergenic_5p_end = next_exon_start - 1
                    intergenic_5p_end = int(math.ceil((intergenic_5p_start + doubled_intergenic_5p_end) / 2))
                    intergenic_5p_length = intergenic_5p_end - intergenic_5p_start + 1

            results.write(exon_useability + '\t' + name1 + '\t' + name2 + '_ex' + str(k + 1) + '/' + str(num_exons_in_transcript) + '\t' + chr_num + '\t' + strand + '\t' + str(UCSC_start) + '\t' + str(UCSC_end) + '\t' + str(this_exon_start) + '\t' + str(this_exon_end) + '\t' + str(this_exon_length) + '\t' + str(intergenic_5p_start) + '\t' + str(intergenic_5p_end) + '\t' + str(intergenic_5p_length) + '\t' + str(intergenic_3p_start) + '\t' + str(intergenic_3p_end) + '\t' + str(intergenic_3p_length) + '\n')            

#             print(previous_exon_end)
#             print(this_exon_start)
#             print(this_exon_end)
#             print(next_exon_start)
#             print('\n')

#     print(positive_strand_assembly)
#     print(negative_strand_assembly)
    
results.close()




#This section references the all intergenic lengths file and matches up the 248 homologous exons to get the intergenic lengths for them.
homo_exons_file = 'C:/Users/Kellen/Downloads/primate_dfe_project/chimp_and_vervet_data/all_homo_exons_reference_exons_only.txt'
chimp_intergenic_file = 'C:/Users/Kellen/Downloads/primate_dfe_project/chimp_and_vervet_data/chimp_all_intergenic_lengths.txt'
results_file = 'C:/Users/Kellen/Downloads/primate_dfe_project/chimp_and_vervet_data/chimp_homo_intergenic_lengths.txt'

results = open(results_file, 'w+')

with open(homo_exons_file, 'r') as g:
    homo_exon_data = g.readlines()[1 : ]
    
with open(chimp_intergenic_file, 'r') as f:
    chimp_intergenic_data = f.readlines()
results.write(chimp_intergenic_data[0])
del chimp_intergenic_data[0]

tally = 0
for i in homo_exon_data:
    found = False
    split = i.strip('\n').split()
    transcript_name = split[12]
    gene_name = split[13]
    chr_num = split[14]
    strand = split[15]
    start = split[16]
    end = split[17]
    for j in chimp_intergenic_data:
        split2 = j.strip('\n').split()
        inter_transcript_name = split2[1]
        inter_gene_name = split2[2]
        inter_chr_num = split2[3]
        inter_strand = split2[4]
        inter_start = split2[5]
        inter_end = split2[6]

        if transcript_name == inter_transcript_name and gene_name == inter_gene_name and chr_num == inter_chr_num and strand == inter_strand and start == inter_start and end == inter_end:
            results.write(j)
            tally += 1
            found = True
            break
    if found == False:
        print(i)
        print(end)
        print(inter_end)
    
results.close()
print(tally)

