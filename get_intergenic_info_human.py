


#This Python script is designed to find the human intergenic region lengths for all homologous exons. While these were 
#known from the get-go, now that we are only focusing on the 248 found to be homologous with vervets and chimps, we first 
#need to pull out the lengths for those 248. Additionally, Parul had originally calculated these to span from one exon 
#all the way to the next exon. However, for our purposes, we are considering the flanking intergenic region for an exon 
#to only extend from that exon to the midpoint between it and the next exon. Therefore, the 248 intergenic lengths given 
#initially in Parul's file also need to be cut in half.

#This script is currently designed to be run locally, and will require path changes on the user's part to the variables:
#'homo_exons_file' with input file "all_homo_exons_reference_exons_only.txt"
#'human_intergenic_file' with input file "converted_hg38.txt"
#'results_file' with output file "DELETE_human_homo_intergenic_lengths.txt"
#'intergenic_file' with input file "DELETE_human_homo_intergenic_lengths.txt"
#'result_file' with output file "human_homo_intergenic_lengths.txt"




import math
import os




#This section reads in the original list of 465 exons we started with to get, matches the listings with the 248 we 
#have now, and then outputs those lengths in a temporary file that is then deleted afterwards.
homo_exons_file = 'C:/Users/Kellen/Downloads/primate_dfe_project/chimp_and_vervet_data/all_homo_exons_reference_exons_only.txt'
human_intergenic_file = 'C:/Users/Kellen/Downloads/primate_dfe_project/chimp_and_vervet_data/converted_hg38.txt'
results_file = 'C:/Users/Kellen/Downloads/primate_dfe_project/chimp_and_vervet_data/DELETE_human_homo_intergenic_lengths.txt'

results = open(results_file, 'w+')


with open(homo_exons_file, 'r') as g:
    homo_exon_data = g.readlines()[1 :]
    
with open(human_intergenic_file, 'r') as f:
    human_intergenic_data = f.readlines()[1 :]
results.write('#' + '\t' + 'chr' + '\t' + 'strand' + '\t' + 'exon_start' + '\t' + 'exon_end' + '\t' + 'exon_len' + '\t' + 'intergenic_5p_start' + '\t' + 'intergenic_5p_end' + '\t' + 'intergenic_5p_len' + '\t' + 'intergenic_3p_start' + '\t' + 'intergenic_3p_end' + '\t' + 'intergenic_3p_len' + '\n')
    
for i in homo_exon_data:
    split = i.strip('\n').split()
    number = split[0]
    gene = split[1]
    chr_num = split[2]
    strand = split[3]
    start = split[6]
    end = split[7]
    for j in human_intergenic_data:
        split2 = j.strip('\n').split()
        inter_number = split2[0]
        inter_gene = split2[1]
        inter_chr_num = split2[2]
        inter_strand = split2[4]
        inter_start = split2[6]
        inter_end = split2[7]
        
        inter_exon_len = split2[5]
        inter_5p_start = split2[8]
        inter_5p_end = split2[9]
        inter_5p_len = split2[12]
        inter_3p_start = split2[10]
        inter_3p_end = split2[11]
        inter_3p_len = split2[13]
        
        if  gene == inter_gene and chr_num == inter_chr_num and strand == inter_strand and start == inter_start and end == inter_end:
#             print(j)
            results.write(number + '\t' + chr_num + '\t' + strand + '\t' + inter_start + '\t' + inter_end + '\t' + inter_exon_len + '\t' + inter_5p_start + '\t' + inter_5p_end + '\t' + inter_5p_len + '\t' + inter_3p_start + '\t' + inter_3p_end + '\t' + inter_3p_len + '\n')
            break
results.close()




#This section takes all of the intergenic lengths and cuts them in half. It also deletes the intermediate file.
"""For HUMAN HOMO intergenic regions"""
intergenic_file = 'C:/Users/Kellen/Downloads/primate_dfe_project/chimp_and_vervet_data/DELETE_human_homo_intergenic_lengths.txt'
with open(intergenic_file, 'r') as g:
    intergenic_data = g.readlines()
    
result_file = 'C:/Users/Kellen/Downloads/primate_dfe_project/chimp_and_vervet_data/human_homo_intergenic_lengths.txt'
results = open(result_file, 'w+')
results.write(intergenic_data[0])
del intergenic_data[0]

for i in intergenic_data:
    split = i.strip('\n').split()
    strand = split[2]
    start_5p = split[6]
    end_5p = split[7]
    start_3p = split[9]
    end_3p = split[10]
    
#     print(strand)
    if strand == '+':
        if start_5p == 'NA' or end_5p == 'NA':
            new_5p_start = 'NA'
            new_5p_len = 'NA'
        else:
            start_5p = int(start_5p)
            end_5p = int(end_5p)
            midpoint_5p = (start_5p + end_5p) / 2
            new_5p_start = int(math.ceil(midpoint_5p))
            new_5p_len = end_5p - new_5p_start + 1
        
        if start_3p == 'NA' or end_3p == 'NA':
            new_3p_end = 'NA'
            new_3p_len = 'NA'
        else:
            start_3p = int(start_3p)
            end_3p = int(end_3p)
            midpoint_3p = (start_3p + end_3p) / 2
            new_3p_end = int(math.floor(midpoint_3p))
            new_3p_len = new_3p_end - start_3p + 1
        
#         print(str(new_5p_start) +' - ' + str(end_5p))
#         print(new_5p_len)
#         print(str(start_3p) + ' - ' + str(new_3p_end))
#         print(new_3p_len)
        
        split[6] = str(new_5p_start)
        split[8] = str(new_5p_len)
        split[10] = str(new_3p_end)
        split[11] = str(new_3p_len)
        new_line = '\t'.join(split)
        results.write(new_line + '\n')
        
    elif strand == '-':
        if start_5p == 'NA' or end_5p == 'NA':
            new_5p_end = 'NA'
            new_5p_len = 'NA'
        else:
            start_5p = int(start_5p)
            end_5p = int(end_5p)
            midpoint_5p = (start_5p + end_5p) / 2
            new_5p_end = int(math.floor(midpoint_5p))
            new_5p_len = new_5p_end - start_5p + 1

        if start_3p == 'NA' or end_3p == 'NA':
            new_3p_start = 'NA'
            new_3p_len = 'NA'
        else:
            start_3p = int(start_3p)
            end_3p = int(end_3p)
            midpoint_3p = (start_3p + end_3p) / 2
            new_3p_start = int(math.ceil(midpoint_3p))
            new_3p_len = end_3p - new_3p_start + 1

#         print(str(start_5p) +' - ' + str(new_5p_end))
#         print(new_5p_len)
#         print(str(new_3p_start) + ' - ' + str(end_3p))
#         print(new_3p_len)
    
        split[7] = str(new_5p_end)
        split[8] = str(new_5p_len)
        split[9] = str(new_3p_start)
        split[11] = str(new_3p_len)
        new_line = '\t'.join(split)
        results.write(new_line + '\n')
        
results.close()
os.remove(intergenic_file)

