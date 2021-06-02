


#For various purposes, I created a Python script that goes through these post alignment files and 
#finds and outputs all positions within the post alignment files where one or all of the three nucleotides 
#needed to calculate the ancestral nucleotide are different.

#This script is currently designed to be run on a cluster and will require path changes on the user's part to the variables:
#'output_file' with the output file "all_sites_that_differ.txt"
#'epo_file' with the input files 'reformatted_alignment_chr' + str(chr_num) + '_' + str(segment_num) + '.txt'




maf_files = {'1' : 13, '2' : 12, '3' : 10, '4' : 10, '5' : 9, '6' : 9, '7' : 9, '8' : 8, '9' : 7,
             '10' : 7, '11' : 7, '12' : 7, '13' : 6, '14' : 5, '15' : 5, '16' : 5, '17' : 6,
             '18' : 4, '19' : 5, '20' : 3, '21' : 3, '22' : 3}

possible_human_chr = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22']

output_file = open('/scratch/kriall/primate_dfe/13_primates_EPO/all_sites_that_differ.txt', 'w+')
output_file.write('position' + '\t' + 'vervet' + '\t' + 'vervet-human' + '\t' + 'vervet-human-marmoset' + '\t' + 'human' + '\n')

for chr_num in possible_human_chr:
    num_of_segments = maf_files.get(str(chr_num))
#     print(num_of_segments)
#     print('\n')
    for segment_num in range(1, num_of_segments + 1):
#         epo_file = 'C:/Users/Kellen/Downloads/filtered_alignment_chr' + str(chr_num) + '_' + str(segment_num) + '.txt'
        epo_file = '/scratch/kriall/primate_dfe/13_primates_EPO/reformatted_aligns/reformatted_alignment_chr' + str(chr_num) + '_' + str(segment_num) + '.txt'
        with open(epo_file, 'r') as f:
            content = f.readlines()
        
        total_num_lines = len(content)
        for i in range(0, total_num_lines):
            split = content[i].strip('\n').split()
            vervet_pos = split[0]
            if 'chr' in vervet_pos:
                v_nuc = split[1].upper()
#                 vh_nuc = split[2].upper()
#                 vhm_nuc = split[3].upper()
                h_nuc = split[4].upper()
                
                if h_nuc != '-':
                    if v_nuc != h_nuc:
                        output_file.write(content[i])
                        
output_file.close()

