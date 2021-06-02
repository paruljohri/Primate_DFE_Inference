


#This is a Python script designed to add information to the output files created in the 
#"find_vervet_alignments.py" script. It goes through each of those files and makes a copy 
#while also adding new columns that include the calculated ancestral allele, whether or not 
#that vervet position is diverged, and whether or not that site is callable from the masking information.

#This script is currently designed to be run on a cluster and will require path changes on the user's part to the variables:
#'mask_file' with input files 'final_10ind_AGM_mappings_mask_chr' + chr_num + '.fasta'
#'output_file' with output files 'reformatted_alignment_with_callability_chr' + str(chr_num) + '_' + str(segment_num) + '.txt'
#'epo_file' with input files 'reformatted_alignment_chr' + str(chr_num) + '_' + str(segment_num) + '.txt'




#This section reads in all of the data from the vervet mask files, allowing us to determine if 
#a site is callable. The number of callable sites is then used in the denominator of the divergence calculation.

possible_vervet_chr = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29']

all_mask_data = {}
for chr_num in possible_vervet_chr:
    mask_file = '/scratch/kriall/primate_dfe/to_drive/Parul_and_Kellen/vervet/mask/final_10ind_AGM_mappings_mask_chr' + chr_num + '.fasta'
    with open(mask_file) as f:
        mask_data = f.readlines()[1 :]
    all_mask_data[chr_num] = mask_data[0]




maf_files = {'1' : 13, '2' : 12, '3' : 10, '4' : 10, '5' : 9, '6' : 9, '7' : 9, '8' : 8, '9' : 7,
             '10' : 7, '11' : 7, '12' : 7, '13' : 6, '14' : 5, '15' : 5, '16' : 5, '17' : 6,
             '18' : 4, '19' : 5, '20' : 3, '21' : 3, '22' : 3}

possible_human_chr = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22']

for chr_num in possible_human_chr:
    num_of_segments = maf_files.get(str(chr_num))
#     print(num_of_segments)
#     print('\n')
    for segment_num in range(1, num_of_segments + 1):
        output_file = open('/scratch/kriall/primate_dfe/to_drive/reformatted_aligns_with_callability/reformatted_alignment_with_callability_chr' + str(chr_num) + '_' + str(segment_num) + '.txt', 'w+')
        output_file.write('position' + '\t' + 'vervet' + '\t' + 'vervet-human' + '\t' + 'vervet-human-marmoset' + '\t' + 'human' + '\t' + 'ancestral' + '\t' + 'diverged' + '\t' + 'callable' + '\n')

        epo_file = '/scratch/kriall/primate_dfe/13_primates_EPO/reformatted_aligns/reformatted_alignment_chr' + str(chr_num) + '_' + str(segment_num) + '.txt'
        with open(epo_file, 'r') as f:
            content = f.readlines()[1 :]
        
        total_num_lines = len(content)
        for i in range(0, total_num_lines):
            split = content[i].strip('\n').split()
            position = split[0]
            vervet_chr_num = position[3 : position.find('_')]
            vervet_bp_position = int(position[position.find('_') + 1 : ])
            v_nuc = split[1]
            v_h_nuc = split[2]
            v_h_m_nuc = split[3]
            h_nuc = split[4]
             
            seq_a = v_h_nuc.upper()
            seq_b = v_nuc.upper()
            seq_c = v_h_m_nuc.upper()
            
            if seq_a == 'X' and seq_c == 'X':
                ancestral = '.'
            elif seq_a == '-' and seq_c == '-':
                ancestral = '-'
            elif seq_a == seq_b == seq_c:
                ancestral = seq_a.upper()
            elif seq_c == 'X' and seq_a == seq_b:
                ancestral = seq_a.lower()
            elif seq_b == '-' and seq_a == seq_c:
                ancestral = seq_a.lower()
            elif (seq_a == seq_b and seq_a != seq_c) or (seq_a == seq_c and seq_a != seq_b):
                ancestral = seq_a.lower()
            elif seq_a != seq_b and seq_a != seq_c:
                ancestral = 'N'
                
            if ancestral.lower() in ['a', 't', 'g', 'c']:
                if v_nuc.lower() != ancestral.lower():
                    is_diverged = '1'
                elif v_nuc.lower() == ancestral.lower():
                    is_diverged = '0'
            else:
                is_diverged = 'N'
                
            callability_value = all_mask_data[vervet_chr_num][vervet_bp_position - 1]
            if callability_value == '0':
                is_callable = 'Y'
            else:
                is_callable = 'N'
                
            to_join_list = [position, v_nuc, v_h_nuc, v_h_m_nuc, h_nuc, ancestral, is_diverged, is_callable]
            output_line = '\t'.join(to_join_list) + '\n'
            output_file.write(output_line)
                        
        output_file.close()

