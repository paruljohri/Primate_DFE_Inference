


#This is a Python script designed to output various characteristics about all of the aligned 
#regions that encompass the 248 homologous vervet exon regions. It works by taking the start 
#and end position of the exon, finds which vervet sequence(s) these exons are a part of in the 
#alignment files (13_primates.epo files), and then outputs the start and end positions of those
#vervet regions, the start and end positions of the corresponding aligned human sequences, and 
#the number of human, vervet, and marmoset listings in that alignment group. Note that some of
#the information output by this script goes unused, but I see no region to spend time changing 
#it if I can easily extract the particular things I need.

#This script is currently designed to be run on a cluster and will require path changes on the user's part to the variables:
#'output_file' with output file "248_exon_alignment_check.txt"
#'exon_data' with input file "248_homo_exons_reference_with_intergenic_lengths.txt"
# 'epo_file' with input files '13_primates.epo.' + str(human_chr) + '_' + str(segment_num) + '.maf'




import itertools

maf_files = {'1' : 13, '2' : 12, '3' : 10, '4' : 10, '5' : 9, '6' : 9, '7' : 9, '8' : 8, '9' : 7,
             '10' : 7, '11' : 7, '12' : 7, '13' : 6, '14' : 5, '15' : 5, '16' : 5, '17' : 6,
             '18' : 4, '19' : 5, '20' : 3, '21' : 3, '22' : 3}


possible_human_chr = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22']
possible_vervet_chr = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29']

output_file = open('/scratch/kriall/primate_dfe/13_primates_EPO/248_exon_alignment_check.txt', 'w+')
output_file.write('human_range' + '\t' 'vervet_range' + '\t' + '#_blocks' + '\t' + 'overlap?' + '\t' + 'human_block(s)_positions' + '\t' + 'vervet_block(s)_positions' + '\t' + 'gene_duplicates' + '\t' + 'source_listings' + '\n')

exon_data = ('/scratch/kriall/primate_dfe/13_primates_EPO/248_homo_exons_reference_with_intergenic_lengths.txt')
with open(exon_data, 'r') as f:
    content = f.readlines()[1 : ]
for exon in content:
    split = exon.split()
    human_chr = int(split[2].replace('chr', ''))
    human_strand = split[3]
    h_start = int(split[6])
    h_end = int(split[7])

    vervet_chr = int(split[34].replace('chr', ''))
    vervet_strand = split[35]
    v_start = int(split[38])
    v_end = int(split[39])

        
    information_blocks = []
    human_is_overlap = False
    human_duplicates_exist = False
    human_chr_oddity = False
    vervet_is_overlap = False
    vervet_duplicates_exist = False
    vervet_chr_oddity = False
    
#     print(human_chr)
#     print(human_strand)
#     print(h_start)
#     print(h_end)
#     print('\n')
#     print(vervet_chr)
#     print(vervet_strand)
#     print(v_start)
#     print(v_end)
    
    num_of_segments = maf_files.get(str(human_chr))
    print(num_of_segments)
    print('\n')
    
    for segment_num in range(1, num_of_segments + 1):
        needed_sources = ['vervet', 'vervet-human', 'vervet-human-marmoset', 'human']
        sources_dictionary = {}


#         epo_file = 'E:/13_primates_EPO/unzipped/13_primates.epo.' + str(human_chr) + '_' + str(segment_num) + '.maf'
        epo_file = '/scratch/kriall/primate_dfe/13_primates_EPO/raw_aligns/13_primates.epo.' + str(human_chr) + '_' + str(segment_num) + '.maf'
        with open(epo_file, 'r') as f:
            content = f.readlines()

        paragraph_pairs = []
        total_line_num = len(content)

        start_found = False
        end_found = False
        for i in range(0, total_line_num):
            split = content[i].split()
            if len(split) >= 2 and split[1] == 'tree:':
                start_line = i
                start_found = True
            elif len(split) == 0:
                end_line = i
                end_found = True

            if start_found == True and end_found == True:
                pair = [start_line, end_line]
                paragraph_pairs.append(pair)
                start_found = False
                end_found = False
#         print(paragraph_pairs)
#         print(len(paragraph_pairs))

        for pair in paragraph_pairs:
            human_alignment_match_found = False
            vervet_alignment_match_found = False
            human_alignment_blocks = []
            vervet_alignment_blocks = []
            listings = []
            gene_duplicates = []
           
            start_line = pair[0]
            end_line = pair[1]
            tree_seq = content[start_line].split()[2]
            id_num = content[start_line + 1].split()[2]

            for j in range(start_line + 3, end_line):
                generic_split = content[j].split()
                generic_source_full = generic_split[1]
                

                if 'homo_sapiens' in generic_source_full:
                    ###Find human sequence
                    h_key = 'human'
                    h_split = content[j].split()
#                     print(h_split[0:6])
                    h_listing = h_split[1]
                    h_chr = h_listing[h_listing.find('.') + 1 : ]
                    h_given_start = int(h_split[2])
                    h_align_len = int(h_split[3])
                    h_strand = h_split[4]
                    h_chr_len = int(h_split[5])

                    listings.append(h_listing + h_strand)
                    gene_duplicates.append('H')
                    
                    if h_strand == '+':
                        h_start = h_given_start + 1
                        h_end = h_start + h_align_len - 1
                    elif h_strand == '-':
                        h_end = h_chr_len - h_given_start
                        h_start = h_end - h_align_len + 1

                    if h_start <= h_start <= h_end or h_start <= h_end <= h_end:
#                         print(h_start)
#                         print(h_start)
#                         print(h_end)
#                         print(h_end)
#                         print(segment_num)

                        human_alignment_match_found = True
                        human_alignment_blocks = [str(h_start), str(h_end)]

                elif 'chlorocebus_sabaeus' in generic_source_full:
                    ###Find vervet sequence
                    v_key = 'vervet'
                    v_split = content[j].split()
#                     print(v_split[0:6])
                    v_listing = v_split[1]
                    v_chr = v_listing[v_listing.find('.') + 1 : ]
                    v_given_start = int(v_split[2])
                    v_align_len = int(v_split[3])
                    v_strand = v_split[4]
                    v_chr_len = int(v_split[5])

                    listings.append(v_listing + v_strand)
                    gene_duplicates.append('V')

                                        
                    if v_strand == '+':
                        v_start = v_given_start + 1
                        v_end = v_start + v_align_len - 1
                    if v_strand == '-':
                        v_end = v_chr_len - v_given_start
                        v_start = v_end - v_align_len + 1

                    if v_start <= v_start <= v_end or v_start <= v_end <= v_end:
#                         print(v_start)
#                         print(v_start)
#                         print(v_end)
#                         print(v_end)
#                         print(segment_num)

                        vervet_alignment_match_found = True
                        vervet_alignment_blocks = [str(v_start), str(v_end)]
                        
                elif 'callithrix_jacchus' in generic_source_full:
                    ###Find vervet-human-marmoset ancestral sequence
                    c_key = 'marmoset'
                    c_split = content[j].split()
#                     print(v_split[0:6])
                    c_listing = c_split[1]
                    c_chr = c_listing[c_listing.find('.') + 1 : ]
                    c_given_start = int(c_split[2])
                    c_align_len = int(c_split[3])
                    c_strand = c_split[4]
                    c_chr_len = int(c_split[5])
                    
                    listings.append(c_listing + c_strand)
                    gene_duplicates.append('C')
   
            if human_alignment_match_found == True:
                information_blocks.append([human_alignment_blocks, vervet_alignment_blocks, gene_duplicates, listings])
   
    for boop in information_blocks:
        print(boop)
    print(information_blocks)
        
    is_overlap = False
    if len(information_blocks) > 1:
        all_alignment_blocks = []
        for i in range(0, len(information_blocks)):
            all_alignment_blocks.append(information_blocks[i][0])
        combinations = list(itertools.combinations(all_alignment_blocks, 2))
        for combo in combinations:
            set1 = range(int(combo[0][0]) , int(combo[0][1]) + 1)
            set2 = range(int(combo[1][0]) , int(combo[1][1]) + 1)
            overlap = set(set1).intersection(set2)
            if len(overlap) > 0:
                is_overlap = True
        if is_overlap == True:
            overlap = 'overlap'
        elif is_overlap == False:
            overlap = 'continuous'
    else:
        overlap = 'NA'
    
    human_range = str(h_start) + '-' + str(h_end)
    vervet_range = str(v_start) + '-' + str(v_end)
    printout = [human_range, vervet_range, str(len(information_blocks)), overlap]
    for x in range(0, 4):
        temp = ''
        for y in range(0, len(information_blocks)):
            join = ', '.join(information_blocks[y][x])
            item = '(' + join + ')'
            temp = temp + item + ', '
        temp = temp[0 : -2]
        printout.append(temp)
    final_printout = '\t'.join(printout) + '\n'
    output_file.write(final_printout)
    
output_file.close()
    

