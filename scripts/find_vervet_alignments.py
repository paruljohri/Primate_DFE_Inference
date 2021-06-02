#A Python script capable of reading through a 13 primates EPO file. It uses the Phylo module from the 
#BioPython package, documentation for which can be found at https://biopython.org/wiki/Phylo.

#This script iterates through each of the alignment files. For every alignment file, it reads it in 
#and first determines the line numbers that delimit the various alignment groups in that particular 
#alignment file. Using a module from BioPython called Phylo, it then goes through each of these alignment 
#groups and builds the phylogenetic tree for that group. It then determines if a human sequence, a vervet 
#sequence, and a marmoset sequence exist in that tree. If multiple listings for either the human or vervet
#are found, the first listing it finds that contains a sequence from an actual chromosome is taken. Now 
#knowing which species are present in this group, the names of the vervet-human ancestral sequence 
#(assuming both human and vervet are present) and the vervet-human-marmoset ancestral sequence 
#(assuming all three species are present) in the tree are found by finding the most recent common ancestor
#between them. With the listing names of the human sequence, vervet sequence, vervet-human ancestral 
#sequence, and vervet-human-marmoset ancestral sequence is this alignment group, it finds these four 
#sequences and aligns them. If the vervet sequence is on the negative strand, it reverse complements 
#all four sequences to maintain the alignment. Any sequences that may be missing in the alignment group 
#are simply given a sequence of “XXXXXX…” that is the same length of the other sequences. Finally, 
#with the four sequences aligned, the script outputs every vervet position given on an individual 
#line along with the vervet nucleotide and other three aligned nucleotides at that position.

#This script is currently designed to be run locally, and will require path changes on the user's part to the variables:
#'output_file' with output files 'reformatted_alignment_chr' + str(chr_num) + '_' + str(segment_num) + '.txt'
#'epo_file' with input files '13_primates.epo.' + str(chr_num) + '_' + str(segment_num) + '.maf'


from Bio import Phylo
from io import StringIO
import sys

def find_reverse_complement(sequence):
    rev_comp = ''
    reverse_seq = sequence[::-1]
    for site in reverse_seq:
#         print(site)
        if site == 'X':
            complement = site
        elif site == '-':
            complement = site
        elif site == 'N' or site == 'n':
            complement = site
        else:
            if site.isupper() == True:
                if site == 'A':
                    complement = 'T'
                if site == 'T':
                    complement = 'A'
                if site == 'G':
                    complement = 'C'
                if site == 'C':
                    complement = 'G'
            elif site.islower() == True:
                if site == 'a':
                    complement = 't'
                if site == 't':
                    complement = 'a'
                if site == 'g':
                    complement = 'c'
                if site == 'c':
                    complement = 'g'
        rev_comp = rev_comp + complement
    return rev_comp

maf_files = {'1' : 13, '2' : 12, '3' : 10, '4' : 10, '5' : 9, '6' : 9, '7' : 9, '8' : 8, '9' : 7,
             '10' : 7, '11' : 7, '12' : 7, '13' : 6, '14' : 5, '15' : 5, '16' : 5, '17' : 6,
             '18' : 4, '19' : 5, '20' : 3, '21' : 3, '22' : 3}

possible_human_chr = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22']
possible_vervet_chr = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29']

for chr_num in possible_human_chr:
    num_of_segments = maf_files.get(str(chr_num))
    for segment_num in range(1, num_of_segments + 1):
        print(str(chr_num) + ' - ' + str(segment_num))
        print('\n')
        needed_sources = ['vervet', 'vervet-human', 'vervet-human-marmoset', 'human']
        sources_dictionary = {}
        content = []

        output_file = open('E:/13_primates_EPO/reformatted_aligns/reformatted_alignment_chr' + str(chr_num) + '_' + str(segment_num) + '.txt', 'w+')
        output_file_title = 'position'
        for source in needed_sources:
            output_file_title = output_file_title + '\t' + source
        output_file.write(output_file_title + '\n') 

        epo_file = 'E:/13_primates_EPO/unzipped/13_primates.epo.' + str(chr_num) + '_' + str(segment_num) + '.maf'
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
        # print(paragraph_pairs)

        for pair in paragraph_pairs:
            base_pair_count = 'X'
            human_found = False
            vervet_found = False
            marmoset_found = False
            what_is_vervet_strand_sign = 'unknown' 
            for source in needed_sources:
                sources_dictionary[source] = []

            start_line = pair[0]
            end_line = pair[1]
            raw_tree = content[start_line].split()[2]
            id_num = content[start_line + 1].split()[2]

            handle = StringIO(raw_tree)
            read_in_tree = Phylo.read(handle, "newick")
            str_tree = str(read_in_tree)

            is_human_present = False
            is_vervet_present = False
            is_marmoset_present = False
            human_tree_name = 'HUMAN TREE ABSENT'
            human_align_name = 'HUMAN ALIGN ABSENT'
            vervet_tree_name = 'VERVET TREE ABSENT'
            vervet_align_name = 'VERVET ALIGN ABSENT'
            marmoset_tree_name = 'MARMOSET TREE ABSENT'
            marmoset_align_name = 'MARMOSET ALIGN ABSENT'
            v_h_tree_name = 'VH TREE ABSENT'
            v_h_align_name = 'VH ALIGN ABSENT'
            v_h_m_tree_name = 'VHM TREE ABSENT'
            v_h_m_align_name = 'VHM ALIGN ABSENT'

            human_instances = [i for i in range(len(raw_tree)) if raw_tree.startswith('Hsap', i)]
            if len(human_instances) > 0:
                for human_finding in human_instances:
                    trimmed_tree = raw_tree[human_finding : ]
                    find_end_name = trimmed_tree.find('[')
                    human_tree_name = trimmed_tree[ : find_end_name]
                    chromosome_match = '_' + str(chr_num) + '_'
                    if chromosome_match in human_tree_name:
                        is_human_present = True
                        human_align_name = 'homo_sapiens.' + str(chr_num) 
                        break
        #     print(human_tree_name)
        #     print(human_align_name)

            vervet_finding = raw_tree.find('Csab')
            if vervet_finding != -1:
                trimmed_tree = raw_tree[vervet_finding : ]
                find_end_name = trimmed_tree.find('[')
                vervet_tree_name = trimmed_tree[ : find_end_name]
                vervet_chromosome_num = vervet_tree_name[5 : 5 + vervet_tree_name[5 : ].find('_') ]
                vervet_align_name = 'chlorocebus_sabaeus.' + str(vervet_chromosome_num)
                if vervet_chromosome_num in possible_vervet_chr:
                    is_vervet_present = True
        #     print(vervet_tree_name)
        #     print(vervet_align_name)

            marmoset_finding = raw_tree.find('Cjac')
            if marmoset_finding != -1:
                trimmed_tree = raw_tree[marmoset_finding : ]
                find_end_name = trimmed_tree.find('[')
                marmoset_tree_name = trimmed_tree[ : find_end_name]
                marmoset_chromosome_num = marmoset_tree_name[5 : 5 + marmoset_tree_name[5 : ].find('_') ]
                marmoset_align_name = 'callithrix_jacchus.' + str(marmoset_chromosome_num)
                is_marmoset_present = True
        #     print(marmoset_tree_name)
        #     print(marmoset_align_name)

        #     print(is_human_present)
        #     print(is_vervet_present)
        #     print(is_marmoset_present)
        
            if is_human_present == True and is_vervet_present == True:
                v_h_tree_name = str(read_in_tree.common_ancestor({'name' : human_tree_name}, {'name' : vervet_tree_name}))
                trim1 = v_h_tree_name[0 : v_h_tree_name.rfind('_')]
                trim2 = v_h_tree_name[0: trim1.rfind('_')]
                v_h_align_name = 'ancestral_sequences.' + trim2[5 :]
            if is_human_present == True and is_vervet_present == True and is_marmoset_present == True:
                v_h_m_tree_name = str(read_in_tree.common_ancestor({'name' : human_tree_name}, {'name' : vervet_tree_name},  {'name' : marmoset_tree_name}))
                trim1 = v_h_m_tree_name[0 : v_h_m_tree_name.rfind('_')]
                trim2 = v_h_m_tree_name[0: trim1.rfind('_')]
                v_h_m_align_name = 'ancestral_sequences.' + trim2[5 :]
                
        #     print(v_h_tree_name)
        #     print(v_h_align_name)
        #     print(v_h_m_tree_name)
        #     print(v_h_m_align_name)
        #     print('\n')


            for j in range(start_line + 3, end_line):
                generic_split = content[j].split()
                generic_source_full = generic_split[1]
                generic_align_length = len(generic_split[6])

                if vervet_align_name in generic_source_full and is_vervet_present == True and vervet_found == False:
                    ###Find vervet sequence
                    v_key = 'vervet'
                    v_split = content[j].split()
                    v_source_full = v_split[1]
                    v_chr_num = 'chr' + v_source_full[v_source_full.find('.') + 1 : ]
                    v_given_start = int(v_split[2])
                    v_strand = v_split[4]
                    v_chr_len = int(v_split[5])
                    v_sequence = v_split[6]
                    v_align_len = len(v_sequence.replace('-', ''))
                    sources_dictionary[v_key] = [v_chr_num, v_strand, v_sequence]

                    if v_strand == '+':
                        what_is_vervet_strand_sign = '+'
                        v_start_position = v_given_start + 1
                        v_end_position = v_start_position + v_align_len - 1
                        base_pair_count = v_start_position
                    if v_strand == '-':
                        what_is_vervet_strand_sign = '-'
                        v_end_position = v_chr_len - v_given_start
                        v_start_position = v_end_position - v_align_len + 1
                        base_pair_count = v_start_position

                    if v_chr_num.replace('chr', '') in possible_vervet_chr:
                        vervet_found = True

                elif v_h_align_name in generic_source_full:
                    ###Find vervet-human ancestral sequence
                    v_h_key = 'vervet-human'
                    v_h_split = content[j].split()
                    v_h_source_full = v_h_split[1]
                    v_h_sequence = v_h_split[6]
                    sources_dictionary[v_h_key] = ['NA', 'NA', v_h_sequence]

                elif v_h_m_align_name in generic_source_full: 
                    ###Find vervet-human-marmoset ancestral sequence
                    v_h_m_key = 'vervet-human-marmoset'
                    v_h_m_split = content[j].split()
                    v_h_m_source_full = v_h_m_split[1]
                    v_h_m_sequence = v_h_m_split[6]
                    sources_dictionary[v_h_m_key] = ['NA', 'NA', v_h_m_sequence]

                elif human_align_name in generic_source_full and human_found == False:
                    ###Find human sequence
                    h_key = 'human'
                    h_split = content[j].split()
                    h_source_full = h_split[1]
                    h_chr_num = 'chr' + h_source_full[h_source_full.find('.') + 1 : ]

                    h_sequence = h_split[6]
                    sources_dictionary[h_key] = ['NA', 'NA', h_sequence]
                    if h_chr_num.replace('chr', '') in possible_human_chr:
                        human_found = True


            for source in needed_sources:
                if sources_dictionary[source] == []:
                    sources_dictionary[source] = ['X', 'X', 'X' * generic_align_length]

            if what_is_vervet_strand_sign == '-':
                for source in sources_dictionary:
                    info = sources_dictionary[source]
                    sequence = info[2]
                    rev_comp_sequence = find_reverse_complement(sequence)
                    info[2] = rev_comp_sequence
                    sources_dictionary[source] = info

            if base_pair_count != 'X':
                current_pos = base_pair_count 
                for site in range(0, generic_align_length):
                    total_entry = ''
                    base_pair_position = 'NA'
                    for source in needed_sources:
                        chromosome = sources_dictionary[source][0]
                        strand = sources_dictionary[source][1]
                        nucleotide = sources_dictionary[source][2][site]
        #                 print(current_pos)
        #                 print(chromosome)
        #                 print(nucleotide)

                        if source == 'vervet' and (nucleotide != '-' and nucleotide != 'X'):
                            base_pair_position = chromosome + '_' + str(current_pos)
                            current_pos += 1                              
                        total_entry = total_entry + '\t' + nucleotide
                    total_entry = base_pair_position + total_entry
                    if 'NA' not in total_entry:
                        output_file.write(total_entry + '\n')

        output_file.close()

