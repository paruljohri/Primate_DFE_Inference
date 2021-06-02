


#This is a Python script that scans through all of the post alignment files with the ancestral allele already 
#calculated and the callability of that site already calculated. Any sites that fall within any of
#the 5' or 3' intergenic regions flanking the 248 identified homologous vervet exons are output. This is done for 
#two reasons. The first is that, unlike for the exonic regions, if part or all of a flanking intergenic region 
#is missing alignment data, the entire exon set does not get immediately thrown out. Thus, we need to know what
#sites are or are not present in the alignment. This requires that all sites in a region be read in. Second,
#rather than read through all of the reformatted alignment files which would take forever and consume an 
#enormous amount of memory, this file is created to be much easier and quicker to read in.

#This script is currently designed to be run on a cluster and will require path changes on the user's part to the variables:
#'exon_file' with input file "248_homo_exons_reference_with_intergenic_lengths.txt"
#'output_file' with output file "pull_out_intergenic_aligned_sites.txt"
#'reformatted_alignment_file' with input file 'reformatted_alignment_with_callability_chr' + str(chr_num) + '_' + str(segment_num) + '.txt'




maf_files = {'1' : 13, '2' : 12, '3' : 10, '4' : 10, '5' : 9, '6' : 9, '7' : 9, '8' : 8, '9' : 7,
             '10' : 7, '11' : 7, '12' : 7, '13' : 6, '14' : 5, '15' : 5, '16' : 5, '17' : 6,
             '18' : 4, '19' : 5, '20' : 3, '21' : 3, '22' : 3}

possible_human_chr = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22']

human_chr_match_windows_dict = {}
for chr_num in possible_human_chr:
    human_chr_match_windows_dict[chr_num] = []
    
exon_file = ('/scratch/kriall/primate_dfe/13_primates_EPO/248_homo_exons_reference_with_intergenic_lengths.txt')
with open(exon_file, 'r') as f:
    exon_data = f.readlines()[1 : ]
for exon in exon_data:
    split = exon.split()
    human_chr = split[2].replace('chr', '')
    vervet_chr = split[34].replace('chr', '')
    vervet_strand = split[35]

    start_5p = int(split[41]) 
    end_5p = int(split[42])
    start_3p = int(split[44])
    end_3p =int(split[45])
    
    build_window_5p = [vervet_chr, start_5p, end_5p]
    build_window_3p = [vervet_chr, start_3p, end_3p]
    human_chr_match_windows_dict[human_chr].append(build_window_5p)
    human_chr_match_windows_dict[human_chr].append(build_window_3p)

output_file = open('/scratch/kriall/primate_dfe/13_primates_EPO/pull_out_intergenic_aligned_sites.txt', 'w+')
output_file.write('position' + '\t' + 'vervet' + '\t' + 'vervet-human' + '\t' + 'vervet-human-marmoset' + '\t' + 'human' + '\t' + 'ancestral' + '\t' + 'diverged' + '\t' + 'callable' + '\n')

for chr_num in possible_human_chr:
    num_of_segments = maf_files.get(str(chr_num))
#     print(num_of_segments)
#     print('\n')
    for segment_num in range(1, num_of_segments + 1):
        reformatted_alignment_file = '/scratch/kriall/primate_dfe/to_drive/reformatted_aligns_with_callability_unzipped/reformatted_alignment_with_callability_chr' + str(chr_num) + '_' + str(segment_num) + '.txt'

        with open(reformatted_alignment_file, 'r') as f:
            content = f.readlines()[1 : ]
        
        for line in content:
            split = line.strip('\n').split()
            total_position = split[0]
            vervet_chr = total_position[3 : total_position.find('_')]
            bp_position = int(total_position[total_position.find('_') + 1 : ])
            
            windows_to_work_with = human_chr_match_windows_dict[chr_num]
            for window in windows_to_work_with:
                window_chr = window[0]
                start = window[1]
                end = window[2]
                
                if vervet_chr == window_chr:
                    if start <= bp_position <= end:
                        print(bp_position)
                        output_file.write(line)

