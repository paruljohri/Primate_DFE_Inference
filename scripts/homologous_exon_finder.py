#This is a Python script designed to find specific chimp and vervet exons that are contained within Susanne’s given 
#homologous regions, thus allowing us to find sets of homologous exons wherein the human, chimp, and vervet exons are 
#certain to be homologous. It first reads in the positions given by Susanne in the “human_exons_human.txt”, 
#“human_exons_pt2.txt”, and “human_exons_chlSab2.txt” files. It then reads in all chimp and vervet exons from their 
#respective UCSC exon annotation files. With all of this information, it begins to go through each chimp and vervet 
#region given by Susanne. This region is then checked against every single exon in the same chromosome according to the 
#annotation looking for exons with either the same start position (in actuality with a start of 1bp further downstream 
#likely due to an indexing change somewhere along the line, so, for example, 53 in the annotation would be a match to 
#a given position by Susanne of 52) or the same end position (which took exact positional matches). This process is done 
#independently for both chimps and vervets, giving you the exon that is found for each of those regions. For both species,
#there were 444 regions from Susanne to check this against, and thus a maximum of 444 possible matching exons in both 
#chimps and vervets. However, as you will see, not all 444 regions have matches for both species. After all matching 
#regions for both species are found, the script checks to see which regions have matches in both species to the homologous
#human region, and then outputs these. The way this script is set up, it finds 248 exons that are homologous between all
#three species, but alterations to the finding process could potentially change this amount.

#This script is currently designed to be run locally, and will require path changes on the user's part to the variables:
#'chimp_exon_positions_table' with input file "human_exons_pt2.txt"
#'vervet_exon_positions_table' with input file "human_exons_chlSab2.txt"
#'human_exon_positions_table' with input file  "homolgous_human_exons.txt"
#'chimp_exon_table' with input file "panTro2.ensGenes.txt"
#'vervet_exon_table' with input file "chlSab2.ensGenes.txt"
#'chimp_file' with output file "chimp_homologous_exons.txt"
#'chimp_file_relevant' with output file "chimp_homologous_exons_with_relevant_exons.txt"
#'chimp_file_formatted' with output file "chimp_homologous_exons_formatted.txt"
#'vervet_file' with output file "vervet_homologous_exons.txt"
#'vervet_file_relevant' with output file "vervet_homologous_exons_with_relevant_exons.txt"
#'vervet_file_formatted' with output file "vervet_homologous_exons_formatted.txt"
#'all_primates_file' with output file "all_homo_exons_reference_exons_only.txt"


#This section reads in all of the supposedly homologous regions given by Susanne.
chimp_exon_positions_table = 'C:/Users/Kellen/Downloads/primate_dfe_project/chimp_and_vervet_data/human_exons_pt2.txt'
with open(chimp_exon_positions_table, 'r+') as f:
    chimp_exon_positions = f.readlines()

vervet_exon_positions_table = 'C:/Users/Kellen/Downloads/primate_dfe_project/chimp_and_vervet_data/human_exons_chlSab2.txt'
with open(vervet_exon_positions_table, 'r+') as f:
    vervet_exon_positions = f.readlines()
    
human_exon_positions_table = 'C:/Users/Kellen/Downloads/primate_dfe_project/chimp_and_vervet_data/homolgous_human_exons.txt'
with open(human_exon_positions_table, 'r+') as f:
    human_exon_positions = f.readlines()
    del human_exon_positions[0]


#This section reads in the exons contained in the chimp UCSC exon annotation.
chimp_exon_table = 'C:/Users/Kellen/Downloads/primate_dfe_project/chimp_and_vervet_data/panTro2.ensGenes/panTro2.ensGenes.txt'
with open(chimp_exon_table, 'r+') as f:
    chimp_exon_data = f.readlines()
del chimp_exon_data[0]

chimp_chr1_exon_info = []
chimp_chr2a_exon_info = []
chimp_chr2b_exon_info = []
chimp_chr3_exon_info = []
chimp_chr4_exon_info = []
chimp_chr5_exon_info = []
chimp_chr6_exon_info = []
chimp_chr7_exon_info = []
chimp_chr8_exon_info = []
chimp_chr9_exon_info = []
chimp_chr10_exon_info = []
chimp_chr11_exon_info = []
chimp_chr12_exon_info = []
chimp_chr13_exon_info = []
chimp_chr14_exon_info = []
chimp_chr15_exon_info = []
chimp_chr16_exon_info = []
chimp_chr17_exon_info = []
chimp_chr18_exon_info = []
chimp_chr19_exon_info = []
chimp_chr20_exon_info = []
chimp_chr21_exon_info = []
chimp_chr22_exon_info = []

index = 0
for a in chimp_exon_data:
    exons = a.split()
    transcript_name = exons[1]
    exon_name = exons[12]
    chr_num = exons[2]
    strand = exons[3]
    exon_starts = exons[9].split(',')[0 : -1]
    exon_ends = exons[10].split(',')[0 : -1]
    exon_info = [exon_starts, exon_ends, index, exon_name, chr_num, strand, transcript_name]
#     print(exon_info)
    
    if chr_num == 'chr1':
        chimp_chr1_exon_info.append(exon_info)
    elif chr_num == 'chr2a':
        chimp_chr2a_exon_info.append(exon_info)
    elif chr_num == 'chr2b':
        chimp_chr2b_exon_info.append(exon_info)
    elif chr_num == 'chr3':
        chimp_chr3_exon_info.append(exon_info)
    elif chr_num == 'chr4':
        chimp_chr4_exon_info.append(exon_info)
    elif chr_num == 'chr5':
        chimp_chr5_exon_info.append(exon_info)
    elif chr_num == 'chr6':
        chimp_chr6_exon_info.append(exon_info)
    elif chr_num == 'chr7':
        chimp_chr7_exon_info.append(exon_info)
    elif chr_num == 'chr8':
        chimp_chr8_exon_info.append(exon_info)
    elif chr_num == 'chr9':
        chimp_chr9_exon_info.append(exon_info)
    elif chr_num == 'chr10':
        chimp_chr10_exon_info.append(exon_info)
    elif chr_num == 'chr11':
        chimp_chr11_exon_info.append(exon_info)
    elif chr_num == 'chr12':
        chimp_chr12_exon_info.append(exon_info)
    elif chr_num == 'chr13':
        chimp_chr13_exon_info.append(exon_info)
    elif chr_num == 'chr14':
        chimp_chr14_exon_info.append(exon_info)
    elif chr_num == 'chr15':
        chimp_chr15_exon_info.append(exon_info)
    elif chr_num == 'chr16':
        chimp_chr16_exon_info.append(exon_info)
    elif chr_num == 'chr17':
        chimp_chr17_exon_info.append(exon_info)
    elif chr_num == 'chr18':
        chimp_chr18_exon_info.append(exon_info)
    elif chr_num == 'chr19':
        chimp_chr19_exon_info.append(exon_info)
    elif chr_num == 'chr20':
        chimp_chr20_exon_info.append(exon_info)
    elif chr_num == 'chr21':
        chimp_chr21_exon_info.append(exon_info)
    elif chr_num == 'chr22':
        chimp_chr22_exon_info.append(exon_info)
    
    index += 1


#This section reads in the exons contained in the vervet UCSC exon annotation.
vervet_exon_table = 'C:/Users/Kellen/Downloads/primate_dfe_project/chimp_and_vervet_data/chlSab2.ensGenes/chlSab2.ensGenes.txt'
with open(vervet_exon_table, 'r+') as f:
    vervet_exon_data = f.readlines()
del vervet_exon_data[0]

vervet_chr1_exon_info = []
vervet_chr2_exon_info = []
vervet_chr3_exon_info = []
vervet_chr4_exon_info = []
vervet_chr5_exon_info = []
vervet_chr6_exon_info = []
vervet_chr7_exon_info = []
vervet_chr8_exon_info = []
vervet_chr9_exon_info = []
vervet_chr10_exon_info = []
vervet_chr11_exon_info = []
vervet_chr12_exon_info = []
vervet_chr13_exon_info = []
vervet_chr14_exon_info = []
vervet_chr15_exon_info = []
vervet_chr16_exon_info = []
vervet_chr17_exon_info = []
vervet_chr18_exon_info = []
vervet_chr19_exon_info = []
vervet_chr20_exon_info = []
vervet_chr21_exon_info = []
vervet_chr22_exon_info = []
vervet_chr23_exon_info = []
vervet_chr24_exon_info = []
vervet_chr25_exon_info = []
vervet_chr26_exon_info = []
vervet_chr27_exon_info = []
vervet_chr28_exon_info = []
vervet_chr29_exon_info = []

index = 0
for a in vervet_exon_data:
    exons = a.split()
    transcript_name = exons[1]
    exon_name = exons[12]
    chr_num = exons[2]
    strand = exons[3]
    exon_starts = exons[9].split(',')[0 : -1]
    exon_ends = exons[10].split(',')[0 : -1]
    exon_info = [exon_starts, exon_ends, index, exon_name, chr_num, strand, transcript_name]
#     print(exon_info)
    
    if chr_num == 'chr1':
        vervet_chr1_exon_info.append(exon_info)
    elif chr_num == 'chr2':
        vervet_chr2_exon_info.append(exon_info)
    elif chr_num == 'chr3':
        vervet_chr3_exon_info.append(exon_info)
    elif chr_num == 'chr4':
        vervet_chr4_exon_info.append(exon_info)
    elif chr_num == 'chr5':
        vervet_chr5_exon_info.append(exon_info)
    elif chr_num == 'chr6':
        vervet_chr6_exon_info.append(exon_info)
    elif chr_num == 'chr7':
        vervet_chr7_exon_info.append(exon_info)
    elif chr_num == 'chr8':
        vervet_chr8_exon_info.append(exon_info)
    elif chr_num == 'chr9':
        vervet_chr9_exon_info.append(exon_info)
    elif chr_num == 'chr10':
        vervet_chr10_exon_info.append(exon_info)
    elif chr_num == 'chr11':
        vervet_chr11_exon_info.append(exon_info)
    elif chr_num == 'chr12':
        vervet_chr12_exon_info.append(exon_info)
    elif chr_num == 'chr13':
        vervet_chr13_exon_info.append(exon_info)
    elif chr_num == 'chr14':
        vervet_chr14_exon_info.append(exon_info)
    elif chr_num == 'chr15':
        vervet_chr15_exon_info.append(exon_info)
    elif chr_num == 'chr16':
        vervet_chr16_exon_info.append(exon_info)
    elif chr_num == 'chr17':
        vervet_chr17_exon_info.append(exon_info)
    elif chr_num == 'chr18':
        vervet_chr18_exon_info.append(exon_info)
    elif chr_num == 'chr19':
        vervet_chr19_exon_info.append(exon_info)
    elif chr_num == 'chr20':
        vervet_chr20_exon_info.append(exon_info)
    elif chr_num == 'chr21':
        vervet_chr21_exon_info.append(exon_info)
    elif chr_num == 'chr22':
        vervet_chr22_exon_info.append(exon_info)
    elif chr_num == 'chr23':
        vervet_chr23_exon_info.append(exon_info)
    elif chr_num == 'chr24':
        vervet_chr24_exon_info.append(exon_info)
    elif chr_num == 'chr25':
        vervet_chr25_exon_info.append(exon_info)
    elif chr_num == 'chr26':
        vervet_chr26_exon_info.append(exon_info)
    elif chr_num == 'chr27':
        vervet_chr27_exon_info.append(exon_info)
    elif chr_num == 'chr28':
        vervet_chr28_exon_info.append(exon_info)
    elif chr_num == 'chr29':
        vervet_chr29_exon_info.append(exon_info)
    
    index += 1



#This function takes in the region given by Susanne, looks through all of exons from the annotation on the same 
#chromosome to see if any have starting or ending positions that match the start or end of Susanne's region, and 
#then marks and outputs the matching exon if one is found.
def find_orthologs(primate, primate_chrXXX_exon_info, primate_chrXXX_hits, primate_start, primate_end, susanne_pos_index):
    start_offset = 1 #As in, how many base pairs is the window Susanne gave ahead of the corresponding UCSC annotation window
    end_offset = 0
    for exon in primate_chrXXX_exon_info:
        transcript_name = exon[6]
        num_exons = len(exon[0])
        for pos in range(0, num_exons):
            exon_start_pos = int(exon[0][pos])
            exon_end_pos = int(exon[1][pos])
            already_counted = False
            for window in primate_chrXXX_hits:
                counted_start = window[0]
                counted_end = window[1]
                if primate_start == counted_start or primate_end == counted_end:
                    already_counted = True
            if (primate_start == (exon_start_pos + start_offset) or primate_end == (exon_end_pos + end_offset)) and already_counted == False:
                if primate_start == (exon_start_pos + start_offset):
                    exon_index = exon[0].index(str(exon_start_pos))
                elif primate_end == (exon_end_pos + end_offset):
                    exon_index = exon[1].index(str(exon_end_pos))
                primate_chrXXX_hits.append([primate_start, primate_end])
                annotation_line = int(exon[2])

                if primate == 'chimp':
                    chimp_hit_indexes.append(susanne_pos_index)
                    output = chimp_exon_data[annotation_line]
                    output_relevant = output.strip('\n') + '\t' + str(exon_start_pos) + '\t' + str(exon_end_pos) + '\n'
                    output_formatted = str(susanne_pos_index) + '\t' + transcript_name + '\t'  + exon[3] + '_ex' + str(exon_index + 1) + '/' + str(num_exons) + '\t' + exon[4] + '\t' + exon[5] + '\t' + str(exon_start_pos) + '\t' + str(exon_end_pos) + '\n'
                    chimp_table.write(output)
                    chimp_table_relevant.write(output_relevant)
                    chimp_table_formatted.write(output_formatted)
                if primate == 'vervet':
                    vervet_hit_indexes.append(susanne_pos_index)
                    output = vervet_exon_data[annotation_line]
                    output_relevant = output.strip('\n') + '\t' + str(exon_start_pos) + '\t' + str(exon_end_pos) + '\n'
#                     output_formatted = str(susanne_pos_index) + '\t' + exon[3] + '_ex' + str(exon_index + 1) + '\t' + exon[4] + '\t' + exon[5] + '\t' + str(exon_start_pos) + '\t' + str(exon_end_pos) + '\n'
                    output_formatted = str(susanne_pos_index) + '\t' + transcript_name + '\t'  + exon[3] + '_ex' + str(exon_index + 1) + '/' + str(num_exons) + '\t' + exon[4] + '\t' + exon[5] + '\t' + str(exon_start_pos) + '\t' + str(exon_end_pos) + '\n'
                    vervet_table.write(output)
                    vervet_table_relevant.write(output_relevant)
                    vervet_table_formatted.write(output_formatted)



#This section creates all of the new output files to be written to, and then loops through each of Susanne's 
#regions for both chimps and vervets, sending that region data to the "find_orthologos" function to look for an 
#exon that matches.
chimp_file = 'C:/Users/Kellen/Downloads/primate_dfe_project/chimp_and_vervet_data/chimp_homologous_exons.txt'
chimp_file_relevant = 'C:/Users/Kellen/Downloads/primate_dfe_project/chimp_and_vervet_data/chimp_homologous_exons_with_relevant_exons.txt'
chimp_file_formatted = 'C:/Users/Kellen/Downloads/primate_dfe_project/chimp_and_vervet_data/chimp_homologous_exons_formatted.txt'
vervet_file = 'C:/Users/Kellen/Downloads/primate_dfe_project/chimp_and_vervet_data/vervet_homologous_exons.txt'
vervet_file_relevant = 'C:/Users/Kellen/Downloads/primate_dfe_project/chimp_and_vervet_data/vervet_homologous_exons_with_relevant_exons.txt'
vervet_file_formatted = 'C:/Users/Kellen/Downloads/primate_dfe_project/chimp_and_vervet_data/vervet_homologous_exons_formatted.txt'

chimp_table = open(chimp_file, 'w+')
chimp_table_relevant = open(chimp_file_relevant, 'w+')
chimp_table_formatted = open(chimp_file_formatted, 'w+')
vervet_table = open(vervet_file, 'w+')
vervet_table_relevant = open(vervet_file_relevant, 'w+')
vervet_table_formatted = open(vervet_file_formatted, 'w+')

chimp_table.write('#bin	name	chrom	strand	txStart	txEnd	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds	score	name2	cdsStartStat	cdsEndStat	exonFrames'+'\n')
chimp_table_relevant.write('#bin	name	chrom	strand	txStart	txEnd	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds	score	name2	cdsStartStat	cdsEndStat	exonFrames	relevantExonStart	relevantExonEnd'+'\n')
chimp_table_formatted.write('susannes_file_line_index	transcript	exon	chrom	exonStart	exonEnd'+'\n')
vervet_table.write('#bin	name	chrom	strand	txStart	txEnd	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds	score	name2	cdsStartStat	cdsEndStat	exonFrames'+'\n')
vervet_table_relevant.write('#bin	name	chrom	strand	txStart	txEnd	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds	score	name2	cdsStartStat	cdsEndStat	exonFrames	relevantExonStart	relevantExonEnd'+'\n')
vervet_table_formatted.write('susannes_file_line_index	transcript	exon	chrom	exonStart	exonEnd'+'\n')

chimp_chr1_hits = []
chimp_chr2a_hits = []
chimp_chr2b_hits = []
chimp_chr3_hits = []
chimp_chr4_hits = []
chimp_chr5_hits = []
chimp_chr6_hits = []
chimp_chr7_hits = []
chimp_chr8_hits = []
chimp_chr9_hits = []
chimp_chr10_hits = []
chimp_chr11_hits = []
chimp_chr12_hits = []
chimp_chr13_hits = []
chimp_chr14_hits = []
chimp_chr15_hits = []
chimp_chr16_hits = []
chimp_chr17_hits = []
chimp_chr18_hits = []
chimp_chr19_hits = []
chimp_chr20_hits = []
chimp_chr21_hits = []
chimp_chr22_hits = []

vervet_chr1_hits = []
vervet_chr2_hits = []
vervet_chr3_hits = []
vervet_chr4_hits = []
vervet_chr5_hits = []
vervet_chr6_hits = []
vervet_chr7_hits = []
vervet_chr8_hits = []
vervet_chr9_hits = []
vervet_chr10_hits = []
vervet_chr11_hits = []
vervet_chr12_hits = []
vervet_chr13_hits = []
vervet_chr14_hits = []
vervet_chr15_hits = []
vervet_chr16_hits = []
vervet_chr17_hits = []
vervet_chr18_hits = []
vervet_chr19_hits = []
vervet_chr20_hits = []
vervet_chr21_hits = []
vervet_chr22_hits = []
vervet_chr23_hits = []
vervet_chr24_hits = []
vervet_chr25_hits = []
vervet_chr26_hits = []
vervet_chr27_hits = []
vervet_chr28_hits = []
vervet_chr29_hits = []

chimp_hit_indexes = []
vervet_hit_indexes = []

num_exons = len(human_exon_positions)
# num_exons = 50
for i in range(0, num_exons):
    susanne_pos_index = i
    
    vervet = vervet_exon_positions[i]
    vervet_colon = vervet.find(':') + 1
    vervet_hyphen = vervet.find('-')
    vervet_chr = vervet[: vervet_colon - 1]
    vervet_start = int(vervet[vervet_colon : vervet_hyphen])
    vervet_end = int(vervet[vervet_hyphen + 1 : -1])
    
    chimp = chimp_exon_positions[i]
    chimp_colon = chimp.find(':') + 1
    chimp_hyphen = chimp.find('-')
    chimp_chr = chimp[: chimp_colon - 1]
    chimp_start = int(chimp[chimp_colon : chimp_hyphen])
    chimp_end = int(chimp[chimp_hyphen + 1 : -1])

    if chimp_chr == 'chr1':
        find_orthologs('chimp', chimp_chr1_exon_info, chimp_chr1_hits, chimp_start, chimp_end, susanne_pos_index)
    elif chimp_chr == 'chr2a':
        find_orthologs('chimp', chimp_chr2a_exon_info, chimp_chr2a_hits, chimp_start, chimp_end, susanne_pos_index)
    elif chimp_chr == 'chr2b':
        find_orthologs('chimp', chimp_chr2b_exon_info, chimp_chr2b_hits, chimp_start, chimp_end, susanne_pos_index)
    elif chimp_chr == 'chr3':
        find_orthologs('chimp', chimp_chr3_exon_info, chimp_chr3_hits, chimp_start, chimp_end, susanne_pos_index)
    elif chimp_chr == 'chr4':
        find_orthologs('chimp', chimp_chr4_exon_info, chimp_chr4_hits, chimp_start, chimp_end, susanne_pos_index)
    elif chimp_chr == 'chr5':
        find_orthologs('chimp', chimp_chr5_exon_info, chimp_chr5_hits, chimp_start, chimp_end, susanne_pos_index)
    elif chimp_chr == 'chr6':
        find_orthologs('chimp', chimp_chr6_exon_info, chimp_chr6_hits, chimp_start, chimp_end, susanne_pos_index)
    elif chimp_chr == 'chr7':
        find_orthologs('chimp', chimp_chr7_exon_info, chimp_chr7_hits, chimp_start, chimp_end, susanne_pos_index)
    elif chimp_chr == 'chr8':
        find_orthologs('chimp', chimp_chr8_exon_info, chimp_chr8_hits, chimp_start, chimp_end, susanne_pos_index)
    elif chimp_chr == 'chr9':
        find_orthologs('chimp', chimp_chr9_exon_info, chimp_chr9_hits, chimp_start, chimp_end, susanne_pos_index)
    elif chimp_chr == 'chr10':
        find_orthologs('chimp', chimp_chr10_exon_info, chimp_chr10_hits, chimp_start, chimp_end, susanne_pos_index)
    elif chimp_chr == 'chr11':
        find_orthologs('chimp', chimp_chr11_exon_info, chimp_chr11_hits, chimp_start, chimp_end, susanne_pos_index)
    elif chimp_chr == 'chr12':
        find_orthologs('chimp', chimp_chr12_exon_info, chimp_chr12_hits, chimp_start, chimp_end, susanne_pos_index)
    elif chimp_chr == 'chr13':
        find_orthologs('chimp', chimp_chr13_exon_info, chimp_chr13_hits, chimp_start, chimp_end, susanne_pos_index)
    elif chimp_chr == 'chr14':
        find_orthologs('chimp', chimp_chr14_exon_info, chimp_chr14_hits, chimp_start, chimp_end, susanne_pos_index)
    elif chimp_chr == 'chr15':
        find_orthologs('chimp', chimp_chr15_exon_info, chimp_chr15_hits, chimp_start, chimp_end, susanne_pos_index)
    elif chimp_chr == 'chr16':
        find_orthologs('chimp', chimp_chr16_exon_info, chimp_chr16_hits, chimp_start, chimp_end, susanne_pos_index)
    elif chimp_chr == 'chr17':
        find_orthologs('chimp', chimp_chr17_exon_info, chimp_chr17_hits, chimp_start, chimp_end, susanne_pos_index)
    elif chimp_chr == 'chr18':
        find_orthologs('chimp', chimp_chr18_exon_info, chimp_chr18_hits, chimp_start, chimp_end, susanne_pos_index)
    elif chimp_chr == 'chr19':
        find_orthologs('chimp', chimp_chr19_exon_info, chimp_chr19_hits, chimp_start, chimp_end, susanne_pos_index)
    elif chimp_chr == 'chr20':
        find_orthologs('chimp', chimp_chr20_exon_info, chimp_chr20_hits, chimp_start, chimp_end, susanne_pos_index)
    elif chimp_chr == 'chr21':
        find_orthologs('chimp', chimp_chr21_exon_info, chimp_chr21_hits, chimp_start, chimp_end, susanne_pos_index)
    elif chimp_chr == 'chr22':
        find_orthologs('chimp', chimp_chr22_exon_info, chimp_chr22_hits, chimp_start, chimp_end, susanne_pos_index)
    
    
    
    if vervet_chr == 'chr1':
        find_orthologs('vervet', vervet_chr1_exon_info, vervet_chr1_hits, vervet_start, vervet_end, susanne_pos_index)
    elif vervet_chr == 'chr2':
        find_orthologs('vervet', vervet_chr2_exon_info, vervet_chr2_hits, vervet_start, vervet_end, susanne_pos_index)
    elif vervet_chr == 'chr3':
        find_orthologs('vervet', vervet_chr3_exon_info, vervet_chr3_hits, vervet_start, vervet_end, susanne_pos_index)
    elif vervet_chr == 'chr4':
        find_orthologs('vervet', vervet_chr4_exon_info, vervet_chr4_hits, vervet_start, vervet_end, susanne_pos_index)
    elif vervet_chr == 'chr5':
        find_orthologs('vervet', vervet_chr5_exon_info, vervet_chr5_hits, vervet_start, vervet_end, susanne_pos_index)
    elif vervet_chr == 'chr6':
        find_orthologs('vervet', vervet_chr6_exon_info, vervet_chr6_hits, vervet_start, vervet_end, susanne_pos_index)
    elif vervet_chr == 'chr7':
        find_orthologs('vervet', vervet_chr7_exon_info, vervet_chr7_hits, vervet_start, vervet_end, susanne_pos_index)
    elif vervet_chr == 'chr8':
        find_orthologs('vervet', vervet_chr8_exon_info, vervet_chr8_hits, vervet_start, vervet_end, susanne_pos_index)
    elif vervet_chr == 'chr9':
        find_orthologs('vervet', vervet_chr9_exon_info, vervet_chr9_hits, vervet_start, vervet_end, susanne_pos_index)
    elif vervet_chr == 'chr10':
        find_orthologs('vervet', vervet_chr10_exon_info, vervet_chr10_hits, vervet_start, vervet_end, susanne_pos_index)
    elif vervet_chr == 'chr11':
        find_orthologs('vervet', vervet_chr11_exon_info, vervet_chr11_hits, vervet_start, vervet_end, susanne_pos_index)
    elif vervet_chr == 'chr12':
        find_orthologs('vervet', vervet_chr12_exon_info, vervet_chr12_hits, vervet_start, vervet_end, susanne_pos_index)
    elif vervet_chr == 'chr13':
        find_orthologs('vervet', vervet_chr13_exon_info, vervet_chr13_hits, vervet_start, vervet_end, susanne_pos_index)
    elif vervet_chr == 'chr14':
        find_orthologs('vervet', vervet_chr14_exon_info, vervet_chr14_hits, vervet_start, vervet_end, susanne_pos_index)
    elif vervet_chr == 'chr15':
        find_orthologs('vervet', vervet_chr15_exon_info, vervet_chr15_hits, vervet_start, vervet_end, susanne_pos_index)
    elif vervet_chr == 'chr16':
        find_orthologs('vervet', vervet_chr16_exon_info, vervet_chr16_hits, vervet_start, vervet_end, susanne_pos_index)
    elif vervet_chr == 'chr17':
        find_orthologs('vervet', vervet_chr17_exon_info, vervet_chr17_hits, vervet_start, vervet_end, susanne_pos_index)
    elif vervet_chr == 'chr18':
        find_orthologs('vervet', vervet_chr18_exon_info, vervet_chr18_hits, vervet_start, vervet_end, susanne_pos_index)
    elif vervet_chr == 'chr19':
        find_orthologs('vervet', vervet_chr19_exon_info, vervet_chr19_hits, vervet_start, vervet_end, susanne_pos_index)
    elif vervet_chr == 'chr20':
        find_orthologs('vervet', vervet_chr20_exon_info, vervet_chr20_hits, vervet_start, vervet_end, susanne_pos_index)
    elif vervet_chr == 'chr21':
        find_orthologs('vervet', vervet_chr21_exon_info, vervet_chr21_hits, vervet_start, vervet_end, susanne_pos_index)
    elif vervet_chr == 'chr22':
        find_orthologs('vervet', vervet_chr22_exon_info, vervet_chr22_hits, vervet_start, vervet_end, susanne_pos_index)
    elif vervet_chr == 'chr23':
        find_orthologs('vervet', vervet_chr23_exon_info, vervet_chr23_hits, vervet_start, vervet_end, susanne_pos_index)
    elif vervet_chr == 'chr24':
        find_orthologs('vervet', vervet_chr24_exon_info, vervet_chr24_hits, vervet_start, vervet_end, susanne_pos_index)
    elif vervet_chr == 'chr25':
        find_orthologs('vervet', vervet_chr25_exon_info, vervet_chr25_hits, vervet_start, vervet_end, susanne_pos_index)
    elif vervet_chr == 'chr26':
        find_orthologs('vervet', vervet_chr26_exon_info, vervet_chr26_hits, vervet_start, vervet_end, susanne_pos_index)
    elif vervet_chr == 'chr27':
        find_orthologs('vervet', vervet_chr27_exon_info, vervet_chr27_hits, vervet_start, vervet_end, susanne_pos_index)
    elif vervet_chr == 'chr28':
        find_orthologs('vervet', vervet_chr28_exon_info, vervet_chr28_hits, vervet_start, vervet_end, susanne_pos_index)
    elif vervet_chr == 'chr29':
        find_orthologs('vervet', vervet_chr29_exon_info, vervet_chr29_hits, vervet_start, vervet_end, susanne_pos_index)

homologous_exons = list(set(chimp_hit_indexes).intersection(vervet_hit_indexes))
homologous_exons.sort()
number_of_homologous_exons = len(homologous_exons)

print('Chimp')
print(len(chimp_hit_indexes))
print(chimp_hit_indexes)

print('Vervet')
print(len(vervet_hit_indexes))
print(vervet_hit_indexes)

print('Overlaps')
print(homologous_exons)
print(number_of_homologous_exons)
print('\n')

chimp_table.close()
chimp_table_relevant.close()
chimp_table_formatted.close()
vervet_table.close()
vervet_table_relevant.close()
vervet_table_formatted.close()




#This section compiles all exons that had matches in both chimps and vervets, and are thus deemed homologous.
all_primates_file = 'C:/Users/Kellen/Downloads/primate_dfe_project/chimp_and_vervet_data/all_homo_exons_reference_exons_only.txt'
all_primates_table = open(all_primates_file, 'w+')
all_primates_table.write('#' + '\t' + 'human_gene' + '\t' + 'human_chr' + '\t' + 'human_strand' + '\t' + 'hg19_exon_start' + '\t' + 'hg19_exon_end' + '\t' + 'hg18_exon_start' + '\t' + 'hg18_exon_end' + '\t' + 'hg18_rec_start' + '\t' + 'hg18_rec_end' + '\t' + 'avg_rec_rate' + '\t' + 'divergence' + '\t' + 'chimp_transcript' + '\t' + 'chimp_gene_and_exon' + '\t' + 'chimp_chr' + '\t' + 'chimp_strand' + '\t' + 'chimp_exon_start' + '\t' + 'chimp_exon_end' + '\t' + 'vervet_transcript' + '\t' + 'vervet_gene_and_exon' + '\t' + 'vervet_chrom' + '\t' + 'vervet_strand' + '\t' + 'vervet_exon_start' + '\t' + 'vervet_exon_end' + '\n')


with open(chimp_file_formatted) as f:
    chimp_matched_exons = f.readlines()[1 :]
with open(vervet_file_formatted) as g:
    vervet_matched_exons = g.readlines()[1 :]

number_count = 1
for match in homologous_exons:
    human_line = '\t'.join(human_exon_positions[match].strip('\n').split()[1 : ])
    
#     human_colon=human.find(':')+1
#     human_hyphen=human.find('-')
#     human_chr=human[:human_colon-1]
#     human_start=int(human[human_colon:human_hyphen])
#     human_end=int(human[human_hyphen+1:-1])
#     human_name=find_human_name(human_chr)
#     human_line=human_chr+'\t'+human_name+'\t'+str(human_start)+'\t'+str(human_end)
#     human_line=human_chr+'\t'+str(human_start)+'\t'+str(human_end)

    for boop in chimp_matched_exons:
        if boop.split()[0] == str(match):
            split = boop.split()[1 :]
            join = '\t'.join(split)
            chimp_line = join
    for boop in vervet_matched_exons:
        if boop.split()[0] == str(match):
            split = boop.split()[1 :]
            join = '\t'.join(split)
            vervet_line = join
    output_line = str(number_count) + '\t' + human_line + '\t' + chimp_line + '\t' + vervet_line + '\n'
    all_primates_table.write(output_line)
    number_count += 1
    
all_primates_table.close()

