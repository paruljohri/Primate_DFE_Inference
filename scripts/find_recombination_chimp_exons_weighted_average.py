


#This is a Python script that calculates recombination rate values for chimp regions by taking 
#weighted averages across the various rate windows in the rate map that a region spans.

#This script is currently designed to be run locally and will require path changes on the user's part to the variables:
#'rec_file' with input files 'chimp_Dec8_Haplotypes_Mar1_chr' + i + '-cleaned.txt'
#'exon_file' with input file "248_homo_exons_reference_with_intergenic_lengths.txt"
#first instance of 'result_file' with output file "chimp_homo_exons_recombination_weighted_average.txt"
#'intergenic_file' with input file "chimp_all_intergenic_lengths.txt"
#second instance of 'result_file' with output file "chimp_all_exons_recombination_weighted_average.txt"




from bisect import bisect_left




#This section reads in the recombination rate map for each chimp chromosome into a list.
all_chr = ['1', '2a', '2b', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22']
chr1_rec = []
chr2a_rec = []
chr2b_rec = []
chr3_rec = []
chr4_rec = []
chr5_rec = []
chr6_rec = []
chr7_rec = []
chr8_rec = []
chr9_rec = []
chr10_rec = []
chr11_rec = []
chr12_rec = []
chr13_rec = []
chr14_rec = []
chr15_rec = []
chr16_rec = []
chr17_rec = []
chr18_rec = []
chr19_rec = []
chr20_rec = []
chr21_rec = []
chr22_rec = []

chr1_rec_pos = []
chr2a_rec_pos = []
chr2b_rec_pos = []
chr3_rec_pos = []
chr4_rec_pos = []
chr5_rec_pos = []
chr6_rec_pos = []
chr7_rec_pos = []
chr8_rec_pos = []
chr9_rec_pos = []
chr10_rec_pos = []
chr11_rec_pos = []
chr12_rec_pos = []
chr13_rec_pos = []
chr14_rec_pos = []
chr15_rec_pos = []
chr16_rec_pos = []
chr17_rec_pos = []
chr18_rec_pos = []
chr19_rec_pos = []
chr20_rec_pos = []
chr21_rec_pos = []
chr22_rec_pos = []

# Position(bp)	Rate(4Ner/kb)	Map(4Ner)	Filtered
for i in all_chr:
    chr_num = 'chr' + i
#     rec_file = 'C:/Users/kelle/Downloads/primate_stuff_test/chimp_data/recombination/chimp_Dec8_Haplotypes_Mar1_chr'+i+'-cleaned.txt'
    rec_file = 'C:/Users/Kellen/Downloads/primate_dfe_project/chimp_and_vervet_data/chimp_recombination/chimp_Dec8_Haplotypes_Mar1_chr' + i + '-cleaned.txt'

    with open(rec_file, 'r') as g:
        rec_data = g.readlines()[1 : ]
    
    for j in rec_data:
        split = j.strip('\n').split('\t')
        position = int(split[0])
        rate_per_kb = float(split[1])
        rate = float(split[2])
        rec_info = [position, rate_per_kb, rate]
    
        if chr_num == 'chr1':
            chr1_rec.append(rec_info)
            chr1_rec_pos.append(position)
        elif chr_num == 'chr2a':
            chr2a_rec.append(rec_info)
            chr2a_rec_pos.append(position)
        elif chr_num == 'chr2b':
            chr2b_rec.append(rec_info)
            chr2b_rec_pos.append(position)
        elif chr_num == 'chr3':
            chr3_rec.append(rec_info)
            chr3_rec_pos.append(position)
        elif chr_num == 'chr4':
            chr4_rec.append(rec_info)
            chr4_rec_pos.append(position)
        elif chr_num == 'chr5':
            chr5_rec.append(rec_info)
            chr5_rec_pos.append(position)
        elif chr_num == 'chr6':
            chr6_rec.append(rec_info)
            chr6_rec_pos.append(position)
        elif chr_num == 'chr7':
            chr7_rec.append(rec_info)
            chr7_rec_pos.append(position)
        elif chr_num == 'chr8':
            chr8_rec.append(rec_info)
            chr8_rec_pos.append(position)
        elif chr_num == 'chr9':
            chr9_rec.append(rec_info)
            chr9_rec_pos.append(position)
        elif chr_num == 'chr10':
            chr10_rec.append(rec_info)
            chr10_rec_pos.append(position)
        elif chr_num == 'chr11':
            chr11_rec.append(rec_info)
            chr11_rec_pos.append(position)
        elif chr_num == 'chr12':
            chr12_rec.append(rec_info)
            chr12_rec_pos.append(position)
        elif chr_num == 'chr13':
            chr13_rec.append(rec_info)
            chr13_rec_pos.append(position)
        elif chr_num == 'chr14':
            chr14_rec.append(rec_info)
            chr14_rec_pos.append(position)
        elif chr_num == 'chr15':
            chr15_rec.append(rec_info)
            chr15_rec_pos.append(position)
        elif chr_num == 'chr16':
            chr16_rec.append(rec_info)
            chr16_rec_pos.append(position)
        elif chr_num == 'chr17':
            chr17_rec.append(rec_info)
            chr17_rec_pos.append(position)
        elif chr_num == 'chr18':
            chr18_rec.append(rec_info)
            chr18_rec_pos.append(position)
        elif chr_num == 'chr19':
            chr19_rec.append(rec_info)
            chr19_rec_pos.append(position)
        elif chr_num == 'chr20':
            chr20_rec.append(rec_info)
            chr20_rec_pos.append(position)
        elif chr_num == 'chr21':
            chr21_rec.append(rec_info)
            chr21_rec_pos.append(position)
        elif chr_num == 'chr22':
            chr22_rec.append(rec_info)
            chr22_rec_pos.append(position)

chr1_rec.append([99999999999999999999, 0, 0])
chr2a_rec.append([99999999999999999999, 0, 0])
chr2b_rec.append([99999999999999999999, 0, 0])
chr3_rec.append([99999999999999999999, 0, 0])
chr4_rec.append([99999999999999999999, 0, 0])
chr5_rec.append([99999999999999999999, 0, 0])
chr6_rec.append([99999999999999999999, 0, 0])
chr7_rec.append([99999999999999999999, 0, 0])
chr8_rec.append([99999999999999999999, 0, 0])
chr9_rec.append([99999999999999999999, 0, 0])
chr10_rec.append([99999999999999999999, 0, 0])
chr11_rec.append([99999999999999999999, 0, 0])
chr12_rec.append([99999999999999999999, 0, 0])
chr13_rec.append([99999999999999999999, 0, 0])
chr14_rec.append([99999999999999999999, 0, 0])
chr15_rec.append([99999999999999999999, 0, 0])
chr16_rec.append([99999999999999999999, 0, 0])
chr17_rec.append([99999999999999999999, 0, 0])
chr18_rec.append([99999999999999999999, 0, 0])
chr19_rec.append([99999999999999999999, 0, 0])
chr20_rec.append([99999999999999999999, 0, 0])
chr21_rec.append([99999999999999999999, 0, 0])
chr22_rec.append([99999999999999999999, 0, 0])

chr1_rec_pos.append(99999999999999999999)
chr2a_rec_pos.append(99999999999999999999)
chr2b_rec_pos.append(99999999999999999999)
chr3_rec_pos.append(99999999999999999999)
chr4_rec_pos.append(99999999999999999999)
chr5_rec_pos.append(99999999999999999999)
chr6_rec_pos.append(99999999999999999999)
chr7_rec_pos.append(99999999999999999999)
chr8_rec_pos.append(99999999999999999999)
chr9_rec_pos.append(99999999999999999999)
chr10_rec_pos.append(99999999999999999999)
chr11_rec_pos.append(99999999999999999999)
chr12_rec_pos.append(99999999999999999999)
chr13_rec_pos.append(99999999999999999999)
chr14_rec_pos.append(99999999999999999999)
chr15_rec_pos.append(99999999999999999999)
chr16_rec_pos.append(99999999999999999999)
chr17_rec_pos.append(99999999999999999999)
chr18_rec_pos.append(99999999999999999999)
chr19_rec_pos.append(99999999999999999999)
chr20_rec_pos.append(99999999999999999999)
chr21_rec_pos.append(99999999999999999999)
chr22_rec_pos.append(99999999999999999999)




#This function takes in the start and end position of the region, that being the 5' start/end and the
#3' start/end depending on the strand. The function then looks in the list for that chromosome and uses 
#bisect to find the given map point closest to the starts and ends. It then verifies it is the closest 
#and proceeds to calculate the rate using weighted averages. This is done by taking the rate for a window
#that an exon overlaps with as the product of the rate with the percentage of the exon that is in the window.
#These weighted rates are then added together to obtain the total rate. For example, imagine an exon spans 
#the region from 1 to 100 and overlaps with the windows 1-20, 21-70, and 71-120. These regions have 
#recombination rates of 1, 2, and 3, respectively. The final rate would then be 0.2 * 1 + 0.5 * 2 + 0.3 *3 = 2.1.

def find_rec_rates(chrXXX_rec, chrXXX_rec_pos, exon_start, exon_end):
    closest_position = bisect_left(chrXXX_rec_pos, exon_start)
#     picked_value=chrXXX_rec_pos[closest_position]
#     print('\n')
#     print(exon_start)
#     print(exon_end)
#     print(closest_position)
#     print(picked_value)
    relevant_pos = []
    relevant_pos.append(list(chrXXX_rec[closest_position - 1]))
    relevant_pos.append(list(chrXXX_rec[closest_position]))
    final_position = closest_position
    for j in range(closest_position + 1, len(chrXXX_rec)):
        if chrXXX_rec[j][0] < exon_end:
            info = list(chrXXX_rec[j])
            relevant_pos.append(info)
            final_position = j
    if final_position + 1 != len(chrXXX_rec_pos):  
        relevant_pos.append(list(chrXXX_rec[final_position + 1]))

    num_pos = len(relevant_pos)
    relevant_pos[0].insert(0, 'X')
    for k in range(1,num_pos): 
        relevant_pos[k].insert(0, relevant_pos[k - 1][1])
    del relevant_pos[0]
#     print(relevant_pos)
    
    rate_per_kb_total = 0
    rate_total = 0
    exon_length = exon_end - exon_start
    all_in_one = False
    for window in relevant_pos:
        window_start = window[0]
        window_end = window[1]
    #     print(window_start)
    #     print(window_end)
        window_length = window_end-window_start
        rate_per_kb = window[2]
        rate = window[3]
        #If the whole exon is entirely within rec window
        if (window_start <= exon_start < window_end) and (window_start < exon_end <= window_end):
            final_rate_per_kb = rate_per_kb
            final_rate = rate
            all_in_one = True
#             print('1')

        #If part of the exon spans a whole rec window
        elif (exon_start <= window_start) and (window_end <= exon_end):
            rate_per_kb_total += rate_per_kb * window_length
            rate_total += rate * window_length
#             print('2')

        #If exon window beginning is partially in one window
        elif (window_start < exon_start < window_end) and (window_end < exon_end):
            rate_per_kb_total += rate_per_kb * (window_end - exon_start)
            rate_total += rate * (window_end - exon_start)
#             print('3')

        #If exon window ending is partially in one window
        elif exon_start < window_start and window_start < exon_end < window_end:
            rate_per_kb_total += rate_per_kb * (exon_end - window_start)
            rate_total += rate * (window_end-exon_start)
#             print('4')

    if all_in_one == False:
        final_rate_per_kb = rate_per_kb_total / exon_length
        final_rate = rate_total / exon_length
    
    return final_rate_per_kb, final_rate




#This section caluclates the recombination rate for the 248 homologous chimp regions.

"""For HOMOLOGOUS exons"""
exon_file = 'C:/Users/Kellen/Downloads/primate_dfe_project/chimp_and_vervet_data/248_homo_exons_reference_with_intergenic_lengths.txt'
with open(exon_file, 'r') as g:
    exon_data = g.readlines()[1 : ]
    
result_file = 'C:/Users/Kellen/Downloads/primate_dfe_project/chimp_and_vervet_data/chimp_homo_exons_recombination_weighted_average.txt'
results = open(result_file, 'w+')
results.write('name' + '\t' + 'exon' + '\t' + 'chr' + '\t' + 'strand' + '\t' + 'start' + '\t' + 'end' + '\t' + '4Ner/kb' + '\t' + '4Ner' + '\n')

total_num_lines = len(exon_data)
for w in range(0, total_num_lines):
    line = exon_data[w].strip('\n')
    split = line.split()
    name = split[16]
    exon = split[17]
    chr_num = split[18]
    strand = split[19]

    start_5p = int(split[25])
    end_5p = int(split[26])
    start_3p = int(split[28])
    end_3p = int(split[29])
    
    if strand == '+':
        exon_start = start_5p
        exon_end = end_3p
    if strand == '-':
        exon_start = start_3p
        exon_end = end_5p
    
    if chr_num == 'chr1':
        final_rate_per_kb, final_rate = find_rec_rates(chr1_rec, chr1_rec_pos, exon_start, exon_end)
    elif chr_num == 'chr2a':
        final_rate_per_kb, final_rate = find_rec_rates(chr2a_rec, chr2a_rec_pos, exon_start, exon_end)
    elif chr_num == 'chr2b':
        final_rate_per_kb, final_rate = find_rec_rates(chr2b_rec, chr2b_rec_pos, exon_start, exon_end)
    elif chr_num == 'chr3':
        final_rate_per_kb, final_rate = find_rec_rates(chr3_rec, chr3_rec_pos, exon_start, exon_end)
    elif chr_num == 'chr4':
        final_rate_per_kb, final_rate = find_rec_rates(chr4_rec, chr4_rec_pos, exon_start, exon_end)
    elif chr_num == 'chr5':
        final_rate_per_kb, final_rate = find_rec_rates(chr5_rec, chr5_rec_pos, exon_start, exon_end)
    elif chr_num == 'chr6':
        final_rate_per_kb, final_rate = find_rec_rates(chr6_rec, chr6_rec_pos, exon_start, exon_end)
    elif chr_num == 'chr7':
        final_rate_per_kb, final_rate = find_rec_rates(chr7_rec, chr7_rec_pos, exon_start, exon_end)
    elif chr_num == 'chr8':
        final_rate_per_kb, final_rate = find_rec_rates(chr8_rec, chr8_rec_pos, exon_start, exon_end)
    elif chr_num == 'chr9':
        final_rate_per_kb, final_rate = find_rec_rates(chr9_rec, chr9_rec_pos, exon_start, exon_end)
    elif chr_num == 'chr10':
        final_rate_per_kb, final_rate = find_rec_rates(chr10_rec, chr10_rec_pos, exon_start, exon_end)
    elif chr_num == 'chr11':
        final_rate_per_kb, final_rate = find_rec_rates(chr11_rec, chr11_rec_pos, exon_start, exon_end)
    elif chr_num == 'chr12':
        final_rate_per_kb, final_rate = find_rec_rates(chr12_rec, chr12_rec_pos, exon_start, exon_end)
    elif chr_num == 'chr13':
        final_rate_per_kb, final_rate = find_rec_rates(chr13_rec, chr13_rec_pos, exon_start, exon_end)
    elif chr_num == 'chr14':
        final_rate_per_kb, final_rate = find_rec_rates(chr14_rec, chr14_rec_pos, exon_start, exon_end)
    elif chr_num == 'chr15':
        final_rate_per_kb, final_rate = find_rec_rates(chr15_rec, chr15_rec_pos, exon_start, exon_end)
    elif chr_num == 'chr16':
        final_rate_per_kb, final_rate = find_rec_rates(chr16_rec, chr16_rec_pos, exon_start, exon_end)
    elif chr_num == 'chr17':
        final_rate_per_kb, final_rate = find_rec_rates(chr17_rec, chr17_rec_pos, exon_start, exon_end)
    elif chr_num == 'chr18':
        final_rate_per_kb, final_rate = find_rec_rates(chr18_rec, chr18_rec_pos, exon_start, exon_end)
    elif chr_num == 'chr19':
        final_rate_per_kb, final_rate = find_rec_rates(chr19_rec, chr19_rec_pos, exon_start, exon_end)
    elif chr_num == 'chr20':
        final_rate_per_kb, final_rate = find_rec_rates(chr20_rec, chr20_rec_pos, exon_start, exon_end)
    elif chr_num == 'chr21':
        final_rate_per_kb, final_rate = find_rec_rates(chr21_rec, chr21_rec_pos, exon_start, exon_end)
    elif chr_num == 'chr22':
        final_rate_per_kb, final_rate = find_rec_rates(chr22_rec, chr22_rec_pos, exon_start, exon_end)
    else:
        continue
        
    results.write(name + '\t' + exon + '\t' + chr_num + '\t' + strand + '\t' + str(exon_start) + '\t' + str(exon_end) + '\t' + str(final_rate_per_kb) + '\t' + str(final_rate) + '\n')
results.close()




#This section caluclates the recombination rate for all regions based on the the intergenic 
#lengths calculated from the chimp UCSC exon annotation.

"""For ALL exons"""
intergenic_file = 'C:/Users/Kellen/Downloads/primate_dfe_project/chimp_and_vervet_data/chimp_all_intergenic_lengths.txt'
with open(intergenic_file, 'r') as g:
    intergenic_data = g.readlines()[1 : ]
    
result_file = 'C:/Users/Kellen/Downloads/primate_dfe_project/chimp_and_vervet_data/chimp_all_exons_recombination_weighted_average.txt'
results = open(result_file, 'w+')
results.write('name' + '\t' + 'exon' + '\t' + 'chr' + '\t' + 'strand' + '\t' + 'start' + '\t' + 'end' + '\t' + '4Ner/kb' + '\t' + '4Ner' + '\n')

    
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
    start_3p = split[13]
    end_3p = split[14]
    
    if strand == '+':
        exon_start = start_5p
        exon_end = end_3p
    if strand == '-':
        exon_start = start_3p
        exon_end = end_5p
        
    if exon_start == 'NA' or exon_end == 'NA':
        final_rate_per_kb = 0
        final_rate = 0
    else:
        exon_start = int(exon_start)
        exon_end = int(exon_end)
        
        if chr_num == 'chr1':
            final_rate_per_kb, final_rate = find_rec_rates(chr1_rec, chr1_rec_pos, exon_start, exon_end)
        elif chr_num == 'chr2a':
            final_rate_per_kb, final_rate = find_rec_rates(chr2a_rec, chr2a_rec_pos, exon_start, exon_end)
        elif chr_num == 'chr2b':
            final_rate_per_kb, final_rate = find_rec_rates(chr2b_rec, chr2b_rec_pos, exon_start, exon_end)
        elif chr_num == 'chr3':
            final_rate_per_kb, final_rate = find_rec_rates(chr3_rec, chr3_rec_pos, exon_start, exon_end)
        elif chr_num == 'chr4':
            final_rate_per_kb, final_rate = find_rec_rates(chr4_rec, chr4_rec_pos, exon_start, exon_end)
        elif chr_num == 'chr5':
            final_rate_per_kb, final_rate = find_rec_rates(chr5_rec, chr5_rec_pos, exon_start, exon_end)
        elif chr_num == 'chr6':
            final_rate_per_kb, final_rate = find_rec_rates(chr6_rec, chr6_rec_pos, exon_start, exon_end)
        elif chr_num == 'chr7':
            final_rate_per_kb, final_rate = find_rec_rates(chr7_rec, chr7_rec_pos, exon_start, exon_end)
        elif chr_num == 'chr8':
            final_rate_per_kb, final_rate = find_rec_rates(chr8_rec, chr8_rec_pos, exon_start, exon_end)
        elif chr_num == 'chr9':
            final_rate_per_kb, final_rate = find_rec_rates(chr9_rec, chr9_rec_pos, exon_start, exon_end)
        elif chr_num == 'chr10':
            final_rate_per_kb, final_rate = find_rec_rates(chr10_rec, chr10_rec_pos, exon_start, exon_end)
        elif chr_num == 'chr11':
            final_rate_per_kb, final_rate = find_rec_rates(chr11_rec, chr11_rec_pos, exon_start, exon_end)
        elif chr_num == 'chr12':
            final_rate_per_kb, final_rate = find_rec_rates(chr12_rec, chr12_rec_pos, exon_start, exon_end)
        elif chr_num == 'chr13':
            final_rate_per_kb, final_rate = find_rec_rates(chr13_rec, chr13_rec_pos, exon_start, exon_end)
        elif chr_num == 'chr14':
            final_rate_per_kb, final_rate = find_rec_rates(chr14_rec, chr14_rec_pos, exon_start, exon_end)
        elif chr_num == 'chr15':
            final_rate_per_kb, final_rate = find_rec_rates(chr15_rec, chr15_rec_pos, exon_start, exon_end)
        elif chr_num == 'chr16':
            final_rate_per_kb, final_rate = find_rec_rates(chr16_rec, chr16_rec_pos, exon_start, exon_end)
        elif chr_num == 'chr17':
            final_rate_per_kb, final_rate = find_rec_rates(chr17_rec, chr17_rec_pos, exon_start, exon_end)
        elif chr_num == 'chr18':
            final_rate_per_kb, final_rate = find_rec_rates(chr18_rec, chr18_rec_pos, exon_start, exon_end)
        elif chr_num == 'chr19':
            final_rate_per_kb, final_rate = find_rec_rates(chr19_rec, chr19_rec_pos, exon_start, exon_end)
        elif chr_num == 'chr20':
            final_rate_per_kb, final_rate = find_rec_rates(chr20_rec, chr20_rec_pos, exon_start, exon_end)
        elif chr_num == 'chr21':
            final_rate_per_kb, final_rate = find_rec_rates(chr21_rec, chr21_rec_pos, exon_start, exon_end)
        elif chr_num == 'chr22':
            final_rate_per_kb, final_rate = find_rec_rates(chr22_rec, chr22_rec_pos, exon_start, exon_end)
        else:
            continue

    results.write(name + '\t' + exon + '\t' + chr_num + '\t' + strand + '\t' + str(exon_start) + '\t' + str(exon_end) + '\t' + str(final_rate_per_kb) + '\t' + str(final_rate) + '\n')
results.close()

