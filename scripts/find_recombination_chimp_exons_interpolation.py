


#This is a Python script that calculates recombination rate values for chimps using 
#linear interpolation with the rate map listings closest to the start and end of the given region.

#This script is currently designed to be run locally and will require path changes on the user's part to the variables:
#'rec_file' with input files 'chimp_Dec8_Haplotypes_Mar1_chr' + i + '-cleaned.txt'
#'exon_file' with input file "248_homo_exons_reference_with_intergenic_lengths.txt"
#first instance of 'result_file' with output file "chimp_homo_exons_recombination_interpolation.txt"
#'intergenic_file' with input file "chimp_all_intergenic_lengths.txt"
#second instance of 'result_file' with output file "chimp_all_recombination_interpolation.txt"




from bisect import bisect_left
    




#This section reads in the chimp recombination rate map for all chromosomes at once and makes a dictionary 
#for each chromosome containing a list of every listed map point and the rate at that point.

all_chr = ['1', '2a', '2b', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22']
chr_positions = {}
for i in all_chr:
    chr_num = 'chr' + i
    chr_positions[chr_num] = []

position_rates = {}
for i in all_chr:
    chr_num = 'chr' + i
    rec_file = 'C:/Users/Kellen/Downloads/primate_dfe_project/chimp_and_vervet_data/chimp_recombination/chimp_Dec8_Haplotypes_Mar1_chr' + i + '-cleaned.txt'

    with open(rec_file, 'r') as g:
        rec_data = g.readlines()[1 :]
    
    for j in rec_data:
        split = j.strip('\n').split('\t')
        position = int(split[0])
        rate_4nerkb = float(split[1])
        rate_4ner = float(split[2])
        
        chr_positions[chr_num].append(position)
        
        position_rates_str = chr_num + '_' + str(position)
        position_rates[position_rates_str] = [rate_4nerkb, rate_4ner]
        
for i in all_chr:
    chr_num = 'chr' + i
    final_pos = chr_positions[chr_num][-1]
    final_pos_str = chr_num + '_' + str(final_pos)
    final_rate = position_rates[final_pos_str]
    chr_positions[chr_num].append(999999999999999999999999999999999999)
    position_rates_str = chr_num + '_' + str(999999999999999999999999999999999999)
    position_rates[position_rates_str] = final_rate
    
    chr_positions[chr_num].sort()
    




#This function takes in the start and end position of the region, that being the 5' start/end and 
#the 3' start/end depending on the strand. The function then looks in the list for that chromosome and 
#uses bisect to find the given map point closest to the starts and ends. It then verifies it's the closest
#and proceeds to calculate the rate using linear interpolation with the two closest points and the given 
#4Ner rate values at those positions. Finally, it checks to see if the region has any inversions based on 
#whether or not there is ever a decrease in the rate. I wrote a basic script to check if any inversions 
#existed at all, and I don't find any, but I left this part in there anyway just in case.

def find_rec_rates(true_start, true_end, chr_num):
    inversions = []

    initial_closest_start_pos = bisect_left(chr_positions[chr_num], true_start)
    start_given = int(chr_positions[chr_num][initial_closest_start_pos])
    start_previous = int(chr_positions[chr_num][initial_closest_start_pos - 1])
    check_start_given = abs(true_start - start_given)
    check_start_previous = abs(true_start - start_previous)
    if check_start_given > check_start_previous:
        start_pos = start_previous
        initial_closest_start_pos -= 1
    elif check_start_given < check_start_previous:
        start_pos = start_given
    elif check_start_given == check_start_previous:
        start_pos = start_given

    initial_closest_end_pos = bisect_left(chr_positions[chr_num], true_end)
    end_given = int(chr_positions[chr_num][initial_closest_end_pos])
    end_previous = int(chr_positions[chr_num][initial_closest_end_pos - 1])
    check_end_given = abs(true_end - end_given)
    check_end_previous = abs(true_end - end_previous)
    if check_end_given > check_end_previous:
        end_pos = end_previous
        initial_closest_end_pos -= 1
    elif check_end_given < check_end_previous:
        end_pos = end_given
    elif check_end_given == check_end_previous:
        end_pos = end_given


    if (start_pos == end_pos):
        position_rates_str = chr_num + '_' + str(start_pos)
        rec_rate = position_rates[position_rates_str][0]
        rec_rate = float(round(rec_rate * 10**9) / 10**9)
    else:
        build_start_str = chr_num + '_' + str(start_pos)
        start_rate = position_rates[build_start_str][1]
        build_end_str = chr_num + '_' + str(end_pos)
        end_rate = position_rates[build_end_str][1]
        y_dif = end_rate - start_rate
        x_dif = end_pos - start_pos
        rec_rate_per_bp = y_dif / x_dif
        rec_rate_per_kb = rec_rate_per_bp * 1000
        rec_rate = float(round(rec_rate_per_kb * 10**9) / 10**9) 

    possible_inversion_pos_and_rates = []
    if initial_closest_start_pos != 0:
        for index in range(initial_closest_start_pos - 1, initial_closest_end_pos + 1):
            pos_at_index = int(chr_positions[chr_num][index])
            build_pos_at_index_str = chr_num + '_' + str(pos_at_index)
            rate_at_pos = position_rates[build_pos_at_index_str][1]
            possible_inversion_pos_and_rates.append([pos_at_index, rate_at_pos])
        list_len = len(possible_inversion_pos_and_rates)
        for i in range(1, list_len):
            current_rate = possible_inversion_pos_and_rates[i][1]
            previous_rate = possible_inversion_pos_and_rates[i - 1][1]
            if previous_rate > current_rate:
                current_pos = possible_inversion_pos_and_rates[i][0]
                inversions.append(str(current_pos))
    if len(inversions) == 0:
        inversions.append('NA')
    
    info = [start_pos, end_pos, rec_rate, inversions]
    return info




#This section runs the script for just the 248 homologous chimp exons

exon_file = 'C:/Users/Kellen/Downloads/primate_dfe_project/chimp_and_vervet_data/248_homo_exons_reference_with_intergenic_lengths.txt'
with open(exon_file, 'r') as g:
    exon_data = g.readlines()[1 :]
    
result_file = 'C:/Users/Kellen/Downloads/primate_dfe_project/chimp_and_vervet_data/chimp_homo_exons_recombination_interpolation.txt'
results = open(result_file, 'w+')
results.write('name' + '\t' + 'exon' + '\t' + 'chr' + '\t' + 'strand' + '\t' + 'closest_rec_start_pos' + '\t' + 'closest_rec_end_pos' + '\t' + '4Ner/kp_rec_rate' + '\t' + 'inversion_sites' + '\n')

total_num_lines = len(exon_data)
for w in range(0, total_num_lines):
    line = exon_data[w].strip('\n')
    split = line.split()
    name = split[16]
    exon = split[17]
    chr_num = split[18]
    strand = split[19]
    
    exon_start = int(split[22])
    exon_end = int(split[23])
    start_5p = int(split[25])
    end_5p = int(split[26])
    start_3p = int(split[28])
    end_3p = int(split[29])
    
    if strand == '+':
        info = find_rec_rates(start_5p, end_3p, chr_num)
    if strand == '-':
        info = find_rec_rates(start_3p, end_5p, chr_num)
    rec_start_pos = info[0]
    rec_end_pos = info[1]
    final_rate = info[2]
    inversions = info[3]
    inversions_str = ', '.join(inversions)
    
    results.write(name + '\t' + exon + '\t' + chr_num + '\t' + strand + '\t' + str(rec_start_pos) + '\t' + str(rec_end_pos) + '\t' + str(final_rate) + '\t' + str(inversions_str)+'\n')    
results.close()




#This section runs the script for just all possible regions based on the chimp UCSC annotation.

intergenic_file = 'C:/Users/Kellen/Downloads/primate_dfe_project/chimp_and_vervet_data/chimp_all_intergenic_lengths.txt'
with open(intergenic_file, 'r') as g:
    intergenic_data = g.readlines()[1 :]
    
result_file = 'C:/Users/Kellen/Downloads/primate_dfe_project/chimp_and_vervet_data/chimp_all_recombination_interpolation.txt'
results = open(result_file, 'w+')
results.write('name' + '\t' + 'exon' + '\t' + 'chr' + '\t' + 'strand' + '\t' + 'closest_rec_start_pos' + '\t' + 'closest_rec_end_pos' + '\t' + '4Ner/kp_rec_rate' + '\t' + 'inversion_sites' + '\n')

total_num_lines = len(intergenic_data)
for w in range(0, total_num_lines):
    line = intergenic_data[w].strip('\n')
    split = line.split()
    name = split[1]
    exon = split[2]
    chr_num = split[3]
    strand = split[4]
    
    exon_start = split[7]
    exon_end = split[8]
    start_5p = split[10]
    end_5p = split[11]
    start_3p = split[13]
    end_3p = split[14]
    
    if strand == '+':
        if start_5p == 'NA' or end_3p == 'NA':
            info = ['NA','NA','NA','NA']
        else:
            info = find_rec_rates(int(start_5p), int(end_3p), chr_num)
    if strand == '-':
        if start_3p == 'NA' or end_5p == 'NA':
            info = ['NA','NA','NA','NA']
        else:
            info = find_rec_rates(int(start_3p), int(end_5p), chr_num)
    rec_start_pos = info[0]
    rec_end_pos = info[1]
    final_rate = info[2]
    inversions = info[3]
    inversions_str = ', '.join(inversions)
    if inversions_str == 'N, A':
        inversions_str = 'NA'
    
    results.write(name + '\t' + exon + '\t' + chr_num + '\t' + strand + '\t' + str(rec_start_pos) + '\t' + str(rec_end_pos) + '\t' + str(final_rate) + '\t' + str(inversions_str) + '\n')    
results.close()

