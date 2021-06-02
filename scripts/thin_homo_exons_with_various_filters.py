


#This is a Python script that works to thin out the 248 homologous exons by applying various filters, leaving us with a new 
#outputed list of exons.

#This script is currently designed to be run locally and will require path changes on the user's part to the variables:
#'homo_exon_file' with input file "248_homo_exons_reference_all_info.txt"
#'alignment_check_file' with input file "248_exon_alignment_check.txt"
#'chimp_ensembl_file' with input file "ensembl_biomart_homo_exons_chimp.txt"
#'vervet_ensembl_file' with input file "ensembl_biomart_homo_exons_vervet.txt"
#'new_filtered_homo_exon_file' with output file "FILTERED_homo_exons_reference_all_info.txt"




#This is a list that has the lines with exons (each corresponding to one of the 248 exons) that 
#fail to pass the filter (for EITHER the human, chimp, or vervet exon) added to it. At the end, 
#any lines NOT in it are printed to the final set of exons.

avoid_lines = []




#This section reads in all relevant information about the 248 homologous exons.

homo_exon_file = 'C:/Users/Kellen/Downloads/primate_dfe_project/chimp_and_vervet_data/248_homo_exons_reference_all_info.txt'
with open(homo_exon_file, 'r') as f:
    homo_exon_data = f.readlines()
homo_exon_data_column_names = homo_exon_data[0]
del homo_exon_data[0]




#This section applies a series of small filters detailed further below above each if statement.

counter = 0
for line in homo_exon_data:
    avoid = False
    split = line.split()
    
    human_5p_start = split[9]
    human_5p_end = split[10]
    human_5p_len = split[11]
    human_3p_start = split[12]
    human_3p_end = split[13]
    human_3p_len = split[14]
    human_exon_div = split[15]
    human_5p_div = split[16]
    human_3p_div = split[17]
    
    chimp_useability = split[21]
    chimp_5p_start = split[31]
    chimp_5p_end = split[32]
    chimp_5p_len = split[33]
    chimp_3p_start = split[34]
    chimp_3p_end = split[35]
    chimp_3p_len = split[36]
    chimp_exon_div = split[37]
    chimp_5p_div = split[39]
    chimp_3p_div = split[40]
    chimp_inversions = split[47]
    
    vervet_useability = split[48]
    vervet_5p_start = split[58]
    vervet_5p_end = split[59]
    vervet_5p_len = split[60]
    vervet_3p_start = split[61]
    vervet_3p_end = split[62]
    vervet_3p_len = split[63]
    vervet_exon_div = split[64]
    vervet_5p_div = split[66]
    vervet_3p_div = split[67]
    vervet_inversions = split[74]

    relevant_values = [human_5p_start, human_5p_end, human_5p_len, human_3p_start, human_3p_end, human_3p_len, human_exon_div, human_5p_div, human_3p_div, chimp_5p_start, chimp_5p_end, chimp_5p_len, chimp_3p_start, chimp_3p_end, chimp_3p_len, chimp_exon_div, chimp_5p_div, chimp_3p_div, vervet_5p_start, vervet_5p_end, vervet_5p_len, vervet_3p_start, vervet_3p_end, vervet_3p_len, vervet_exon_div, vervet_5p_div, vervet_3p_div]
    #This eleminates any exon group where one at least on of the species has a missing exon length, integrenic length, or divergence value noted by an "NA" in the spot.
    if 'NA' in relevant_values:
        avoid = True
    #This eliminates any exon group where either the vervet or chimp exon there had duplicates in the UCSC annotation
    if chimp_useability == 'N' or vervet_useability == 'N':
        avoid = True
    #This eliminates any exon group where either the vervet or chimp exon there had an inversion according to the recombination map, though none do currently. 
    if chimp_inversions != 'NA' or vervet_inversions != 'NA':
        avoid = True
        
    if avoid == True:
        avoid_lines.append(counter)
    counter += 1
    
print(avoid_lines)




#This section filters based on intergenic lengths of the chimp and vervet, only taking those where the 
#intergenic region is at least X% of the corresponding human intergenic length (matching 5' with 5' and 3' with 3')

acceptable_length_percentage = 0.75

counter = 0
for line in homo_exon_data:
    avoid = False
    split = line.split()
    human_5p_len = split[11]
    human_3p_len = split[14]
    chimp_5p_len = split[33]
    chimp_3p_len = split[36]
    vervet_5p_len = split[60]
    vervet_3p_len = split[63]
    
    if human_5p_len == 'NA' or human_3p_len == 'NA' or chimp_5p_len == 'NA' or chimp_3p_len == 'NA' or vervet_5p_len == 'NA' or vervet_3p_len == 'NA':
        avoid_lines.append(counter)
        counter += 1
        continue
    
    human_5p_len = float(human_5p_len) * acceptable_length_percentage
    human_3p_len = float(human_3p_len) * acceptable_length_percentage
    chimp_5p_len = float(chimp_5p_len)
    chimp_3p_len = float(chimp_3p_len)
    vervet_5p_len = float(vervet_5p_len)
    vervet_3p_len = float(vervet_3p_len)
    
    if chimp_5p_len >= human_5p_len:
        chimp_5p_viable = True
    else:
        chimp_5p_viable = False
        print('chimp_5p: ' + str(human_5p_len) + ' > ' + str(chimp_5p_len))

    if chimp_3p_len >= human_3p_len:
        chimp_3p_viable = True
    else:
        chimp_3p_viable = False
        print('chimp_3p: ' + str(human_3p_len) + ' > ' + str(chimp_3p_len))
    
    if vervet_5p_len >= human_5p_len:
        vervet_5p_viable = True
    else:
        vervet_5p_viable = False
        print('vervet_5p: ' + str(human_5p_len) + ' > ' + str(vervet_5p_len))

    if vervet_3p_len >= human_3p_len:
        vervet_3p_viable = True
    else:
        vervet_3p_viable = False
        print('vervet_3p: ' + str(human_3p_len) + ' > ' + str(vervet_3p_len))
        
    if chimp_5p_viable == False or chimp_3p_viable == False or vervet_5p_viable == False or vervet_3p_viable == False:
        avoid_lines.append(counter)
    counter += 1

print(avoid_lines)




#This section removes any vervet exon that has some sort of issue with the alignment. 
#These alignment issues are detected using the ouput of the "find_gene_duplicates_and_incomplete_alignments.py" script 
#and need to be identified because they may cause issues when calculating divergence. Though 
#divergence values have already been obtained for these problematic regions by this point, those 
#problematic regions are removed here. These issues include having no vervet alignment data for 
#the vervet exon region, having incomplete alignment data that only spans part of the exon region, 
#having the region align to multiple different human, vervet, or marmoset gene duplicates, or having 
#no marmoset alignment in the exon region and thus having no vervet-human-marmoset ancestral sequence 
#to use in the calculation of the ancestral allele.

alignment_check_file = 'C:/Users/Kellen/Downloads/primate_dfe_project/chimp_and_vervet_data/248_exon_alignment_check.txt'
with open(alignment_check_file, 'r') as f:
    alignment_check = f.readlines()[1 : ]

counter = 0    
for line in alignment_check:
    split = line.strip('\n').split('\t')
#     print(split)
    if len([i for i in split if i != '']) < 5:
        avoid_lines.append(counter)
#         print(counter)
#         print('A')
        counter += 1
    else:
        human_range = split[0]
        vervet_range = split[1]
        num_blocks = split[2]
        human_blocks = split[4]
        vervet_blocks = split[5]
        gene_duplicates = split[6]
        source_listings = split[7]
        
        vervet_blocks = vervet_blocks.split('),')
        find_empty = []
        for z in vervet_blocks:
            if z == '(' or z == ' ()' or z == ' (':
                find_empty.append(vervet_blocks.index(z))       
        temp_vervet_blocks = []
        
        for i in range(0, len(vervet_blocks)):
            if i not in find_empty:
                temp_vervet_blocks.append(vervet_blocks[i])        
        vervet_blocks = list(temp_vervet_blocks)
        gene_duplicates = gene_duplicates.split('),')
        block_inclusions = []
        for q in gene_duplicates:
            q = q.replace('(', '')
            q = q.replace(')', '')
            q = q.replace(',', '')
            q = q.replace(' ', '')
            block_inclusions.append(q)

        temp_block_inclusions = []
        for i in range(0, len(block_inclusions)):
            if i not in find_empty:
                temp_block_inclusions.append(block_inclusions[i])        
        block_inclusions = list(temp_block_inclusions)

        all_blocks_okay = True
        for r in block_inclusions:
            H_count = r.count('H')
            V_count = r.count('V')
            C_count = r.count('C')
            if H_count != 1 or V_count != 1 or C_count != 1:
                all_blocks_okay = False

        if all_blocks_okay == True:
            vervet_start = int(vervet_range[0 : vervet_range.find('-')])
            vervet_end = int(vervet_range[vervet_range.find('-') + 1 : ])

            windows = []
            for v in vervet_blocks:
                start = v[ : v.find(',')].replace('(', '').replace(' ', '')
                end = v[v.find(',') + 1 : ].replace(')', '').replace(' ', '')
                pair = [int(start), int(end)]
                windows.append(pair)

            in_exon = set([x for x in range(vervet_start, vervet_end + 1)])
            in_alignment_region = []
            for pair in windows:
                start = pair[0]
                end = pair[1]
                for y in range(start, end + 1):
                    in_alignment_region.append(y)
            in_alignment_region = set(in_alignment_region)

            if in_exon.issubset(in_alignment_region) == False:
                avoid_lines.append(counter)
#                 print(counter)
#                 print('C')
                counter += 1
            else:
                counter += 1

        else:
            avoid_lines.append(counter)
#             print(counter)
#             print('B')
            counter += 1
            
alignment_check = []
print(avoid_lines)




#This section reads in information from Ensembl's BioMart database that we use to confirm if 
#exons are truly homologous. It filters out exons based on that based on that information, only 
#taking exons where the chimp and vervet exon have a homology type of 1 to 1, an orthology 
#confidence of 1, and a gene order conservation (GOC) score of 100

chimp_ensembl_file = 'C:/Users/Kellen/Downloads/primate_dfe_project/chimp_and_vervet_data/ensembl_biomart_homo_exons_chimp.txt'
with open(chimp_ensembl_file, 'r') as f:
    chimp_ensembl_data = f.readlines()[1 :]
    
vervet_ensembl_file = 'C:/Users/Kellen/Downloads/primate_dfe_project/chimp_and_vervet_data/ensembl_biomart_homo_exons_vervet.txt'
with open(vervet_ensembl_file, 'r') as f:
    vervet_ensembl_data = f.readlines()[1 :]

    

chimp_chr_nums = ['1', '2A', '2B', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22']
chimp_chr_dic_gene = {}
for num in chimp_chr_nums:
    chimp_chr_dic_gene[num] = []
    
chimp_chr_dic_gene_remove_dups = {}
for num in chimp_chr_nums:
    chimp_chr_dic_gene_remove_dups[num] = []
    
for line in chimp_ensembl_data:
    line = line.strip('\n')
    if '\t\t' in line:
        line = line.replace('\t\t', '\tNA\t')
    split = line.split()
    if ('ENSP' in split[-1] or 'ENST' in split[-1]) and (isinstance(split[-1], str) == True):
        split.append('NA')

    chr_num = split[7] 
#     print(chr_num)
    if chr_num in chimp_chr_nums:
        gene_name = split[4]
        homology_type = split[10]
        goc_score = split[11]
        orthology_confidence = split[12]
        info = [gene_name, homology_type, goc_score, orthology_confidence]
        chimp_chr_dic_gene[chr_num].append(info)
    else:
        continue
        
for num in chimp_chr_nums:        
    for w in chimp_chr_dic_gene[num]: 
        if w not in chimp_chr_dic_gene_remove_dups[num]: 
            chimp_chr_dic_gene_remove_dups[num].append(w)          
     
    
    
vervet_chr_nums = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22','23','24','25','26','27','28','29']
vervet_chr_dic_gene = {}
for num in vervet_chr_nums:
    vervet_chr_dic_gene[num] = []
    
vervet_chr_dic_gene_remove_dups = {}
for num in vervet_chr_nums:
    vervet_chr_dic_gene_remove_dups[num] = []
    
for line in vervet_ensembl_data:
    line = line.strip('\n')
    if '\t\t' in line:
        line = line.replace('\t\t', '\tNA\t')
    split = line.split()
    if ('ENST' in split[-1]) and (isinstance(split[-1], str) == True):
        split.append('NA')

    chr_num = split[5] 
#     print(chr_num)
    if chr_num in vervet_chr_nums:
        gene_name = split[12]
        homology_type = split[8]
        goc_score = split[9]
        orthology_confidence = split[10]
        info = [gene_name, homology_type, goc_score, orthology_confidence]
        vervet_chr_dic_gene[chr_num].append(info)
    else:
        continue
        
for num in vervet_chr_nums:        
    for w in vervet_chr_dic_gene[num]: 
        if w not in vervet_chr_dic_gene_remove_dups[num]: 
            vervet_chr_dic_gene_remove_dups[num].append(w)

            
            
num_exons = len(homo_exon_data)
for p in range(0, num_exons): 
    split = homo_exon_data[p].strip('\n').split()
    
    chimp_chr_num = split[24]
    chimp_chr_num = chimp_chr_num.replace('chr','')
    if chimp_chr_num == '2a':
        chimp_chr_num = '2A'
    if chimp_chr_num == '2b':
        chimp_chr_num = '2B'
    chimp_gene_name = split[23]
    chimp_find = chimp_gene_name.find('_')
    chimp_gene_name = chimp_gene_name[0 : chimp_find]
    
    vervet_chr_num = split[51]
    vervet_chr_num = vervet_chr_num.replace('chr','')
    vervet_gene_name = split[50]
    vervet_find = vervet_gene_name.find('.')
    vervet_gene_name = vervet_gene_name[0 : vervet_find]
    
#     print(chimp_chr_num)
#     print(chimp_gene_name)
#     print(vervet_chr_num)
#     print(vervet_gene_name)
    
    is_chimp_good = False
    for entry in chimp_chr_dic_gene_remove_dups[chimp_chr_num]:
        entry_name = entry[0]
#         print(entry_name)
        if chimp_gene_name == entry_name:
            chimp_homology_type = entry[1]
            chimp_goc_score = entry[2]
            chimp_orthology_confidence = entry[3]
#             print(homology_type)
#             print(chimp_goc_score)
#             print(orthology_confidence)
            if ('one2one' in chimp_homology_type) and chimp_goc_score == '100' and chimp_orthology_confidence == '1':
                is_chimp_good = True
    
    is_vervet_good = False
    for entry in vervet_chr_dic_gene_remove_dups[vervet_chr_num]:
        entry_name = entry[0]
        if vervet_gene_name == entry_name:
            vervet_homology_type = entry[1]
            vervet_goc_score = entry[2]
            vervet_orthology_confidence = entry[3]
            if ('one2one' in vervet_homology_type) and vervet_goc_score == '100' and vervet_orthology_confidence == '1':
                is_vervet_good = True
                
    if is_chimp_good == False or is_vervet_good == False:
        avoid_lines.append(p)
    
print(avoid_lines)




#This section copies any exons that passed through all of the filters to a new file.

avoid_lines = list(set(avoid_lines))
avoid_lines.sort()
print(str(len(homo_exon_data) - len(avoid_lines)) + ' left')

new_filtered_homo_exon_file = 'C:/Users/Kellen/Downloads/primate_dfe_project/chimp_and_vervet_data/FILTERED_homo_exons_reference_all_info.txt'
results_file = open(new_filtered_homo_exon_file, 'w+')
results_file.write(homo_exon_data_column_names)

new_number = 1
for q in range(0, len(homo_exon_data)):
    if q not in avoid_lines:
        split = homo_exon_data[q].strip('\n').split()[1 : ]
        new_line = str(new_number) + '\t' + '\t'.join(split) + '\n'
        results_file.write(new_line)
        new_number += 1

results_file.close()

