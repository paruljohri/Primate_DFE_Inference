


#This is a Python script that is used to graph both exonic and intergenic divergence correlations for the list of homologous 
#exons we have between humans, vervets, and chimps.

#This script is currently designed to be run locally and will require path changes on the user's part to the variables:
#'homo_exons_file' with input file "FILTERED_homo_exons_reference_all_info.txt"
#'human_vs_chimp_exon_div_graph'with input file 'human-chimp_exon_divergence_graph' + file_tag + '.pdf'
#'human_vs_vervet_exon_div_graph' with input file 'human-vervet_exon_divergence_graph' + file_tag + '.pdf'
#'chimp_vs_vervet_exon_div_graph' with input file 'chimp-vervet_exon_divergence_graph' + file_tag + '.pdf'
#'human_vs_chimp_5p_div_graph' with input file 'human-chimp_5p_divergence_graph' + file_tag + '.pdf'
#'human_vs_vervet_5p_div_graph' with input file 'human-vervet_5p_divergence_graph' + file_tag + '.pdf'
#'chimp_vs_vervet_5p_div_graph' with input file 'chimp-vervet_5p_divergence_graph' + file_tag + '.pdf'
#'human_vs_chimp_3p_div_graph' with input file 'human-chimp_3p_divergence_graph' + file_tag + '.pdf'
#'human_vs_vervet_3p_div_graph' with input file 'human-vervet_3p_divergence_graph' + file_tag + '.pdf'
#'chimp_vs_vervet_3p_div_graph' with input file 'chimp-vervet_3p_divergence_graph' + file_tag + '.pdf'




import matplotlib.pyplot as plt
import sys




#Note the variable minimum_exon_length. By changing this value, you can change the minimum 
#size an exon needs to be for it to be considered in the plot. Currently, either a 
#minimum size of 0bp or 4000bp are possible, but others can be easily added.

homo_exons_file = 'C:/Users/Kellen/Downloads/primate_dfe_project/chimp_and_vervet_data/FILTERED_homo_exons_reference_all_info.txt'
with open(homo_exons_file, 'r') as f:
    homo_exons_data = f.readlines()[1 :]

human_exon_div_data = []
human_5p_div_data = []
human_3p_div_data = []
chimp_exon_div_data = []
chimp_5p_div_data = []
chimp_3p_div_data = []
vervet_exon_div_data = []
vervet_5p_div_data = []
vervet_3p_div_data = []

minimum_exon_len = 0
if minimum_exon_len == 0:
    file_tag = ''
elif minimum_exon_len == 4000:
    file_tag = '_gt4kb'
else:
    print('ERROR! There\'s no file name tag for the current minimum exon length filter\n\n\n')
    sys.exit()

for i in homo_exons_data[:]:
    split = i.strip('\n').split('\t')
    
    human_exon_div = float(split[15])
#     human_5p_div=float(split[16])
#     human_3p_div=float(split[17])
    
    chimp_exon_div = float(split[35])
    chimp_5p_div = float(split[37])
    chimp_3p_div = float(split[39])
    
    vervet_exon_div = float(split[61])
    vervet_5p_div = float(split[63])
    vervet_3p_div = float(split[65])
    
    human_exon_len = int(split[8]) + 1
    if human_exon_len >= minimum_exon_len:
        human_exon_div_data.append(human_exon_div)
#         human_5p_div_data.append(human_5p_div)
#         human_3p_div_data.append(human_3p_div)
        
        chimp_exon_div_data.append(chimp_exon_div)
        chimp_5p_div_data.append(chimp_5p_div)
        chimp_3p_div_data.append(chimp_3p_div)
        
        vervet_exon_div_data.append(vervet_exon_div)
        vervet_5p_div_data.append(vervet_5p_div)
        vervet_3p_div_data.append(vervet_3p_div)




primary_folder = 'C:/Users/Kellen/Downloads/primate_dfe_project/chimp_and_vervet_data/raw_plots/'
human_vs_chimp_exon_div_graph = primary_folder + 'human-chimp_exon_divergence_graph' + file_tag + '.pdf'
human_vs_vervet_exon_div_graph = primary_folder + 'human-vervet_exon_divergence_graph' + file_tag + '.pdf'
chimp_vs_vervet_exon_div_graph = primary_folder + 'chimp-vervet_exon_divergence_graph' + file_tag + '.pdf'

human_vs_chimp_5p_div_graph = primary_folder + 'human-chimp_5p_divergence_graph' + file_tag + '.pdf'
human_vs_vervet_5p_div_graph = primary_folder + 'human-vervet_5p_divergence_graph' + file_tag + '.pdf'
chimp_vs_vervet_5p_div_graph = primary_folder + 'chimp-vervet_5p_divergence_graph' + file_tag + '.pdf'

human_vs_chimp_3p_div_graph = primary_folder + 'human-chimp_3p_divergence_graph' + file_tag + '.pdf'
human_vs_vervet_3p_div_graph = primary_folder + 'human-vervet_3p_divergence_graph' + file_tag + '.pdf'
chimp_vs_vervet_3p_div_graph = primary_folder + 'chimp-vervet_3p_divergence_graph' + file_tag + '.pdf'

fig = plt.figure(figsize = (20, 11.3))
ax = fig.add_subplot(111)
plt.rc('xtick', labelsize = 40)
ax.xaxis.set_tick_params(length = 8, width = 3)
ax.xaxis.set_tick_params(length = 5, width = 2, which = 'minor')
plt.rc('ytick', labelsize = 40)
ax.yaxis.set_tick_params(length = 8, width = 3)
ax.yaxis.set_tick_params(length = 5, width = 2, which = 'minor')
plt.close()

do = [1, 2, 3, 6, 9]
# do=[9]
for i in do:
    if i == 1:
        fig = plt.figure(figsize = (20, 15))
        ax = fig.add_subplot(111)
        ax.scatter(human_exon_div_data, chimp_exon_div_data)
        ax.set_xlim([0, 0.012])
        ax.set_ylim([0, 0.012])
        ax.set_xlabel('Human Exon Divergence', fontsize = 40)
        ax.set_ylabel('Chimp Exon Divergence', fontsize = 40)
        fig.savefig(human_vs_chimp_exon_div_graph)
    elif i == 2:
        fig = plt.figure(figsize = (20, 15))
        ax = fig.add_subplot(111)
        ax.scatter(human_exon_div_data, vervet_exon_div_data)
        ax.set_xlim([0, 0.006])
        ax.set_ylim([0, 0.05])
        ax.set_xlabel('Human Exon Divergence', fontsize = 40)
        ax.set_ylabel('Vervet Exon Divergence', fontsize = 40)
        fig.savefig(human_vs_vervet_exon_div_graph)
    elif i == 3:
        fig = plt.figure(figsize = (20, 15))
        ax = fig.add_subplot(111)
        ax.scatter(chimp_exon_div_data, vervet_exon_div_data)
        ax.set_xlim([0, 0.006])
        ax.set_ylim([0, 0.05])
#         ax.set_xlim([0, 0.05])
#         ax.set_ylim([0, 0.05])
        ax.set_xlabel('Chimp Exon Divergence', fontsize = 40)
        ax.set_ylabel('Vervet Exon Divergence', fontsize = 40)
        fig.savefig(chimp_vs_vervet_exon_div_graph)  
#     elif i == 4:
#         fig = plt.figure(figsize = (20, 15))
#         ax = fig.add_subplot(111)
#         ax.scatter(human_5p_div_data, chimp_5p_div_data)
#         ax.set_xlim([0, 0.006])
#         ax.set_ylim([0, 0.006])
#         ax.set_xlabel('Human 5\' Divergence', fontsize = 40)
#         ax.set_ylabel('Chimp 5\' Divergence', fontsize = 40)
#         fig.savefig(human_vs_chimp_5p_div_graph)
#     elif i == 5:
#         fig = plt.figure(figsize = (20, 15))
#         ax = fig.add_subplot(111)
#         ax.scatter(human_5p_div_data, vervet_5p_div_data)
#         ax.set_xlim([0, 0.006])
#         ax.set_ylim([0, 0.05])
#         ax.set_xlabel('Human 5\' Divergence', fontsize = 40)
#         ax.set_ylabel('Vervet 5\' Divergence', fontsize = 40)
#         fig.savefig(human_vs_vervet_5p_div_graph)
    elif i == 6:
        fig = plt.figure(figsize = (20, 15))
        ax = fig.add_subplot(111)
        ax.scatter(chimp_5p_div_data, vervet_5p_div_data)
        ax.set_xlim([0, 0.005])
        ax.set_ylim([0, 0.05])
        ax.set_xlabel('Chimp 5\' Divergence', fontsize = 40)
        ax.set_ylabel('Vervet 5\' Divergence', fontsize = 40)
        fig.savefig(chimp_vs_vervet_5p_div_graph)
#     if i == 7:
#         fig = plt.figure(figsize = (20, 15))
#         ax = fig.add_subplot(111)
#         ax.scatter(human_3p_div_data, chimp_3p_div_data)
#         ax.set_xlim([0, 0.006])
#         ax.set_ylim([0, 0.006])
#         ax.set_xlabel('Human 3\' Divergence', fontsize = 40)
#         ax.set_ylabel('Chimp 3\' Divergence', fontsize = 40)
#         fig.savefig(human_vs_chimp_3p_div_graph)
#     elif i == 8:
#         fig = plt.figure(figsize = (20, 15))
#         ax = fig.add_subplot(111)
#         ax.scatter(human_3p_div_data, vervet_3p_div_data)
#         ax.set_xlim([0, 0.006])
#         ax.set_ylim([0, 0.05])
#         ax.set_xlabel('Human 3\' Divergence', fontsize = 40)
#         ax.set_ylabel('Vervet 3\' Divergence', fontsize = 40)
#         fig.savefig(human_vs_vervet_3p_div_graph)
    elif i == 9:
        fig = plt.figure(figsize = (20, 15))
        ax = fig.add_subplot(111)
        ax.scatter(chimp_3p_div_data, vervet_3p_div_data)
        ax.set_xlim([0, 0.005])
        ax.set_ylim([0, 0.05])
        ax.set_xlabel('Chimp 3\' Divergence', fontsize = 40)
        ax.set_ylabel('Vervet 3\' Divergence', fontsize = 40)
        fig.savefig(chimp_vs_vervet_3p_div_graph)

