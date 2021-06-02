


#This is a Python script that is used to graph recombination rate correlations for the list of homologous 
#exons we have between humans, vervets, and chimps.

#This script is currently designed to be run locally and will require path changes on the user's part to the variables:
#'all_primates_file' with input file "FILTERED_homo_exons_reference_all_info.txt"
#'human_vs_chimp_rec_graph' with output file 'human-chimp_recombination_graph' + file_tag + '.pdf'
#'human_vs_vervet_rec_graph' with output file 'human-vervet_recombination_graph' + file_tag + '.pdf'
#'chimp_vs_vervet_rec_graph' with output file 'chimp-vervet_recombination_graph' + file_tag + '.pdf'




import matplotlib.pyplot as plt
import sys




#Note the variable minimum_exon_length. By changing this value, you can change the minimum 
#size an exon needs to be for it to be considered in the plot. Currently, either a 
#minimum size of 0bp or 4000bp are possible, but others can be easily added.

homo_exons_file = 'C:/Users/Kellen/Downloads/primate_dfe_project/chimp_and_vervet_data/FILTERED_homo_exons_reference_all_info.txt'
with open(homo_exons_file, 'r') as f:
    homo_exons_data = f.readlines()[1:]

human_rec_data = []
chimp_rec_data = []
vervet_rec_data = []

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
    
    human_rec = float(split[20])
    chimp_rec = float(split[45])
    vervet_rec = float(split[72])

    human_exon_len = int(split[8]) + 1
    if human_exon_len >= minimum_exon_len:
        human_rec_data.append(human_rec)      
        chimp_rec_data.append(chimp_rec)      
        vervet_rec_data.append(vervet_rec)




primary_folder = 'C:/Users/Kellen/Downloads/primate_dfe_project/chimp_and_vervet_data/raw_plots/'

human_vs_chimp_rec_graph = primary_folder + 'human-chimp_recombination_graph' + file_tag + '.pdf'
human_vs_vervet_rec_graph = primary_folder + 'human-vervet_recombination_graph' + file_tag + '.pdf'
chimp_vs_vervet_rec_graph = primary_folder + 'chimp-vervet_recombination_graph' + file_tag + '.pdf'

fig = plt.figure(figsize = (20, 11.3))
ax = fig.add_subplot(111)
plt.rc('xtick', labelsize = 40)
ax.xaxis.set_tick_params(length = 8, width = 3)
ax.xaxis.set_tick_params(length = 5, width = 2, which = 'minor')
plt.rc('ytick', labelsize = 40)
ax.yaxis.set_tick_params(length = 8, width = 3)
ax.yaxis.set_tick_params(length = 5, width = 2, which = 'minor')
plt.close()

do = [1, 2, 3]
# do=[9]
for i in do:
    if i == 1:
        fig = plt.figure(figsize = (20, 15))
        ax = fig.add_subplot(111)
        ax.scatter(human_rec_data, chimp_rec_data)
        ax.set_xlim([0, 10])
        ax.set_ylim([0, 3])
        ax.set_xlabel('Human Recombination', fontsize = 40)
        ax.set_ylabel('Chimp Recombination (4Ner/kb)', fontsize = 40)
        fig.savefig(human_vs_chimp_rec_graph)
    elif i == 2:
        fig = plt.figure(figsize = (20, 15))
        ax = fig.add_subplot(111)
        ax.scatter(human_rec_data, vervet_rec_data)
        ax.set_xlim([0, 10])
        ax.set_ylim([0, 3])
        ax.set_xlabel('Human Recombination', fontsize = 40)
        ax.set_ylabel('Vervet Recombination (4Ner/kb)', fontsize = 40)
        fig.savefig(human_vs_vervet_rec_graph)
    elif i == 3:
        fig = plt.figure(figsize = (20, 15))
        ax = fig.add_subplot(111)
        ax.scatter(chimp_rec_data, vervet_rec_data)
        ax.set_xlim([0, 3])
        ax.set_ylim([0, 3])
        ax.set_xlabel('Chimp Recombination (4Ner/kb)', fontsize = 40)
        ax.set_ylabel('Vervet Recombination (4Ner/kb)', fontsize = 40)
        fig.savefig(chimp_vs_vervet_rec_graph)

