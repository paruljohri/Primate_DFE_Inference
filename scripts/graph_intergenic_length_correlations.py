


#This is a Python script that is used to graph intergenic length correlations for the list of homologous 
#exons we have between humans, vervets, and chimps.

#This script is currently designed to be run locally and will require path changes on the user's part to the variables:
#'all_primates_file' with input file "FILTERED_homo_exons_reference_all_info.txt"
#'human_chimp_5p_graph' with output file "human_vs_chimp_5p_len_graph.pdf"
#'human_chimp_3p_graph' with output file "human_vs_chimp_3p_len_graph.pdf"
#'human_vervet_5p_graph' with output file "human_vs_vervet_5p_len_graph.pdf"
#'human_vervet_3p_graph' with output file "human_vs_vervet_3p_len_graph.pdf"
#'chimp_vervet_5p_graph' with output file "chimp_vs_vervet_5p_len_graph.pdf"
#'chimp_vervet_3p_graph' with output file "chimp_vs_vervet_3p_len_graph.pdf"




import matplotlib.pyplot as plt




all_primates_file = 'C:/Users/Kellen/Downloads/primate_dfe_project/chimp_and_vervet_data/FILTERED_homo_exons_reference_all_info.txt'
with open(all_primates_file) as g:
    homo_exon_data = g.readlines()[1 :]
    
human_5p = []
human_3p = []
chimp_5p = []
chimp_3p = []
vervet_5p = []
vervet_3p = []

for i in homo_exon_data:
    split = i.split()

    human_5p_len = float(split[11])
    human_3p_len = float(split[14])
    chimp_5p_len = float(split[33])
    chimp_3p_len = float(split[36])
    vervet_5p_len = float(split[60])
    vervet_3p_len = float(split[63])

    human_5p.append(human_5p_len)
    human_3p.append(human_3p_len)
    chimp_5p.append(chimp_5p_len)
    chimp_3p.append(chimp_3p_len)
    vervet_5p.append(vervet_5p_len)
    vervet_3p.append(vervet_3p_len)




human_chimp_5p_graph = 'C:/Users/Kellen/Downloads/primate_dfe_project/chimp_and_vervet_data/raw_plots/human_vs_chimp_5p_len_graph.pdf'
human_chimp_3p_graph = 'C:/Users/Kellen/Downloads/primate_dfe_project/chimp_and_vervet_data/raw_plots/human_vs_chimp_3p_len_graph.pdf'
human_vervet_5p_graph = 'C:/Users/Kellen/Downloads/primate_dfe_project/chimp_and_vervet_data/raw_plots/human_vs_vervet_5p_len_graph.pdf'
human_vervet_3p_graph = 'C:/Users/Kellen/Downloads/primate_dfe_project/chimp_and_vervet_data/raw_plots/human_vs_vervet_3p_len_graph.pdf'
chimp_vervet_5p_graph = 'C:/Users/Kellen/Downloads/primate_dfe_project/chimp_and_vervet_data/raw_plots/chimp_vs_vervet_5p_len_graph.pdf'
chimp_vervet_3p_graph = 'C:/Users/Kellen/Downloads/primate_dfe_project/chimp_and_vervet_data/raw_plots/chimp_vs_vervet_3p_len_graph.pdf'

for i in range(0, 6):
# for i in [1]:
    fig = plt.figure(figsize = (23, 13))
    ax = fig.add_subplot(111)
    plt.rc('xtick', labelsize = 40)
    ax.xaxis.set_tick_params(length = 8, width = 3)
    ax.xaxis.set_tick_params(length = 5, width = 2, which = 'minor')
    plt.rc('ytick', labelsize = 40)
    ax.yaxis.set_tick_params(length = 8, width = 3)
    ax.yaxis.set_tick_params(length = 5, width = 2, which = 'minor')

    if i == 0:
        ax.scatter(human_5p, chimp_5p)
        ax.set_xlim([0, 100000])
        ax.set_ylim([0, 100000])
        ax.set_xlabel('Human 5\' intergenic length', fontsize = 40)
        ax.set_ylabel('Chimp 5\' intergenic length', fontsize = 40)
        plt.show()
        fig.savefig(human_chimp_5p_graph)
    if i == 1:
        ax.scatter(human_3p, chimp_3p)
        ax.set_xlim([0, 100000])
        ax.set_ylim([0, 100000])
        ax.set_xlabel('Human 3\' intergenic length', fontsize = 40)
        ax.set_ylabel('Chimp 3\' intergenic length', fontsize = 40)
        plt.show()
        fig.savefig(human_chimp_3p_graph)
    if i == 2:
        ax.scatter(human_5p, vervet_5p)
        ax.set_xlim([0, 100000])
        ax.set_ylim([0, 100000])
        ax.set_xlabel('Human 5\' intergenic length', fontsize = 40)
        ax.set_ylabel('Vervet 5\' intergenic length', fontsize = 40)
        plt.show()
        fig.savefig(human_vervet_5p_graph)
    if i == 3:
        ax.scatter(human_3p, vervet_3p)
        ax.set_xlim([0, 100000])
        ax.set_ylim([0, 100000])
        ax.set_xlabel('Human 3\' intergenic length', fontsize = 40)
        ax.set_ylabel('Vervet 3\' intergenic length', fontsize = 40)
        plt.show()
        fig.savefig(human_vervet_3p_graph)
    if i == 4:
        ax.scatter(chimp_5p, vervet_5p)
        ax.set_xlim([0, 100000])
        ax.set_ylim([0, 100000])
        ax.set_xlabel('Chimp 5\' intergenic length', fontsize = 40)
        ax.set_ylabel('Vervet 5\' intergenic length', fontsize = 40)
        plt.show()
        fig.savefig(chimp_vervet_5p_graph)
    if i == 5:
        ax.scatter(chimp_3p, vervet_3p)
        ax.set_xlim([0, 100000])
        ax.set_ylim([0, 100000])
        ax.set_xlabel('Chimp 3\' intergenic length', fontsize = 40)
        ax.set_ylabel('Vervet 3\' intergenic length', fontsize = 40)
        plt.show()
        fig.savefig(chimp_vervet_3p_graph)






