


#This is a Python script that is used to graph exon length correlations for the list of homologous 
#exons we have between humans, vervets, and chimps, as well as exon length distributions for all 
#exons for vervets and chimps.

#This script is currently designed to be run locally and will require path changes on the user's part to the variables:
#'all_primates_file' with input file "FILTERED_homo_exons_reference_all_info.txt"
#'human_chimp_graph' with output file "human_vs_chimp_exon_len_graph.pdf"
#'human_vervet_graph' with output file "human_vs_vervet_exon_len_graph.pdf"
#'chimp_vervet_graph' with output file "chimp_vs_vervet_exon_graph.pdf"




import matplotlib.pyplot as plt




all_primates_file = 'C:/Users/Kellen/Downloads/primate_dfe_project/chimp_and_vervet_data/FILTERED_homo_exons_reference_all_info.txt'
with open(all_primates_file) as g:
    homo_exon_data = g.readlines()[1 : ]
    
human_exons = []
chimp_exons = []
vervet_exons = []

for i in homo_exon_data:
    split = i.split()

    human_len = float(split[8])
    chimp_len = float(split[30])
    vervet_len = float(split[57])

    human_exons.append(human_len)
    chimp_exons.append(chimp_len)
    vervet_exons.append(vervet_len)




human_chimp_graph = 'C:/Users/Kellen/Downloads/primate_dfe_project/chimp_and_vervet_data/raw_plots/human_vs_chimp_exon_len_graph.pdf'
human_vervet_graph = 'C:/Users/Kellen/Downloads/primate_dfe_project/chimp_and_vervet_data/raw_plots/human_vs_vervet_exon_len_graph.pdf'
chimp_vervet_graph = 'C:/Users/Kellen/Downloads/primate_dfe_project/chimp_and_vervet_data/raw_plots/chimp_vs_vervet_exon_graph.pdf'

for i in range(0, 3):
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
        ax.scatter(human_exons, chimp_exons)
        ax.set_xlim([0, 6000])
        ax.set_ylim([0, 6000])
        ax.set_xlabel('Human exon length', fontsize = 40)
        ax.set_ylabel('Chimp exon length', fontsize = 40)
        plt.show()
        fig.savefig(human_chimp_graph)

    if i == 1:
        ax.scatter(human_exons, vervet_exons)
        ax.set_xlim([0, 6000])
        ax.set_ylim([0, 6000])
        ax.set_xlabel('Human exon length', fontsize = 40)
        ax.set_ylabel('Vervet exon length', fontsize = 40)
        plt.show()
        fig.savefig(human_vervet_graph)

    if i == 2:
        ax.scatter(chimp_exons, vervet_exons)
        ax.set_xlim([0, 6000])
        ax.set_ylim([0, 6000])
        ax.set_xlabel('Chimp exon length', fontsize = 40)
        ax.set_ylabel('Vervet exon length', fontsize = 40)
        plt.show()
        fig.savefig(chimp_vervet_graph)






