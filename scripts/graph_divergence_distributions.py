


#This is a Python script that is used to graph divergence distributions for the list of human exons we have
#identified, as well for all exons in the chimp and vervet UCSC exon annotations.

#This script is currently designed to be run locally and will require path changes on the user's part to the variables:
#'all_chimp_div_file' with input file "chimp_all_exons_divergence.txt"
#'all_vervet_div_file' with input file "vervet_all_exons_divergence.txt"
#'homo_div_file' with input file "FILTERED_homo_exons_reference_all_info.txt"
#'chimp_div_graph' with input file "raw_plots/chimp_exon_divergence_dist.pdf"
#'vervet_div_graph' with input file "raw_plots/vervet_exon_divergence_dist.pdf"
#'human_div_graph' with input file "raw_plots/human_exon_divergence_dist.pdf"




import matplotlib.pyplot as plt
import numpy as np




main_folder = 'C:/Users/Kellen/Downloads/primate_dfe_project/chimp_and_vervet_data'
homo_human_div = []
homo_chimp_div = []
homo_vervet_div = []
all_chimp_div = []
all_vervet_div = []

all_chimp_div_file = main_folder + '/' + 'chimp_all_exons_divergence.txt'
all_vervet_div_file = main_folder + '/' + 'vervet_all_exons_divergence.txt'
homo_div_file = main_folder + '/' + 'FILTERED_homo_exons_reference_all_info.txt'

chimp_div_graph = main_folder + '/' + 'raw_plots/chimp_exon_divergence_dist.pdf'
vervet_div_graph = main_folder + '/' + 'raw_plots/vervet_exon_divergence_dist.pdf'
human_div_graph = main_folder + '/' + 'raw_plots/human_exon_divergence_dist.pdf'

with open(homo_div_file, 'r') as f:
    homo_div_data = f.readlines()[1 :]
    for i in homo_div_data:
        split = i.strip('\n').split()
        homo_human_div_value = float(split[15])
        homo_chimp_div_value = float(split[37])
        homo_vervet_div_value = float(split[64])
        
        homo_human_div.append(homo_human_div_value)
        homo_chimp_div.append(homo_chimp_div_value)
        homo_vervet_div.append(homo_vervet_div_value)
            
with open(all_chimp_div_file, 'r') as f:
    all_chimp_div_data = f.readlines()[1 :]
    for i in all_chimp_div_data:
        split = i.strip('\n').split()
        if split[5] != 'NA':
            all_chimp_div_value = float(split[5])
            all_chimp_div.append(all_chimp_div_value)
    
with open(all_vervet_div_file, 'r') as f:
    all_vervet_div_data = f.readlines()[1 :]
    for i in all_vervet_div_data:
        split = i.strip('\n').split()
        if split[5] != 'NA':
            all_vervet_div_value = float(split[5])
            all_vervet_div.append(all_vervet_div_value)




for i in range(0, 3):
    if i == 0:
        print('Chimp-ALL')
        print('min: ' + str(min(all_chimp_div)))
        print('max: ' + str(max(all_chimp_div)))
        print('median: ' + str(np.median(all_chimp_div)))
        print('mean: ' + str(np.mean(all_chimp_div)))
        print('std dev: ' + str(np.std(all_chimp_div)))
        print('Chimp-HOMO')
        print('min: ' + str(min(homo_chimp_div)))
        print('max: ' + str(max(homo_chimp_div)))
        print('median: ' + str(np.median(homo_chimp_div)))
        print('mean: ' + str(np.mean(homo_chimp_div)))
        print('std dev: ' + str(np.std(homo_chimp_div)))
        
        bins = []
        start = 0
        bin_size = 0.001
        max_x_range = 0.1
        num_bins = int(max_x_range / bin_size)
        for i in range(0, num_bins + 1):
            bins.append(start)
            start = start + bin_size
        
        plt.hist(all_chimp_div, bins, density = True, stacked = True, fill = False, alpha = 0.6, ec = 'blue')
        plt.hist(homo_chimp_div, bins, density = True, stacked = True,  fill = False, alpha = 0.5, ec = 'green')
        plt.gca().set(xlabel = 'Chimp Divergence', ylabel = 'Probability \n (÷ by bin width of ' + str(bin_size) + ')', ylim = [0, 100]);
        plt.savefig(chimp_div_graph)
        plt.show()
        
    if i == 1:        
        print('Vervet-ALL')
        print('min: ' + str(min(all_vervet_div)))
        print('max: ' + str(max(all_vervet_div)))
        print('median: ' + str(np.median(all_vervet_div)))
        print('mean: ' + str(np.mean(all_vervet_div)))
        print('std dev: ' + str(np.std(all_vervet_div)))
        print('Vervet-HOMO')
        print('min: ' + str(min(homo_vervet_div)))
        print('max: ' + str(max(homo_vervet_div)))
        print('median: ' + str(np.median(homo_vervet_div)))
        print('mean: ' + str(np.mean(homo_vervet_div)))
        print('std dev: ' + str(np.std(homo_vervet_div)))
        
        bins = []
        start = 0
        bin_size = 0.001
        max_x_range = 0.1
#         bin_size=0.001
#         max_x_range=0.35
        num_bins = int(max_x_range / bin_size)
        for i in range(0, num_bins + 1):
            bins.append(start)
            start = start + bin_size
            
        plt.hist(all_vervet_div, bins, density = True, stacked = True, fill = False, alpha = 0.6, ec = 'blue')
        plt.hist(homo_vervet_div, bins, density = True, stacked = True, fill = False, alpha = 0.5, ec = 'green')
        plt.gca().set(xlabel = 'Vervet Divergence', ylabel = 'Probability \n (÷ by bin width of ' + str(bin_size) + ')', ylim = [0, 100]);
        plt.savefig(vervet_div_graph)
        plt.show()
        
    if i == 2:
        print('Human-HOMO')
        print('min: ' + str(min(homo_human_div)))
        print('max: ' + str(max(homo_human_div)))
        print('median: ' + str(np.median(homo_human_div)))
        print('mean: ' + str(np.mean(homo_human_div)))
        print('std dev: ' + str(np.std(homo_human_div)))
        
        bins = []
        start = 0
        bin_size = 0.0001
        max_x_range = 0.005
        num_bins = int(max_x_range / bin_size)
        for i in range(0, num_bins + 1):
            bins.append(start)
            start = start + bin_size
            
        plt.hist(homo_human_div, bins, density = True, stacked = True, fill = False, alpha = 0.5, ec = 'green')
        plt.gca().set(xlabel = 'Human Divergence', ylabel = 'Probability \n (÷ by bin width of ' + str(bin_size) + ')');
        plt.savefig(human_div_graph)
        plt.show()






