


#This is a Python script that is used to graph recombination distributions for the list of human exons we have
#identified, as well for all exons in the chimp and vervet UCSC exon annotations.

#This script is currently designed to be run locally and will require path changes on the user's part to the variables:
#'all_chimp_rec_file' with input file "chimp_all_recombination.txt"
#'all_vervet_rec_file' with input file "vervet_all_recombination.txt"
#'homo_rec_file' with input file "FILTERED_homo_exons_reference_all_info.txt"
#'chimp_rec_graph' with output file "chimp_recombination_dist.pdf"
#'vervet_rec_graph' with output file "vervet_recombination_dist.pdf"
#'human_rec_graph' with output file "human_recombination_dist.pdf"




import matplotlib.pyplot as plt
import numpy as np




main_folder = 'C:/Users/Kellen/Downloads/primate_dfe_project/chimp_and_vervet_data'
homo_human_rec = []
homo_chimp_rec = []
homo_vervet_rec = []
all_chimp_rec = []
all_vervet_rec = []

all_chimp_rec_file = main_folder + '/chimp_all_recombination.txt'
all_vervet_rec_file = main_folder + '/vervet_all_recombination.txt'
homo_rec_file = main_folder + '/FILTERED_homo_exons_reference_all_info.txt'

chimp_rec_graph = main_folder + '/raw_plots/chimp_recombination_dist.pdf'
vervet_rec_graph = main_folder + '/raw_plots/vervet_recombination_dist.pdf'
human_rec_graph = main_folder + '/raw_plots/human_recombination_dist.pdf'

with open(homo_rec_file, 'r') as f:
    homo_rec_data = f.readlines()[1 :]
    for i in homo_rec_data:
        split = i.strip('\n').split()
        homo_human_rec_value = float(split[18])
        homo_chimp_rec_value = float(split[43])
        homo_vervet_rec_value = float(split[69])
        
        homo_human_rec.append(homo_human_rec_value)
        homo_chimp_rec.append(homo_chimp_rec_value)
        homo_vervet_rec.append(homo_vervet_rec_value)
            
with open(all_chimp_rec_file, 'r') as f:
    all_chimp_rec_data = f.readlines()[1 :]
    for i in all_chimp_rec_data:
        split = i.strip('\n').split()
        if split[6] != 'NA':
            all_chimp_rec_value = float(split[6])
            all_chimp_rec.append(all_chimp_rec_value)
    
with open(all_vervet_rec_file, 'r') as f:
    all_vervet_rec_data = f.readlines()[1 :]
    for i in all_vervet_rec_data:
        split = i.strip('\n').split()
        if split[6] != 'NA':
            all_vervet_rec_value = float(split[6])
            all_vervet_rec.append(all_vervet_rec_value)




for i in range(0, 3):
    if i == 0:
        print('Chimp-ALL')
        print('min: ' + str(min(all_chimp_rec)))
        print('max: ' + str(max(all_chimp_rec)))
        print('median: ' + str(np.median(all_chimp_rec)))
        print('mean: ' + str(np.mean(all_chimp_rec)))
        print('std dev: ' + str(np.std(all_chimp_rec)))
        print('Chimp-HOMO')
        print('min: ' + str(min(homo_chimp_rec)))
        print('max: ' + str(max(homo_chimp_rec)))
        print('median: ' + str(np.median(homo_chimp_rec)))
        print('mean: ' + str(np.mean(homo_chimp_rec)))
        print('std dev: ' + str(np.std(homo_chimp_rec)))
        
        bins = []
        start = 0
        bin_size = 0.1
        max_x_range = 5
        num_bins = int(max_x_range / bin_size)
        for i in range(0, num_bins + 1):
            bins.append(start)
            start = start + bin_size
            
        plt.hist(all_chimp_rec, bins, density = True, stacked = True, fill = False, alpha = 0.6, ec = 'blue')
        plt.hist(homo_chimp_rec, bins, density = True, stacked = True, fill = False, alpha = 0.5, ec = 'green')
        plt.gca().set(xlabel = 'Chimp Recombination', ylabel = 'Probability \n (÷ by bin width of ' + str(bin_size) + ')', ylim = [0, 3]);
        plt.savefig(chimp_rec_graph)
        plt.show()
        
    if i == 1:
        print('Vervet-ALL')
        print('min: ' + str(min(all_vervet_rec)))
        print('max: ' + str(max(all_vervet_rec)))
        print('median: ' + str(np.median(all_vervet_rec)))
        print('mean: ' + str(np.mean(all_vervet_rec)))
        print('std dev: ' + str(np.std(all_vervet_rec)))
        print('Vervet-HOMO')
        print('min: ' + str(min(homo_vervet_rec)))
        print('max: ' + str(max(homo_vervet_rec)))
        print('median: ' + str(np.median(homo_vervet_rec)))
        print('mean: ' + str(np.mean(homo_vervet_rec)))
        print('std dev: ' + str(np.std(homo_vervet_rec)))
        
        bins = []
        start = 0
        bin_size = 0.1
        max_x_range = 5
        num_bins = int(max_x_range / bin_size)
        for i in range(0, num_bins + 1):
            bins.append(start)
            start = start + bin_size
            
        plt.hist(all_vervet_rec, bins, density = True, stacked = True, fill = False, alpha = 0.6, ec = 'blue')
        plt.hist(homo_vervet_rec, bins, density = True, stacked = True, fill = False, alpha = 0.5, ec = 'green')
        plt.gca().set(xlabel = 'Vervet Recombination', ylabel = 'Probability \n (÷ by bin width of ' + str(bin_size) + ')', ylim = [0, 3]);
        plt.savefig(vervet_rec_graph)
        plt.show()
        
    if i == 2:
        print('Human-HOMO')
        print('min: ' + str(min(homo_human_rec)))
        print('max: ' + str(max(homo_human_rec)))
        print('median: ' + str(np.median(homo_human_rec)))
        print('mean: ' + str(np.mean(homo_human_rec)))
        print('std dev: ' + str(np.std(homo_human_rec)))
        
        bins = []
        start = 0
        bin_size = 0.1
        max_x_range = 10
        num_bins = int(max_x_range / bin_size)
        for i in range(0, num_bins + 1):
            bins.append(start)
            start = start + bin_size
            
        plt.hist(homo_human_rec, bins, density = True, stacked = True, fill = False, alpha = 0.5, ec = 'green')
        plt.gca().set(xlabel = 'Human Recombination', ylabel = 'Probability \n (÷ by bin width of ' + str(bin_size) + ')');
        plt.savefig(human_rec_graph)
        plt.show()






