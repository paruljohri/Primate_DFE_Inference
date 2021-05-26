


#This is a Python script that graphs the correlations between the recombination rates calculated 
#for chimps and vervets using the weighted average method and the linear interpolation method.

#This script is currently designed to be run locally, and will require path changes on the user's part to the variables:
#'filtered_exon_file' with input file  "FILTERED_homo_exons_reference_all_info.txt"
#'chimp_vs_rec_graph' with output file "chimp_weighted_avg_vs_interpolation_recombination.pdf"
#'vervet_vs_rec_graph' with output file "vervet_weighted_avg_vs_interpolation_recombination.pdf"




import matplotlib.pyplot as plt




filtered_exon_file = 'C:/Users/Kellen/Downloads/primate_dfe_project/chimp_and_vervet_data/FILTERED_homo_exons_reference_all_info.txt'
with open(filtered_exon_file, 'r') as f:
    filtered_exon_data = f.readlines()[1 : ]

chimp_interpolated_rec_data = []
chimp_weighted_avg_rec_data = []
vervet_interpolated_rec_data = []
vervet_weighted_avg_rec_data = []

for i in filtered_exon_data:
    split = i.strip('\n').split('\t')
    chimp_interpolated_rec = float(split[45])
    chimp_weighted_avg_rec = float(split[46])
    vervet_interpolated_rec = float(split[72])
    vervet_weighted_avg_rec = float(split[73])

    chimp_interpolated_rec_data.append(chimp_interpolated_rec)
    chimp_weighted_avg_rec_data.append(chimp_weighted_avg_rec)
    vervet_interpolated_rec_data.append(vervet_interpolated_rec)
    vervet_weighted_avg_rec_data.append(vervet_weighted_avg_rec)




primary_folder = 'C:/Users/Kellen/Downloads/primate_dfe_project/chimp_and_vervet_data/raw_plots'

chimp_vs_rec_graph = primary_folder + '/chimp_weighted_avg_vs_interpolation_recombination.pdf'
vervet_vs_rec_graph = primary_folder + '/vervet_weighted_avg_vs_interpolation_recombination.pdf'

fig = plt.figure(figsize = (20, 11.3))
ax = fig.add_subplot(111)
plt.rc('xtick', labelsize = 40)
ax.xaxis.set_tick_params(length = 8, width = 3)
ax.xaxis.set_tick_params(length = 5, width = 2, which = 'minor')
plt.rc('ytick', labelsize = 40)
ax.yaxis.set_tick_params(length = 8, width = 3)
ax.yaxis.set_tick_params(length = 5, width = 2, which = 'minor')
plt.close()

fig = plt.figure(figsize = (20, 15))
ax = fig.add_subplot(111)
ax.scatter(chimp_weighted_avg_rec, chimp_interpolated_rec)
ax.set_xlim([0, 3])
ax.set_ylim([0, 3])
ax.set_xlabel('Chimp Weighted Avg Recombination (4Ner/kb)', fontsize = 40)
ax.set_ylabel('Chimp Interpolation Recombination (4Ner/kb)', fontsize = 40)
fig.savefig(chimp_vs_rec_graph)

fig = plt.figure(figsize = (20, 15))
ax = fig.add_subplot(111)
ax.scatter(vervet_weighted_avg_rec_data, vervet_interpolated_rec_data)
ax.set_xlim([0, 3])
ax.set_ylim([0, 3])
ax.set_xlabel('Vervet Weighted Avg Recombination (4Ner/kb)', fontsize = 40)
ax.set_ylabel('Vervet Interpolation Recombination (4Ner/kb)', fontsize = 40)
fig.savefig(vervet_vs_rec_graph)
    






