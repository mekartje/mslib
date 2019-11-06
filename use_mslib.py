##Python script implementing mslib, prints some summaries of msms output to STDOUT
##1st argument -- .txt file containing raw msms output

import sys
import mslib

sim_obj = mslib.simulation(sys.argv[1])

nuc_divers = []
watt_theta = []
#win_fst = []
Zns = []

for rep in range(sim_obj.reps):
    print('Replicate ' + str(rep) + ' in progress.')
    nuc_divers.append(mslib.get_nuc_divers(sim_obj.replicate(rep)))
    watt_theta.append(mslib.get_watt_theta(sim_obj.replicate(rep), sim_obj.samplesize, sim_obj.segsites(rep)))
    #win_fst.append(mslib.get_fst_rep(sim_obj.replicate(rep), sim_obj.pop_sizes[0], sim_obj.pop_sizes[1]))
    Zns.append(mslib.Zns(sim_obj.replicate(rep)))

avg_nuc_divers = sum(nuc_divers)/len(nuc_divers)
avg_watt_theta = sum(watt_theta)/len(watt_theta)
#avg_fst = sum(win_fst)/len(win_fst)
avg_Zns = sum(Zns)/len(Zns)
print('Average Total Nucleotide Diversity: ' + str(avg_nuc_divers))
print("Average Watterson's Theta: " + str(avg_watt_theta))
print('Average Window Fst: ' + str(avg_fst))
print('Average Zns: ' + str(avg_Zns))
