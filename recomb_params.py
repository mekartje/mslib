#Script using mslib to obtain a handful of summaries. returns a table # rows = # replicates
##1st argument -- path to ms or ms_recomb results file
##2nd argument -- outfile path

import mslib, sys

sim = mslib.simulation(sys.argv[1])

with open(sys.argv[2], 'w') as outf:
    #write header line
    outf.writelines('rep\thaplo_num\tr2_mean\tr2_var\thap_freqs\n')
    rep = 0
    for rep in range(sim.reps):
        haplo_freqs = mslib.seq_freqs(sim.replicate(rep))
        haplo_num = len(haplo_freqs)
        freq_line = ''
        for i in haplo_freqs.values():
            freq_line += str(i) + ','
        freq_line = freq_line[:-1]
        r2 = mslib.get_r2_meanVar(sim.replicate(rep))
        outf.writelines(str(rep) + '\t' + str(haplo_num) + '\t' + str(r2['mean']) + '\t' + str(r2['var']) + '\t' + freq_line + '\n')
        rep += 1
