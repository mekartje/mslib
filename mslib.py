import sys, math

"""
Parses msms output
"""

class simulation(object):
    def __init__(self, path):
        with open(path, 'r') as inf:
            headerline = inf.readline().rstrip()
        #additional switches
        recom = False
        subpop = False
        if '-r' in headerline:
            recom = True
        if '-I' in headerline:
            subpop = True
        if recom is True:
            self.recomb = float(headerline.split(' ')[headerline.split(' ').index('-r') + 1])
        if subpop is True:
            self.subpop = int(headerline.split(' ')[headerline.split(' ').index('-I') + 1])
            self.pop_sizes = []
            for pop in range(int(self.subpop)):
                self.pop_sizes.append(int(headerline.split(' ')[headerline.split(' ').index('-I') + 2 + pop]))
        headerline = headerline.split('-')
        self.path = path
        self.samplesize = int(headerline[0].split(' ')[1])
        self.reps = int(headerline[0].split(' ')[2])
        self.theta = float(headerline[1].split(' ')[1])

    """
    list of sequences (i.e., seqence length x sample size matrix)
    (don't code adaptive allele any differently
    """
    def replicate(self, number):
        rep_ls = []
        rep_count = -1
        with open(self.path, 'r') as inf:
            for line in inf:
                if line.startswith('//'):
                    rep_count += 1
                    if rep_count == number:
                        inf.readline()#read off segsites and position lines
                        inf.readline()
                        samp = 1
                        while samp <= self.samplesize:
                            rep_ls.append(inf.readline().rstrip())
                            samp += 1
                        break
        return(rep_ls)

    """
    retun segsites for a replicate, number
    """
    def segsites(self, number):
        rep_count = 0
        segs = 0
        with open(self.path, 'r') as inf:
            for line in inf:
                if line.startswith('//'):
                    if rep_count == number:
                        segs = int(inf.readline().rstrip().split(' ')[1])
                        break
                    rep_count += 1
            return(segs)

"""
function returning total nucleotide diveristy, or within populatino nucleotide diversity (PER-SITE)
arguments:
    seqls -- list of sequences as output by SIM_OBJ.replicate(number)
    scope -- 'total' for total nucleotide diversity (default). 'within' for average within populatino nucleotide diversity
    n_pop_1, n_pop_2 -- number of individuals in populations 1 and 2

    ***Make fast option w/ map/reduce architecture***
"""
def get_nuc_divers(seqls, scope = 'total', n_pop_1 = 0, n_pop_2 = 0):
    snp_heterozygosity = []
    snp_id = []
    #need to convert to float for python 2 compatibility
    n_pop_1 = float(n_pop_1)
    n_pop_2 = float(n_pop_2)
    #iterate through segregating sites
    for snp in range(len(seqls[0])):
        #snp_id.clear() have to remove clearing step for python 2.6.6 compatibility. make new list instead
        snp_id = []
        #iterage through samples, appending the current snp to a list
        #need to catch ValueError for letters mixed in w simulated sequence data
        for sample in range(len(seqls)):
            try:
                snp_id.append(float(seqls[sample][snp]))
            except ValueError: #for adaptive varaints (encoded by msms as > 1 in some cases for soft sweeps (i.e., when > 1 adaptive variants are present at the onset of selection)
                snp_id.append(float(1))
        if scope == 'total':
            #obtain allele frequencies and total heterozygosity
            p_tot = float(sum(snp_id))/len(snp_id)
            H_tot = (float(len(snp_id))/(float(len(snp_id)) - 1)) * (1 - (float(p_tot)**2 + (1-float(p_tot))**2))
            snp_heterozygosity.append(H_tot)
        elif scope == 'within':
            p_pop1 = sum(snp_id[:int(n_pop_1)])/n_pop_1
            p_pop2 = sum(snp_id[int(n_pop_1):])/n_pop_2
            H_pop1 = (n_pop_1 / (n_pop_1 - 1)) * (1 - (p_pop1**2 + (1-p_pop1)**2))
            H_pop2 = (n_pop_2 / (n_pop_2 - 1)) * (1 - (p_pop2**2 + (1-p_pop2)**2))
            H_win = (H_pop1 + H_pop2) / float(2)
            snp_heterozygosity.append(H_win)
    nuc_divers = sum(snp_heterozygosity)/len(snp_heterozygosity)
    return(nuc_divers)

"""
function returning Watterson's theta for a replicate, rep (not per-site)
arguments:
    rep: list of sequences as output by SIM_OBJ.replicate(number)
    sampsize: sample size, as output by SIM_OBJ.samplesize
    segsites: number of segregating sites, as output by SIM_OBJ.segsites
"""
def get_watt_theta(rep, sampsize, segsites):
    a = 0
    for i in range(sampsize - 1):
        a += 1/(i + 1)
    watt_theta = segsites/a
    return(watt_theta)

"""
function returning list of snp fsts
arguments:
    seqls -- list of sequences as output by SIM_OBJ.replicate(number)
    n_pop_1, n_pop_2 -- number of individuals in populations 1 and 2
"""
def get_snp_fst(seqls, n_pop_1, n_pop_2):
    #initialize lists for snp fst results and individual snps
    snp_fst = []
    snp_id = []
    #need to convert to float for python 2 compatibility
    n_pop_1 = float(n_pop_1)
    n_pop_2 = float(n_pop_2)
    #iterate through segregating sites
    for snp in range(len(seqls[0])):
        #snp_id.clear() have to remove clearing step for python 2.6.6 compatibility. make new list instead
        snp_id = []
        #iterage through samples, appending the current snp to a list
        #need to catch ValueError for letters mixed in w simulated sequence data
        for sample in range(len(seqls)):
            try:
                snp_id.append(float(seqls[sample][snp]))
            except ValueError: #for adaptive varaints (encoded by msms as > 1 in some cases for soft sweeps (i.e., when > 1 adaptive variants are present at the onset of selection)
                snp_id.append(float(1))
        #omit sites containing at least one 'NA', although this shouldn't be a issue
        if 'NA' in snp_id:
            snp_fst.append('NA')
            continue
        #omit sites with > 2 alleles
        elif sum(i > 1 for i in snp_id) > 0:
            snp_id.append(float(1))
            continue
        else:
            #obtain allele frequencies and within population/total heterozygosity
            p_tot = float(sum(snp_id))/len(snp_id)
            H_tot = (float(len(snp_id))/(float(len(snp_id)) - 1)) * (1 - (float(p_tot)**2 + (1-float(p_tot))**2))
            p_pop1 = sum(snp_id[:int(n_pop_1)])/n_pop_1
            p_pop2 = sum(snp_id[int(n_pop_1):])/n_pop_2
            H_pop1 = (n_pop_1 / (n_pop_1 - 1)) * (1 - (p_pop1**2 + (1-p_pop1)**2))
            H_pop2 = (n_pop_2 / (n_pop_2 - 1)) * (1 - (p_pop2**2 + (1-p_pop2)**2))
            H_win = (H_pop1 + H_pop2) / float(2)
            #compute fst estimate, append to list of snp fst's
            try:
                snp_fst.append((H_tot - H_win)/float(H_tot))
            except ZeroDivisionError:
                snp_fst.append('NA')
    #return snp fst list
    return(snp_fst)

"""
function returning max snp fst for a replicate (i.e., window)
arguments:
    rep: list of sequences, as output by SIM_OBJ.replicate(number)
    nx: population size for pop x
"""
def get_max_snp_fst(rep, n1, n2):
    fst_ls = get_snp_fst(rep, n1, n2)
    return(max(fst_ls))

"""
funciton to obtain average snp fst
given replicate from a simulation object
##need to work out population size inheritance! pop sizes must bet defined in function call currently -- not desireable
arguments:
    rep: list of sequences, as output by SIM_OBJ.replicate(number)
    nx: population size for pop x (number of chromosomes)
"""
def get_fst_rep(rep, n1, n2):
    fst_ls = get_snp_fst(rep, n1, n2)
    return(sum(fst_ls)/len(fst_ls))

"""
function returning Dxy
arguments:
    rep: list of sequenses as output by SIM_OBJ.replicate(number)
    n_pop_x: number of chromosomes sampled from  population x (perhaps passed as SIM_OBJ.pop_sizes[x-1]
"""
def get_dxy(rep, n_pop_1, n_pop_2):
    nsites = len(rep[0])
    npairs = n_pop_1 * n_pop_2
    correction = nsites * npairs
    diffs = 0
    #choose each pair of sequences
    for i in range(n_pop_1):
        for j in range(n_pop_2):
            j += n_pop_1
            for site in range(nsites): #will error out if many adaptive alleles of independent origin (NameError)
                if rep[i][site] != rep[j][site]:
                    diffs += 1
    return(diffs/correction)

"""
function returning Hudson, Slatkin, Madison 'Fst' (1992, i think)
arguments:
    rep: list ofsequences as output by SIM_OBJ.replicate(number)
    n_pop_x: population size forpopulation x (number of chromosomes)
"""
def get_hsm_fst(rep, n_pop_1, n_pop_2):
    pi_win = get_nuc_divers(rep, scope = 'within', n_pop_1 = n_pop_1, n_pop_2 = n_pop_2) #returns list w/ two elements, pi for pop1 and pi for pop2
    dxy = get_dxy(rep, n_pop_1, n_pop_2)
    hsm_fst = 1 - (pi_win / dxy)
    return(hsm_fst)

"""
function returning list of site frequencies
arguments:
    rep: list of sequences as output by SIM_OBJ.replicate(number)
    allele: 'derived' (default) or 'minor' (i.e., return derived or minor allele frequencies.
"""
def get_site_freqs(rep, allele = 'derived'):
    site_freq_ls = []
    snp_id = []
    for snp in range(len(rep[1])):
        snp_id.clear()
        for sample in range(len(rep)):
            try:
                snp_id.append(float(rep[sample][snp]))
            except ValueError:
                print(rep[sample])
        if allele == 'derived':
            site_freq_ls.append(sum(snp_id)/len(snp_id))
        elif allele == 'minor':
            if sum(snp_id) <= len(snp_id)/2: #if derived variant is minor
                site_freq_ls.append(sum(snp_id)/len(snp_id))
            elif sum(snp_id) > len(snp_id)/2:
                site_freq_ls.append((1 - sum(snp_id))/len(snp_id))
    return(site_freq_ls)

"""
fuction returning genotype frequencies between a pair of sites
get counts, rather than freqs, for FET implementation
arguments:
    rep: list of sequences asoutput by SIM_OBJ.replicate(number)
    posx: snp haplotypes to be counted
"""
def get_haplo_counts(rep, pos1, pos2): #positions are 0-based
    haplo_dict = {'anc1_anc2':0, 'anc1_der2':0, 'der1_anc2':0, 'der1_der2':0}
    for samp in range(len(rep)):
        if (int(rep[samp][pos1]) + int(rep[samp][pos2])) == 0:
            haplo_dict['anc1_anc2'] += 1
        elif (int(rep[samp][pos1]) + int(rep[samp][pos2])) == 1:
            if int(rep[samp][pos1]) == 1:
                haplo_dict['der1_anc2'] += 1
            else:
                haplo_dict['anc1_der2'] += 1
        elif (int(rep[samp][pos1]) + int(rep[samp][pos2])) == 2:
            haplo_dict['der1_der2'] += 1
    return(haplo_dict)

"""
function returning D given haplotype counts output by get_haplo_counts, list of site frequencies, and positions of sites
assumes derived allele frequencies (i.e., unfolded SFS)
arguments:
    haplo_count_dict: dictionary output by the get_haplo_counts() function
    site_freq_ls: list of site frequenies output by get_site_freqs() function
    posx: positions between which D will be computed
"""
def get_D(haplo_count_dict, site_freq_ls, pos1, pos2):
    obs_der_der = float(haplo_count_dict['der1_der2'])/sum(haplo_count_dict.values())
    exp_der_der = float(site_freq_ls[int(pos1)]) * float(site_freq_ls[int(pos2)])
    D = obs_der_der - exp_der_der
    return(D)

"""
function returning dictionary of r2 for a replicate
option to print dictionary to STDOUT
#use del rather than .remove(), remove returns nothing (i.e., 'NoneType'), and will remove the first value equivalent to the input value (i.e., not necessarily the same position as the reference list)
#some r2 >>1
#can make this a hell of a lot faster by just computing correlation coef for discrete data (i.e., cov/sqrt(blah)) and squaring. maybe.
arguments:
    rep: list of sequences as output by SIM_OBJ.replicate(number)
"""
def get_r2(rep):
    r2_dict = {}
    site_freq_ls = get_site_freqs(rep)
    for i in range(len(site_freq_ls)):
        for j in range(len(site_freq_ls)):
            if i == j:
                continue
            else:
                obs_haplo_counts = get_haplo_counts(rep, i, j)
                D = get_D(obs_haplo_counts, site_freq_ls, i, j)
                r2 = D**2/(site_freq_ls[i]*site_freq_ls[j]*(1 - site_freq_ls[i])*(1 - site_freq_ls[j]))
                r2_dict[(str(i),str(j))] = r2
    return(r2_dict)

"""
function returning Kelly's Zns (avg r2)
arguments:
    rep: list of sequences, as output by SIM_OBJ.repliate(number)
"""
def Zns(rep):
    r2_dict = get_r2(rep)
    Zns = sum(r2_dict.values())/len(r2_dict)
    return(Zns)

"""
distance matrix for all unique sequences in pooled sample
now, uses Euclidean dist. Modify to take other distance metrics (dictionary w/ sequence pairs as keys, dists as values)
arguments:
    rep: list ofsequences, as output by SIM_OBJ.replicate(number)
"""
def dist_matrix(rep):
    dist_dict = {}
    nsnps = len(rep[0])
    for i in rep:
        for j in rep:
            if i == j:
                continue
            elif (i,j) in dist_dict.keys():
                continue
            elif (j,i) in dist_dict.keys():
                continue
            else:
                dist_dict[(i,j)] = 0
                for snp in range(nsnps):
                    if i[snp] != j[snp]:
                        dist_dict[(i,j)] += 1
    return(dist_dict)

"""
sequence frequencies
returns dictionary with sequences as keys, their frequencies as values
arguments:
    rep: list of sequenes, as output by SIM_OBJ.replicate(number)
    scope: 'total' for pooled population sequence frequencies (default); 'within' for within-population sequence frequences.
        'within' option returns dict with 2 nested pop-specific dictionaries
    n_pop_x: population sizes for population x
"""
def seq_freqs(rep, scope = 'total', n_pop_1 = 0, n_pop_2 = 0):
    seq_freq_dict = {}
    if scope == 'total':
        for unique_seq in rep:
            if unique_seq in seq_freq_dict.keys():
                continue
            else:
                seq_freq_dict[unique_seq] = 0
                for seq in rep:
                    if seq == unique_seq:
                        seq_freq_dict[unique_seq] += 1
        for k,v in seq_freq_dict.items():
            seq_freq_dict[k] = float(v) / len(rep)
    elif scope == 'within':
        seq_freq_dict['pop_1'] = {}
        seq_freq_dict['pop_2'] = {}
        rep_1 = rep[:n_pop_1]
        rep_2 = rep[n_pop_1:]
        for unique_seq in rep_1:
            if unique_seq in seq_freq_dict['pop_1'].keys():
                continue
            else:
                seq_freq_dict['pop_1'][unique_seq] = 0
                for seq in rep_1:
                    if seq == unique_seq:
                        seq_freq_dict['pop_1'][unique_seq] += 1
        for unique_seq in rep_2:
            if unique_seq in seq_freq_dict['pop_2'].keys():
                continue
            else:
                seq_freq_dict['pop_2'][unique_seq] = 0
                for seq in rep_2:
                    if seq == unique_seq:
                        seq_freq_dict['pop_2'][unique_seq] += 1
        for k,v in seq_freq_dict['pop_1'].items():
            seq_freq_dict['pop_1'][k] = float(v) / n_pop_1
        for k,v in seq_freq_dict['pop_2'].items():
            seq_freq_dict['pop_2'][k] = float(v) / n_pop_2
    return(seq_freq_dict)

"""
number of unique haplotypes
returns the number of unique haplotypes in a dictionary output by seq_freqs
"""
def num_unique_haplo(seq_freq_dict):
    return(len(seq_freq_dict))

"""
average number of differences between sequences in a sample (total or between population)
if total is passed, give total seq freqs, if 'within', give pop specific
arguments:
    dist_matrix: dictionary of distances between unique sequences, as output by the function dist_matrix().
    seq_freq_dict: dictionary of sequence frequencies, as output by the function seq_freqs().
        total or within, corresponding to option passed for 'scope' argument.
    scope: 'total' for delta computed from pooled data (default). 'within' to compute average within-population delta.
"""
def get_delta(dist_matrix, seq_freq_dict, scope = 'total'):
    diffs = 0
    if scope == 'total':
        for seq_pair in dist_matrix.keys():
            freqi = seq_freq_dict[seq_pair[0]]
            freqj = seq_freq_dict[seq_pair[1]]
            diffs += freqi * freqj * dist_matrix[seq_pair]
    elif scope == 'within':
        pop_diff_sum = 0
        for pop in seq_freq_dict.keys():
            pop_diffs = 0
            for seq_pair in dist_matrix.keys():
                if seq_pair[0] not in seq_freq_dict[pop].keys() or seq_pair[1] not in seq_freq_dict[pop].keys():
                    continue
                freqi = seq_freq_dict[pop][seq_pair[0]]
                freqj = seq_freq_dict[pop][seq_pair[1]]
                pop_diffs += freqi * freqj * dist_matrix[seq_pair]
            pop_diff_sum += pop_diffs
        diffs = pop_diff_sum / 2
    return(diffs)

"""
'phi_st'
arguments:
    rep: list of sequences, as output by SIM_OBJ.replicate(number)
    n_pop_x: population size for population x
"""
def get_phist(rep, n_pop_1, n_pop_2):
    dist_mat = dist_matrix(rep)
    tot_freqs = seq_freqs(rep)
    pop_freqs = seq_freqs(rep, scope = 'within', n_pop_1 = n_pop_1, n_pop_2 = n_pop_2)
    tot_delta = get_delta(dist_mat, tot_freqs)
    win_delta = get_delta(dist_mat, pop_freqs, scope = 'within')
    phist = (tot_delta - win_delta) / tot_delta
    return(phist)

"""
Tajima's D
#can pull sampsize and segsites from SIM_OBJ
"""
def get_tajima_d(rep, segsites, sampsize):
    if segsites == 0:
        print('replicate contains no polymorphic sites')
    else:
        #watterson's theta here is unscaled, need to divide by number of sites for measure comprable with nucleotide diversity
        wat_theta = get_watt_theta(rep, sampsize = sampsize, segsites = segsites) / segsites
        nuc_divers = get_nuc_divers(rep)
        d = nuc_divers - wat_theta
        #standardization -- from Tajima (1989)
        n = sampsize
        S = segsites
        a1 = 0
        for i in range(n - 1):
            i = i + 1
            a1 += 1 / i
        a2 = 0
        for i in range(n - 1):
            i = i + 1
            a2 += 1 / i ** 2
        b1 = (n + 1) / (3 * (n - 1))
        b2 = (2 * (n**2 + n + 3)) / (9 * n * (n - 1))
        c1 = b1 - (1/a1)
        c2 = b2 - ((n + 2) / (a1 * n)) + (a2 / a1 ** 2)
        e1 = c1 / a1
        e2 = c2 / (a1 ** 2 + a2)
        var_d = (e1 * S + e2 * S * (S - 1)) ** 1/2
        taj_d = d / var_d
        return(d)

"""
Future modules:
    Fay and Wu's H
"""
