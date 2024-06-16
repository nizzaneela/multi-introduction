#!/usr/bin/env python3
'''
Get coalescent times based on transmission network and time tree with labeled internal nodels 
'''
import treeswift
import os
import subprocess
import numpy as np
import pandas as pd
import math
import numpy.random as rng
import random

def tn_to_time_inf_dict(tn_path, tree, ascertained=True):
    # track the number of infected individuals at across each day in the simulation

    if tn_path.endswith('gz'):
        tn = gopen(tn_path)
    else:
        tn = open(tn_path)

    tree_label_dict = {(k.split('|')[1]):k  for k in list(tree.label_to_node().keys()) if k}    
    time_inf_dict = {} # specific people infected based on time
    current_inf_dict = {} # current number infected based on time
    total_inf_dict = {} # total number infected based on time
    for line in tn:
        if tn_path.endswith('gz'):
            row = line.decode().strip().split('\t')
        else:
            row = line.strip().split('\t')

        if row[1] in tree_label_dict:
            label = tree_label_dict[row[1]]

        if row[0] == 'None':
            # time_inf_dict[0] = [label]
            time_inf_dict[0] = []
            current_inf_dict[0] = 1
            total_inf_dict[0] = 1
            prior_day = 0
        else:
            day = math.ceil(float(row[2])*365)
            while day > prior_day + 1:
                prior_day += 1
                time_inf_dict[prior_day] = time_inf_dict[prior_day-1].copy()
                current_inf_dict[prior_day] = current_inf_dict[prior_day-1]
                total_inf_dict[prior_day] = total_inf_dict[prior_day-1]

            if day not in time_inf_dict:
                prior_day = day
                time_inf_dict[day] = time_inf_dict[day-1].copy()
                total_inf_dict[day] = total_inf_dict[day-1]
                current_inf_dict[day] = current_inf_dict[day-1]
            if row[0] == row[1]:
                if row[1] in tree_label_dict:
                    time_inf_dict[day].remove(label)
                current_inf_dict[day] -= 1
            else:
                if row[1] in tree_label_dict:
                    time_inf_dict[day].append(label)
                current_inf_dict[day] += 1
                total_inf_dict[day] += 1
    return time_inf_dict, current_inf_dict, total_inf_dict


def coalescent_timing(time_inf_dict, current_inf_dict, total_inf_dict, tree, num_days=100):
    # determine the stable coalescence of the tree
    num = 0
    for n in tree.traverse_inorder(): # traverse the tree to label internal nodes
        if n.is_leaf() == False:
            label = 'internal_' + (5-len(str(num)))*'0' + str(num)
            n.label = label
            num += 1

    distance_from_zero_dict = {} # record the distance from the root at each node
    for n in tree.distances_from_root():
        node_label = n[0].label
        distance_from_zero_dict[node_label] = n[1]

    times = list(time_inf_dict.keys()) # the number of days in the sim
    used_times = []
    heights = []
    total_inf = []
    current_inf = []
    current_samples = []
    current_labels = []

    current_time_index = 0
    for index, time in enumerate(times):
        if time > num_days: # sometimes gemf goes past the limit but we don't always know when
            break
        labels = time_inf_dict[time] # currently circulating infections
        if time == 0:
            used_times.append(time) # record time
            heights.append(time) # height of tree is 0 at time 0 
            current_time_index += 1 
            total_inf.append(total_inf_dict[time]) 
            current_inf.append(current_inf_dict[time]) 
            current_samples.append(len(labels))
            current_labels.append(labels) # currently circulating individuals
            continue
        if len(labels) == 0: # height of tree is 0 if there are no currently infected individuals
            heights.append(time)
        elif len(labels) == 1: # height of the tree is the current time if only one infected individual
            heights.append(min(time/365, float(labels[0].split('|')[-1])))
        else: # get height of the subtree; the tMRCA
            subtree = tree.extract_tree_with(labels)
            height_subtree = distance_from_zero_dict[subtree.root.label]
            heights.append(height_subtree)
        total_inf.append(total_inf_dict[time])
        current_inf.append(current_inf_dict[time])
        current_samples.append(len(labels))
        current_labels.append(labels)
        used_times.append(time)
        # stop at 50,000 infections
        if total_inf_dict[time] >= 50000: # stop when we reach 50,000 infections, as described in SM
            break
    # get the latest tMRCA, and then walk back back day as long as the tMRCA of the next earlier day
    # is earlier than the ltaest tMRCA by less than one day
    # or if the next earlier has only one sample, i.e. day doesn't have a tMRCA
    end_tMRCA = heights[-1]
    day = -1
    while (0 <= (end_tMRCA - heights[day-1]) <= (1/365)) and current_samples[day-1] != 1:
        day -= 1
    # stable coalescence is the MRCA of that day
    stable_coalescence = tree.mrca(current_labels[day])
    # return the actual stable coalescence
    coalescent_timing_results = [used_times, heights, total_inf, current_inf, current_samples]
    return heights[day], stable_coalescence, coalescent_timing_results


def clade_size(tree):
    # store the size of the clade for each node
    for n in tree.traverse_postorder():
        n.clade_size = 0
        if n.is_leaf():
            n.clade_size += 1
            continue
        for c in n.children:
            n.clade_size += c.clade_size
                
                
def mutation_simulation(n, mutations = [], size=29903, factor = 1, max_mutations = 3):
    # simulate mutations along each branch
    if n.is_root():
        n.mutations = []
        for child in n.children:
            mutation_simulation(child, n.mutations.copy(), factor=factor)
    
    elif n.is_leaf():
        n.mutations = mutations
        edge_length = n.get_edge_length()*factor
        for x in range(rng.poisson(29903 * edge_length)):
            n.mutations.append(random.getrandbits(16))
    
    else:
        n.mutations = mutations
        edge_length = n.get_edge_length()*factor
        for x in range(rng.poisson(29903 * edge_length)):
            n.mutations.append(random.getrandbits(16))
            
        for child in n.children:
            mutation_simulation(child, n.mutations.copy(), factor=factor)

            
def get_mutation_clades(n, num_mutations, root_mutations=0, include_unmutated_leaves=False):
    '''
    Return all clades that have at least num_mutations, along with the remaining leaves (clade size = 1) that never reach that many mutations
    
    Note: to get basally divergent clades, set num_mutations = 1
    '''
    clade_sizes = []
    n.mutations = [x for x in n.mutations if x != 'None']
    num_mutations += root_mutations
    
    if n.is_root():
        if len(n.mutations) >= num_mutations: # a subtree might already have enough mutations -- hence the whole tree is the desired clade
            return [n.clade_size]
        
        for child in n.children:
            clade_sizes += get_mutation_clades(child, num_mutations, include_unmutated_leaves=include_unmutated_leaves)
        return clade_sizes
    
    elif n.is_leaf():
        if len(n.mutations) >= num_mutations:
            return [n.clade_size]
        else:
            if include_unmutated_leaves:
                return [-1] # making it -1 will make for easy identification later
            else:
                return []
    
    else:
        if len(n.mutations) >= num_mutations:
            return [n.clade_size]
        else:
            for child in n.children:
                clade_sizes += get_mutation_clades(child, num_mutations, include_unmutated_leaves=include_unmutated_leaves)
        return clade_sizes


def get_mutation_clade_nodes(n, num_mutations, root_mutations=0, include_unmutated_leaves=True):
    '''
    Return all MRCA nodes of clades that have at least num_mutations, along with the remaining leaves (clade size = 1) that never reach that many mutations
    
    Note: to get basally divergent clades, set num_mutations = 1

    Unmutated leaves are returned with a clade size of -1 to differentiate them from mutated leaves
    '''
    
    clade_mrcas = []
    n.mutations = [x for x in n.mutations if x != 'None']
    num_mutations += root_mutations
    
    if n.is_root():
        if len(n.mutations) >= num_mutations: # a subtree might already have enough mutations -- hence the whole tree is the desired clade
            return [n]
        
        for child in n.children:
            clade_mrcas += get_mutation_clade_nodes(child, num_mutations, include_unmutated_leaves)
        return clade_mrcas
    
    elif n.is_leaf():
        if len(n.mutations) >= num_mutations:
            return [n]
        else:
            if include_unmutated_leaves:
                return [n]
            else:
                return []
    
    else:
        if len(n.mutations) >= num_mutations:
            return [n]
        else:
            for child in n.children:
                clade_mrcas += get_mutation_clade_nodes(child, num_mutations, include_unmutated_leaves)
        return clade_mrcas


# main function
if __name__ == "__main__":
    # parse args
    from sys import stdin,stdout; from gzip import open as gopen; import argparse
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-tn', '--transmission_network', required=True, type=str, help="transmission network")
    parser.add_argument('-tt', '--time_tree', required=True, type=str, help="internally labeled time tree")
    parser.add_argument('-g', '--gemf', required=True, type=str, help="GEMF output")
    parser.add_argument('-n', '--num_days', required=False, type=int, default=100, help='Number of days in simulation')
    parser.add_argument('-s', '--sub_rate', required=True, type=float, help="substitution rate")
    parser.add_argument('-m', '--mutations', required=True, type=int, help="number of mutations for clade analysis")
    parser.add_argument('-id', '--sim_id', required=True, type=str, help="output directory")
    parser.add_argument('-seed', '--random_seed', required=True, type=int, help="RNG seed for reproducibility", default=None)
    args,unknown = parser.parse_known_args()

    # Set the RNG seed
    rng.seed(args.random_seed)

    # identify index case
    if args.gemf.endswith('.gz'):
        gemf = gopen(args.gemf)
    else:
        gemf = open(args.gemf)

    for index, line in enumerate(gemf):
        if args.gemf.endswith('gz'):
            t,rate,v,pre,post,num0,num1,num2,num3,num4,num5,num6,num7,num8,num9,lists = [i.strip() for i in line.decode().split()]
        else:
            t,rate,v,pre,post,num0,num1,num2,num3,num4,num5,num6,num7,num8,num9,lists = [i.strip() for i in line.split()]
        
        index_case = v
        break

    # get subtree that excludes index case
    tree = treeswift.read_tree_newick(args.time_tree)
    subtree_leaves = []
    for n in tree.traverse_leaves():
        if index_case != n.label.split('|')[1]:
            subtree_leaves.append(n.label)
    subtree = tree.extract_tree_with(subtree_leaves)
    subtree.suppress_unifurcations()

    # stable coalescence
    time_inf_dict, current_inf_dict, total_inf_dict = tn_to_time_inf_dict(args.transmission_network, subtree)
    tSC, stable_coalescence, coal_timing = coalescent_timing(time_inf_dict, current_inf_dict, total_inf_dict, subtree, args.num_days)

    # prepare for clade analysis; get the subtree with the stable coalescence (MRCA) root
    subtree_sc = tree.extract_subtree(stable_coalescence)
    subtree_sc.root.edge_length = 0
    subtree_sc.suppress_unifurcations()

    # determine clade sizes,
    clade_size(subtree_sc)
    for resample in range(1000):
        # simulate mutations
        mutation_simulation(subtree_sc.root, factor = args.sub_rate)

        # what are the sizes of 2-mutation clades who parents have no mutations?
        clades_2muts_0parentMuts_polytomy = []
        for n in subtree_sc.traverse_internal():
            if len(n.mutations) >= 2 and len(n.parent.mutations) == 0:
                clades_2muts_0parentMuts_polytomy.append([n.clade_size, sorted(get_mutation_clades(n, args.mutations, root_mutations=len(n.mutations), include_unmutated_leaves=True), reverse=True)])

        # How many 0-mutation tips, 1-or-more-mutation tips, and 1-or-more-mutation internal nodes are there?
        clades_polytomy_CC = []
        for n in subtree_sc.traverse_rootdistorder():
            if n[1].is_root():
                continue
            elif len(n[1].parent.mutations) > 0:
                continue
            elif n[1].is_leaf():
                clades_polytomy_CC.append([n[1].clade_size, [n[1].clade_size]])
            elif not n[1].is_leaf() and len(n[1].mutations) > 0:
                clades_polytomy_CC.append([n[1].clade_size, sorted(get_mutation_clades(n[1], args.mutations, root_mutations=len(n[1].mutations), include_unmutated_leaves=True), reverse=True)])

        #same again but exact
        # what are the sizes of 2-mutation clades who parents have no mutations?
        clades_2muts_0parentMuts_polytomy_exact = []
        for n in subtree_sc.traverse_internal():
            if len(n.mutations) == 2 and len(n.parent.mutations) == 0:
                clades_2muts_0parentMuts_polytomy_exact.append([n.clade_size, sorted(get_mutation_clades(n, args.mutations, root_mutations=len(n.mutations), include_unmutated_leaves=True), reverse=True)])

        # How many 0-mutation tips, 1-or-more-mutation tips, and 1-or-more-mutation internal nodes are there?
        clades_polytomy_CC_exact = []
        for n in subtree_sc.traverse_rootdistorder():
            if n[1].is_root():
                continue
            elif len(n[1].parent.mutations) > 0:
                continue
            elif n[1].is_leaf():
                clades_polytomy_CC_exact.append([n[1].clade_size, [n[1].clade_size]])
            elif len(n[1].mutations) == 1 and len(n[1].parent.mutations) == 0:
                clades_polytomy_CC_exact.append([n[1].clade_size, sorted(get_mutation_clades(n[1], args.mutations, root_mutations=len(n[1].mutations), include_unmutated_leaves=True), reverse=True)])


        # write resample results
        OUTPUT_CC_polytomy_clade = 'resamples/' + str(resample) + '/clade_analyses_CC/' + args.sim_id + '_one_mutation_subclades.txt'
        OUTPUT_linAB_polytomy_clade = 'resamples/' + str(resample) + '/clade_analyses_AB/' + args.sim_id + '_two_mutation_subclades.txt'
        OUTPUT_CC_polytomy_clade_exact = 'resamples/' + str(resample) + '/clade_analyses_CC_exact/' + args.sim_id + '_one_mutation_subclades.txt'
        OUTPUT_linAB_polytomy_clade_exact = 'resamples/' + str(resample) + '/clade_analyses_AB_exact/' + args.sim_id + '_two_mutation_subclades.txt'
        with open(OUTPUT_CC_polytomy_clade, 'w') as f:
            for clade in clades_polytomy_CC:
                f.write('%s %s\n' % (str(clade[0]), str(clade[1])))
        with open(OUTPUT_linAB_polytomy_clade, 'w') as f:
            for clade in clades_2muts_0parentMuts_polytomy:
                f.write('%s %s\n' % (str(clade[0]), str(clade[1])))
        with open(OUTPUT_CC_polytomy_clade_exact, 'w') as f:
            for clade in clades_polytomy_CC_exact:
                f.write('%s %s\n' % (str(clade[0]), str(clade[1])))
        with open(OUTPUT_linAB_polytomy_clade_exact, 'w') as f:
            for clade in clades_2muts_0parentMuts_polytomy_exact:
                f.write('%s %s\n' % (str(clade[0]), str(clade[1])))
        # print root mutations
        OUTPUT_root_muts = 'resamples/' + str(resample) + '/' + args.sim_id + '_root_mutations.txt'
        with open(OUTPUT_root_muts, 'w') as f:
            f.write(str(rng.poisson(tSC*0.00092*29903)))


    # simulation results
    OUTPUT_subtree = 'simulations/' + args.sim_id + '/corrected_subtree.newick'
    OUTPUT_coal = 'simulations/' + args.sim_id + '/corrected_coalData_parameterized.txt'
    subtree.write_tree_newick(OUTPUT_subtree)
    with open(OUTPUT_coal, 'w') as f:
        f.write('time\tcoalescence time\ttotal infected\tcurrently infected\tcurrent samples\n')
        for index, i in enumerate(coal_timing[0]):
            f.write('%s\t%s\t%s\t%s\t%s\n' % (coal_timing[0][index], coal_timing[1][index], coal_timing[2][index], coal_timing[3][index], coal_timing[4][index]))
