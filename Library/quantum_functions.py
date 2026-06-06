import networkx as nx
import numpy as np



####################################
def sort_two_lists (list1, list2):
    return zip(*sorted(zip(list1, list2)))


#####

def measure_auc_entanglement (res):

    x, y = sort_two_lists (res[:, 0], res[:, 1]) 

    auc = 0.0
    old_x = 0.0
    old_y = 0.0
    for i in range(0, len(x)):
        dx =  x[i] - old_x
        dy = 0.5* (y[i] + old_y)
        auc += dy * dx
        old_x = x[i]
        old_y = y[i]
    dx = 1.0 - old_x
    #to avoid problems with networks with GCC < 1
    dy = 0.5* (1.0 + old_y)
    #dy = old_y
    auc += dy * dx
    
    return auc

#####


##################
def schmidt_to_concurrence (lam):
    c = 2 * np.sqrt(lam) * np.sqrt(1 - lam)
    if c > 1:
        c = 1
    if c < 0:
        c = 0
    return c

def concurrence_to_schmidt (c):
    lam = 0.5 + 0.5 * np.sqrt(1 - c * c)
    return lam
#################
def compute_score_path (nr_edges, weight, perc_type):
    score = 1.0
    if perc_type == 1:
        conc = schmidt_to_concurrence (weight)
        score = nr_edges * np.log(conc)
        score = np.exp(score)
        score = concurrence_to_schmidt (score)
    if perc_type == 2:
        prob =  nr_edges * np.log(2.0*(1-weight))
        score = np.exp(prob)
    return score

###############################
###############################
def select_random_pair (G):
    pair = np.random.choice(G.nodes(), size = 2, replace = False)
    source = pair[0]
    target = pair[1]
    if G.degree(source) > G.degree(target):
         target = pair[0]
         source = pair[1]
    return [str(source), str(target)] 
##############################
def find_all_shortest_paths (G, pair):
    source = pair[0]
    target = pair[1]
    try:
        paths = nx.all_shortest_paths(G, source, target, method='bellman-ford')
        return list(paths)
    except:
        return []
##############################
def select_random_path (paths):
    p = np.random.randint(len(paths))
    return paths[p]
##############################
def remove_edges_from_path (G, path):
    removed_edges = []
    for i in range(0, len(path)-1):
        n = path[i]
        m = path[i+1]
        G.remove_edge(n,m)
        removed_edges.append([n,m])
    return removed_edges
##############################
def find_best_neighbor_refined (G, focal_node, other_node, max_distance):

    ##standard version
    if max_distance == 0:
        return focal_node

    ##improved version, max distance == 1
    if max_distance == 1:
        control = -1
        best_node = focal_node
        max_degree = G.degree(focal_node)
        for m in G.neighbors(focal_node):
            if G.degree(m)>max_degree and control < 0:
                max_degree = G.degree(m)
                best_node = m
            if m == other_node:
                best_node = m
                control = 1
                
    ##improved version, max distance == 2
    if max_distance == 2:
        control = -1
        best_node = focal_node
        max_degree = G.degree(focal_node)
        for m in G.neighbors(focal_node):
            if G.degree(m)>max_degree and control < 0:
                max_degree = G.degree(m)
                best_node = m
            if m == other_node:
                best_node = m
                control = 1
            for v in G.neighbors(m):
                if G.degree(v)>max_degree and control < 0:
                    max_degree = G.degree(v)
                    best_node = v
                if v == other_node:
                    best_node = v
                    control = 1            
    return best_node
###############################
def pair_entanglement_refined (G, pair, weight, perc_type, max_distance):

    removed_edges, scores = {}, {}
    nodes = []

    source = pair[0]
    target = pair[1]

    best_source = find_best_neighbor_refined (G, source, target, max_distance)
    best_target = find_best_neighbor_refined (G, target, source, max_distance)

    #print('# source ', source, ' best neig source ', best_source)
    #print('# target ', target, ' best neig target ', best_target)

    
    #case 0
    if (best_source == target) or (best_target == source):
        #print('#case 0')
        removed_edges[1], scores[1] = pair_entanglement (G, [source, target], weight, perc_type)
        nodes = [source, target]
    else:
        #case 1
        if (source != best_source) and (best_source != best_target) and (best_target != target):
            #print('#case 1')
            removed_edges[1], scores[1] = pair_entanglement (G, [source, best_source], weight, perc_type)
            removed_edges[3], scores[3] = pair_entanglement (G, [best_target, target], weight, perc_type)
            removed_edges[2], scores[2] = pair_entanglement (G, [best_source, best_target], weight, perc_type)
            nodes = [source, best_source, best_target, target]
        #case 2
        if (source != best_source) and (best_source != best_target) and (best_target == target):
            #print('#case 2')
            removed_edges[1], scores[1] = pair_entanglement (G, [source, best_source], weight, perc_type)
            removed_edges[2], scores[2] = pair_entanglement (G, [best_source, target], weight, perc_type)
            nodes = [source, best_source, target]
        #case 3
        if (source != best_source) and (best_source == best_target) and (best_target != target):
            #print('#case 3')
            removed_edges[1], scores[1] = pair_entanglement (G, [source, best_source], weight, perc_type)
            removed_edges[2], scores[2] = pair_entanglement (G, [best_source, target], weight, perc_type)
            nodes = [source, best_source, target]
        #case 4
        if (source != best_source) and (best_source == best_target) and (best_target == target):
            #print('#case 4')
            removed_edges[1], scores[1] = pair_entanglement (G, [source, target], weight, perc_type)
            nodes = [source, target]
        #case 5
        if (source == best_source) and (best_source != best_target) and (best_target != target):
            #print('#case 5')
            removed_edges[1], scores[1] = pair_entanglement (G, [source, best_target], weight, perc_type)
            removed_edges[2], scores[2] = pair_entanglement (G, [best_target, target], weight, perc_type)
            nodes = [source, best_target, target]
        #case 6
        if (source == best_source) and (best_source != best_target) and (best_target == target):
            #print('#case 6')
            removed_edges[1], scores[1] = pair_entanglement (G, [source, target], weight, perc_type)
            nodes = [source, target]
        #case 7
        if (source == best_source) and (best_source == best_target) and (best_target != target):
            #print('#case 7')
            removed_edges[1], scores[1] = pair_entanglement (G, [source, target], weight, perc_type)
            nodes = [source, target]
        #case 8 #impossible, source and target would be the same nodes
        if (source == best_source) and (best_source == best_target) and (best_target == target):
            #print('#case 8')
            removed_edges[1], scores[1] = [], []
            nodes = []

    return removed_edges, scores, nodes
####################################
def pair_entanglement (G, pair, weight, perc_type):
    current_score = 1.0
    removed_edges = {}
    scores = {}
    t = -1
    control = -1
    while control < 0:
        paths = find_all_shortest_paths(G, pair)
        if len(paths) == 0:
            control = 1
        if len(paths)>0:
            t += 1
            selected_path = select_random_path (paths)
            #print (selected_path, len(selected_path))
            removed_edges[t] = remove_edges_from_path (G, selected_path)
            scores[t] = compute_score_path (len(selected_path)-1, weight, perc_type)
            if perc_type == 1:
                current_score = current_score * scores[t]
                if current_score < 0.5:
                    control = 1
    return removed_edges, scores
#############################
def restore_edges (G, removed_edges):
    for t in removed_edges:
        for e in removed_edges[t]:
            n = e[0]
            m = e[1]
            G.add_edge(n, m)
#############################
def restore_edges_refined (G, removed_edges):
    for q in removed_edges:
        restore_edges (G, removed_edges[q])
#############################
def quantum_percolation_refined (G, P, weight, perc_type, max_distance):
    results = {}
    for p in range(0, P):
        pair = select_random_pair (G)
        removed_edges, scores, nodes = pair_entanglement_refined (G, pair, weight, perc_type, max_distance)
        restore_edges_refined (G, removed_edges)
        results[p] = [ [pair[0], pair[1]] , removed_edges, scores, nodes]
    return results
############################
def estimate_entanglement_result_refined (results, perc_type):

    #print('\n\n#estimate result\n')
    
    source = results[0][0]
    target = results[0][1]


    #print ('#source ', source)
    #print ('#target ', target)

    max_index = len(results[3])-1
    #print ('# max_index ', max_index)
    #print (results[3])

    
    ##check if edges form a proper path
    control = 1
    total_edges = 0
    ##################
    for i in range(1, max_index+1, 1):
        for q in results[1][i]:
            n = results[1][i][q][0][0]
            m = results[1][i][q][-1][1]
            total_edges +=  len(results[1][i][q])
            if n != results[3][i-1]:
                control = -1
                #print ('#err a ',n)
            if m != results[3][i]:
                control = -1
                #print ('#err b ',m)
    #################


    ent = 0.0
    ###################
    if control > 0:
        score_steps= []
        for i in results[2]:
            #print (i, results[2][i])
            tmp = 1.0
            for q in results[2][i]:
                if perc_type == 1:
                    tmp *= results[2][i][q]
                if perc_type == 2:
                    tmp *= 1.0 - results[2][i][q]
            if perc_type == 1:
                score_steps.append(tmp)
            if perc_type == 2:
                score_steps.append(1-tmp)

        for i in range(0, len(score_steps)):
             if perc_type == 1:
                 if score_steps[i] < 0.5:
                     score_steps[i] = 0.5
            
        tmp = 1.0
        for s in score_steps:
            if perc_type == 1:
                if s < 0.5:
                    s = 0.5
                c = schmidt_to_concurrence (s)
                tmp *= c
            if perc_type == 2:
                tmp *= s
        if perc_type == 1:
            tmp = concurrence_to_schmidt (tmp)
            ent = 2*(1-tmp)
        if perc_type == 2:
            ent = tmp
    ####################        

    if control < 0:
        print('# no path found\n')

    return ent, total_edges
############################
def summarize_results_refined (results, perc_type):
    total = av = avc = ave = connected = 0.0
    for p in results:

        total += 1.0
        ent = 0.0
        ####

        ent, edges = estimate_entanglement_result_refined (results[p], perc_type)

        if ent > 0:
            ave += edges
            connected += 1
            av += ent
            avc += schmidt_to_concurrence (1-0.5*ent)
            
        ####

    av /= total
    avc /= total
    ave /= total
    connected /= total

    return [av, avc, ave, connected]
###########################








###################################
##Equal-weight case dedicated functions
##################################


##############################
def remove_all_possible_paths(G, pair):
    removed_edges = {}
    t = -1
    control = -1
    while control < 0:
        paths = find_all_shortest_paths(G, pair)
        if len(paths) == 0:
            control = 1
        if len(paths)>0:
            t += 1
            selected_path = select_random_path (paths)
            removed_edges[t] = remove_edges_from_path (G, selected_path)
    return removed_edges
#############################
#############################
def restore_edges_and_create_summary (G, pair, removed_edges):
    data = [len(removed_edges)]
    for t in removed_edges:
        data.append(len(removed_edges[t])+1) ##number of nodes in the path
        for e in removed_edges[t]:
            n = e[0]
            m = e[1]
            G.add_edge(n, m)
    degree = [G.degree(pair[0]), G.degree(pair[1])]
    summary = degree + data
    return summary
############################
def quantum_percolation_equal_weights (G, P):
    results = {}
    for p in range(0, P):
        pair = select_random_pair (G)
        removed_edges = remove_all_possible_paths(G, pair)  
        data = restore_edges_and_create_summary (G, pair, removed_edges)
        results[p] = data
    return results
############################
def print_to_file_equal_weights (results, outfile):
    with open(outfile, 'w') as f:
        for p in results:
            print(*results[p], file=f)
############################
def generate_diagram (data, perc_model, e_min, e_max, bins):

    res = []
    
    e_val = np.linspace(e_min, e_max, bins)
    for e in e_val:
        lam = 1.0 - 0.5 *e
        av = compute_entanglement (data, perc_model, lam)
        tmp = [e, av]
        res.append(tmp)
    return np.array(res)
#################
def compute_entanglement (data, perc_model, lam):
    av = 0.0
    total = 0.0
    for i in range(0, len(data)):
        score = 1.0
        total = total + 1
        for j in range(3, len(data[i])):
            nr_edges = data[i][j]-1
            val = compute_score_path (nr_edges, lam, perc_model)
            if perc_model == 1:
                score = score * val
            if perc_model == 2:
                score = score * (1-val)
        if perc_model == 1:
            if score < 0.5:
                score = 0.5
            ent = 2*(1-score)
        if perc_model == 2:
            ent = 1 - score
        av += ent
    return av/total
#################
