from compsci260lib import *
from UltrametricAdditive import *

# Get the aligned DNA sequences:
fasta_dict = get_fasta_dict('CoV.aligned.fasta')

seq_names = {1:"AAL57308.1", 2:"AAK83356.1", 3:"NP_937950.1", 4: "AAL80031", 
             5:"YP_209233.1", 6:"NP_045300.1", 7:"NP_040831.1", 
             8:"NP_598310.1", 9:"BAA06805", 10:"AAP41037.1", 11:"AAV49723"}

def compute_dist(ind1, ind2):
    
    seq1_name = seq_names[ind1]
    seq2_name = seq_names[ind2]
    
    seq1 = fasta_dict[seq1_name]
    seq2 = fasta_dict[seq2_name]
    
    mismatches = 0.0
    total = 0.0
    for i in range(len(seq1)):
        
        if seq1[i] != "-" and seq2[i] != "-":
            total += 1
            if seq1[i] != seq2[i]:
                mismatches += 1
    
    return mismatches/total



def build_tree():
    # Compute the pairwise distances:
    dist = {}
    for i in range(1, len(seq_names.keys())):
        for j in range(i+1, len(seq_names.keys()) + 1):
            # Write the code to complete the compute_dist function
            dist[str(i) +"," + str(j)] = compute_dist(i, j)
    
   
    # uncomment the following line to print out the distance dictionary
    print_dist(dist)
    
    # Check if the distance is ultrametric.
    if is_ultrametric(dist):
        print "\nThe distance is ultrametric."
    else:
        print "\nThe distance is not ultrametric."
    
    # Check if the distance is additive.
    if is_additive(dist):
        print "\nThe distance is additive.\n"
    else:
        print "\nThe distance is not additive.\n"
    
    
    ########### Implementation of the Neighbor Joining algorithm ###########
    # Keeping track of variable names:
    
    # dist              - dictionary containing the computed pair-wise distances between nodes
    # node_list         - array to contain the list of current nodes, while iterating to build the tree
    # avg_dist          - dictionary to contain the averaged distances for all current nodes
    # D                 - dictionary to contain the adjusted distances between nodes
    # newick            - dictionary to maintain the Newick representation of each leaf node or subtree
    
    ########### Initialization: ##############################################
    
    node_list = seq_names.keys() # L = the list of current nodes
    newick = {}
    for i in range(1, len(node_list)+1):
        newick[i] = "" + str(i)
    
    avg_dist = {}  # averaged distances (the r terms) for all current nodes
    
    for i in range(1, len(node_list) + 1):
        avg_dist[i] = 0       
        for j in range(1, len(node_list) + 1):
            if i != j:
                # if i < j, then  avg_dist[i] = dist["i,j"] 
                # else avg_dist[i] = dist["j,i"]
                avg_dist[i] += dist["%d,%d"%(i,j)] if i<j else dist["%d,%d"%(j,i)]
        avg_dist[i] = float(avg_dist[i]) / (len(node_list) -2)
    
    max_node = len(node_list) # the maximum key used for a node 
    print avg_dist
    print dist
    ############## Iteration to build the tree #######################
    
    # As long as the list contains more than 2 nodes, iterate 
    while len(node_list) > 2:
        
        ########### Begin your code ################
        
        # Compute the 'adjusted' distances between nodes using the original distances (from dist)
        # and averaged distances (from avg_dist)  
        
        # Let D be the dict to contain 'adjusted' distances
        D = {}
        
        
        # Loop through each pair of nodes as entered in the dist dict
        # Use the entries from the avg_dist and calculate entries for the D dict   
        
        for key in dist.keys():
            indices = key.split(",")
            D[key] = dist[key] - avg_dist[int(indices[0])] - avg_dist[int(indices[1])]
        
        print "D matrix"
        print_dist(D)
        # Pick the pair i,j in node_list for which adjusted distance D_ij is minimal. 
        # You may find the function two_min_in_dict helpful. 
        (i, j) = two_min_in_dict(D)
        print i,j
          
        
        # Define a new node k and set dist[m,k] for all nodes m in node_list as (dist[i,m] + dist[j,m] - dist[i,j]) / 2
        k = max_node + 1 # max_node had been earlier set to the largest key used for a node 
        max_node += 1
        m = 0
      
        for ind in range(len(node_list)):
            m = node_list[ind]

            if i == m or j == m:
                continue
                
            if i > m:
                pair1 = str(m) + "," + str(i)
            else:
                pair1 = str(i) + "," + str(m)
            
            if j > m:
                pair2 = str(m) + "," + str(j)
            else:
                pair2 = str(j) + "," + str(m)
                
            dist["%d,%d"% (m,k)] = (dist[pair1] + dist[pair2] - dist["%d,%d"%(i,j)]) / 2
       
        print "dist matrix"
        print_dist(dist)
        ########### End your code ################
        
        # Add the new node k to the Newick format representation of the tree with edges of lengths 
        # dik = (dist[i,j] + avg_dist[i] - avg_dist[j])/2
        # djk = dist[i,j]-d[i,k], 
        # joining k to i and j
        
        d_ik = (dist["%d,%d"%(i,j)] + avg_dist[i] - avg_dist[j]) /2
        d_jk = dist["%d,%d"%(i,j)] - d_ik        
        newick[k] = "(" + newick[i] + ":" + "%.4f"%(d_ik) + "," + newick[j] + ":" +"%.4f"%(d_jk) + ")"
        
        # Remove i and j from node_list and add k
        temp = []
        for idx in range(len(node_list)):
            if node_list[idx] != i and node_list[idx] != j:
                temp.append(node_list[idx])
        temp.append(k)
        
        node_list = list(temp)
        
        # Update the r terms
        if len(node_list) > 2:
            avg_dist[k] = 0
            for ind in range(len(node_list) - 1):
                m = node_list[ind]
                avg_dist[m] = avg_dist[m] * (len(node_list) - 1)
                avg_dist[m] -= dist["%d,%d"%(m,i)] if m < i else dist["%d,%d"%(i,m)]
                avg_dist[m] -= dist["%d,%d"%(m,j)] if m < j else dist["%d,%d"%(j,m)]
                avg_dist[m] += dist["%d,%d" % (m,k)]
                
                avg_dist[m] /= (len(node_list) - 2)
                avg_dist[k] += dist["%d,%d"%(m,k)]
                
            avg_dist[k] = avg_dist[k] / (len(node_list) - 2)
        
        # Remove any elements from the dict that contain nodes i or j
        delete_from_dict(dist, i, j)
        delete_from_dict(D, i, j)
        delete_from_dict(newick, i, j)
        delete_from_dict(avg_dist, i, j)
        
    print "The Newick representation for the phylogenetic tree is:\n(" + newick[node_list[0]] + ":" + "%.4f"%(dist.values()[0]) + "," + newick[node_list[1]] + ":0);\n"


def compute_dist(i, j):
    """Returns the distance between two sequences. The distance is computed
    as the ratio between the number of mismatches and the total number
    of matches or mismatches."""
    
    seq1 = fasta_dict[seq_names[i]]
    seq2 = fasta_dict[seq_names[j]]
    
    mismatch = 0
    matchORmismatch = 0
    
    for i in range(0, len(seq1)):
        ########### Begin your code ################ 
        # Fill in your code to count the number of matches and mismatches. Ignore gaps in either sequence.
        
        if seq1[i] != "-" and seq2[i] != "-":
            matchORmismatch += 1
            if seq1[i] != seq2[i]:
                mismatch += 1
  
        ########### End your code ################ 
    return float(mismatch)/matchORmismatch


def is_ultrametric(dist):
    """Return True if the given distance metric is an ultrametric, False otherwise."""
    ########### Begin - your code ################  
    # Fill in your code from Problem 2 UltrametricAdditive.py
    # Remember to use the is_almost_equal function instead of checking for equality
    """Return True if the given distance metric is an ultrametric, False otherwise."""
    
    n = find_n(dist)
    
    for i in range(1,n+1):
        for j in range(i+1, n+1):
            for k in range(j+1, n+1):

                dij = find_dist(dist, i, j)
                dik = find_dist(dist, i, k)
                djk = find_dist(dist, j, k)
                
                found = False
                
                if is_almost_equal(dij, dik) and dik > djk:
                    found = True
                
                if is_almost_equal(djk, dij) and dij > dik:
                    found = True
                  
                if is_almost_equal(dik, djk) and djk > dij:
                    found = True
                
                if not found:
                    return False
    return True
    
    
    
    ########### End - your code ################


def is_additive(dist):
    """Return True if the given distance metric is additive, False otherwise"""
    ########### Begin - your code ################ 
    # Fill in your code from Problem 2 UltrametricAdditive.py
    # Remember to use the is_almost_equal function instead of checking for equality
    
    """Return True if the given distance metric is additive, False otherwise."""
    n = find_n(dist)
    
    for i in range(1,n+1):
        for j in range(i+1, n+1):
            for k in range(j+1, n+1):
                for l in range(k+1, n+1):
                    
                    dij = find_dist(dist, i, j)
                    dik = find_dist(dist, i, k)
                    djk = find_dist(dist, j, k)
                    dil = find_dist(dist, i, l)
                    djl = find_dist(dist, j, l)
                    dkl = find_dist(dist, k, l)
                    
                    d1 = dij + dkl
                    d2 = dik + djl
                    d3 = dil + djk
                    
                    if is_almost_equal(d1, d2):
                        if d2 <= d3:
                            return False
                      
                    elif is_almost_equal(d2, d3):
                        if d3 <= d1:
                            return False
                      
                    elif is_almost_equal(d1, d3):
                        if d1 <= d2:
                            return False
    
    return True
    
    
    ########### End - your code ################   


def is_almost_equal(num_1, num_2):
    """Return true if the difference between the two parameters is negligible
    enough that the parameters can be considered equal."""
    threshold = 0.03
    if abs(num_1 - num_2) < threshold:
        return True
    else:
        return False


def delete_from_dict(dictionary, i, j):
    """Deletes the dict entries with keys that contain i or j.  Doesn't return
    anything; just changes the dictionary."""
    for k in dictionary.keys():
        ks = [int(_) for _ in str(k).split(',')]
        if i in ks or j in ks:
            del dictionary[k]


def min_in_dict(dictionary):
    """Returns the key associated with the minimum value in the dict."""
    import operator
    return min(dictionary.iteritems(), key=operator.itemgetter(1))[0]


def two_min_in_dict(dictionary):
    import operator
    sorted_dict = sorted(dictionary.iteritems(), key=operator.itemgetter(1))
    element = sorted_dict[0][0]  # get the first element of the tuple
    (i, j) = element.split(",")
    return (int(i), int(j))


def print_dist(dist):
    '''print the distance store in the given dict structure'''
    idx = [int(i.split(',')[0]) for i in dist.keys()]
    idx.extend([int(i.split(',')[1]) for i in dist.keys()])
    max_idx = max(idx)
    
    print ' \t',
    for col in range(2, max_idx + 1):
        print '%d\t' % col,
    print 
    
    for row in range(1, max_idx):
        print '%d\t' % row,
        for col in range(2, row+1):
            print ' \t',
        for col in range(row+1, max_idx + 1):
            if (str(row)+','+str(col)) in dist.keys():
                print '%.4f\t' % dist[str(row)+','+str(col)],
        print 


if __name__ == '__main__':
    


    build_tree()
    lst = [(3,4)]
    d = {'team':[(4,5)]}
    d['team'] += lst
    lst = [(8,9)]
    d['team'] += lst
    print d.values()