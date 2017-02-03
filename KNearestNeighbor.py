from compsci260lib import *
from math import sqrt

def solve_k_nearest(k):
    """Find the k nearest neighbors to find a prognosis for the patients"""
    
    input_patients = read_training_data("gene_expression_training_set.txt")
    new_patients = read_test_data("gene_expression_test_set.txt")

    # check to make sure k is odd
    if k % 2 == 0:
        if k != len(input_patients):
            k += 1
        else:
            k -= 1
    
    patient_num = 1 
    
    for n in new_patients:
        data = []
        # make a list of distances to all neighbors
        for i in input_patients:
            dist = compute_dist(i['expression'], n['expression'])
            data.append((dist,i['class']))
        sorted_data = sorted(data)
        
        # index 0 will have the count of class R, index 1 will have the count of class N
        cts = [0,0]
        
        n['five closest'] = []
        for x in range(k):
            
            if x < 5:
                n['five closest'].append(sorted_data[x])
            if sorted_data[x][1] == 'R':
                cts[0] +=1
            else:
                cts[1] +=1
        if cts[0] > cts[1]:
            n['class'] = 'R'
        else:
            n['class'] = 'N'
        
        print "patient number", patient_num, n
        patient_num += 1 
        
def read_training_data(file_name):
    """Read the training gene expression data from a text file. Note: the
       patients in the training data are classified as "R" (responsive to
       treatment) or "N" (non-responsive to treatment).  For example,
       input_patients[0]["class"] = the class of the first patient (R or N)
       input_patients[0]["expression"][0] = the expression of the first 
       gene for the first patient."""
    
    f = open(file_name, "r")
    if f is None:
        sys.stderr("File " + file_name + " does not exist") 
    else:
        input_patients = [] 
        lines = f.readlines()
        for line in lines:
            line = line.strip()
            data = line.split()      # check if you are splitting on "\t" which is what you want to split on
            class_name = data.pop(0)
            float_data = [float(datum) for datum in data]   # convert to float data
            patient = {"class": class_name, "expression": float_data}
            input_patients.append(patient)
        f.close()
        return input_patients

def read_test_data(file_name):
    """Read the test gene expression data from a text file. Note: the
       patients in the test data are not classified."""
    
    f = open(file_name, "r")
    if f is None:
        sys.stderr("File " + file_name + " does not exist") 
    else:
        test_patients = [] 
        lines = f.readlines()
        for line in lines:
            line = line.strip()
            data = line.split()      # check if you are splitting on "\t" which is what you want to split on
            float_data = [float(datum) for datum in data]
            patient = {"class": "unknown", "expression": float_data}
            test_patients.append(patient)
        f.close()
        return test_patients


def compute_dist(tuple_1, tuple_2):
    """Return the Euclidean distance between two points in any number of dimensions."""
    
    if len(tuple_1) != len(tuple_2):
        sys.stderr("Cannot compute Euclidean distance between tuples of different sizes!")
    dist = 0
    for i in range(len(tuple_1)):
        dist += (tuple_1[i] - tuple_2[i]) * (tuple_1[i] - tuple_2[i])
    return sqrt(dist)


if __name__ == '__main__':
    solve_k_nearest(5)