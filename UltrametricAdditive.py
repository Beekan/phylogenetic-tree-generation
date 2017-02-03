from compsci260lib import *
from timeit import itertools


def solve_ultrametric_additive():
    
    dist_1 = {"1,2" : 0.3, "1,3" : 0.6, "1,4" : 0.5,
			  "2,3" : 0.5, "2,4" : 0.6,
              "3,4" : 0.9}
    
	# Fill in D_2 here.
    dist_2 = {"1,2" : 0.8, "1,3" : 0.6, "1,4" : 0.8, "1,5" : 0.6,
              "2,3" : 0.8, "2,4" : 0.4, "2,5" : 0.8,
              "3,4" : 0.8, "3,5": 0.2,
              "4,5" : 0.8}
    
    sodm =  {"1,2" : 0.5, "1,3" : 0.5, "1,4" : 0.1,
              "2,3" : 0.3, "2,4" : 0.5,
              "3,4" : 0.5}

  
    if is_ultrametric(sodm):
        
        print "dist is ultrametric.\n"
    else:
        print "dist is not ultrametric.\n"
        
    
    if is_additive(sodm):
        
        print "dist is additive.\n"
    else:
        print "dist is not additive.\n"

def find_n(dist):
    """Finds the number of sequences in a distance matrix"""
    nums = set()
    for key in dist:
        temp = key.split(",")
        for item in temp:
            nums.add(int(item))
    
    return len(nums)
        
def is_ultrametric(dist):
    """Return True if the given distance metric is an ultrametric, False otherwise."""
    
    n = find_n(dist)
    
    for i in range(1,n+1):
        for j in range(i+1, n+1):
            for k in range(j+1, n+1):

                dij = find_dist(dist, i, j)
                dik = find_dist(dist, i, k)
                djk = find_dist(dist, j, k)
                
                found = False
                
                if dij == dik and dik > djk:
                    found = True
                
                if djk == dij and dij > dik:
                    found = True
                  
                if dik == djk and djk > dij:
                    found = True
                
                if not found:
                    return False
    return True

def find_dist(dist,x,y):
    """Returns distance between two points in the distance matrix."""
    key = str(x) + "," + str(y)
    return dist[key]

def is_additive(dist):
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
                    
                    if d1 == d2:
                        if d2 <= d3:
                            return False
                      
                    elif d2 == d3:
                        if d3 <= d1:
                            return False
                      
                    elif d1 == d3:
                        if d1 <= d2:
                            return False
    
    return True

if __name__ == '__main__':
    solve_ultrametric_additive()