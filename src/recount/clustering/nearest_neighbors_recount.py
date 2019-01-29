#!/usr/bin/python

import numpy as np
import copy
import heapq
import sys
from scipy.spatial import distance
import csv

class node:
    def __init__(self):
        self.index =0
        self.threshold=0.0
        self.left=None
        self.right=None

def euclidian_distance(p1,p2):
    return np.sqrt(sum(np.square(p1-p2)))


def build_tree(lower,upper):

    global points
    global rowinds

    # lower = 0
    # upper = len(points)

    if(upper == lower):
        return None

    lower_node = node()
    lower_node.index = lower

    if upper - lower > 1:
        
        i = np.random.randint(lower,upper)
        
        tmp = copy.deepcopy(points[lower])
        points[lower] = points[i]
        points[i] = tmp

        tmp = copy.deepcopy(rowinds[lower])
        rowinds[lower] = rowinds[i]
        rowinds[i] = tmp

        # points[i], points[lower] = points[lower], points[i]

        median = int((upper+lower)/2)
        # dists = map(lambda x: np.linalg.norm(points[lower] - x),points[lower+1:upper])
        dists = distance.cdist(points[lower+1:upper,:],points[None,lower,:])

        dist_partition = np.argpartition(dists[:,0].tolist(),median-(lower+1))
                
        points[lower+1:upper] = points[lower+1+dist_partition]

        lower_node.threshold = np.linalg.norm(points[lower]-points[median])

        lower_node.index = lower

        lower_node.left = build_tree(lower+1,median)
        lower_node.right = build_tree(median,upper)

    return lower_node
        
def search(node, target, k):

    global tmpResults
    global tau

    # if node != None:
    #     print node.index
    # else:
    #     print -1
    # print tmpResults
    # print tau

    if node == None:
        return

    dist = np.linalg.norm(points[node.index]-target)

    if dist < tau:
        
        if len(tmpResults) == k:
            heapq.heappop(tmpResults)

        heapq.heappush(tmpResults,(-1*dist,node.index))
        if len(tmpResults) == k:
            tau = -1*heapq.nsmallest(1,tmpResults)[0][0]

    if node.left == None and node.right == None:
        return

    if dist < node.threshold:
        if dist - tau <= node.threshold:
            search(node.left,target,k)

        if dist + tau >= node.threshold:
            search(node.right,target,k)

    else:
        if dist + tau >= node.threshold:
            search(node.right,target,k)

        if dist - tau <= node.threshold:
            search(node.left,target,k)
            
if __name__ =='__main__':

    WORK_DIR = '/gpfs/group/su/lhgioia/map/'
    
    points = np.genfromtxt(WORK_DIR + 'results/recount/pca/recount_250_dim_noScaled_noProj_over50_noSingle.csv',delimiter=',')

    with open(WORK_DIR + 'results/recount/pca/recount_over50_noSingle_rownames.txt','r') as in_file:
        rownames = in_file.readlines()

    rowinds = map(lambda x: x,range(len(points)))

    # points = np.array([[0,0],[0,1],[1,0],[2,0]])
    # points = np.array([[2.5,3.3,2.1],[1.0,1.0,1.0],[4.2,4.0,5.0]])
    # points = np.random.rand(100,10)

    tree = build_tree(0,len(points))

    results = []

    for i in range(len(points)):
    # for i in range(5):
        if i % 1000 == 0:
            print i

        tmpResults = []
        tau = sys.float_info.max

        search(tree,points[i],100)

        nearest_neighbors = copy.deepcopy(tmpResults)

        # neighbor_inds = map(lambda x: rowinds[x[1]],nearest_neighbors) # corresponds to original row number
        neighbor_inds = map(lambda x: x[1],nearest_neighbors)

        if not i in neighbor_inds:
            print i
            print rownames[rowinds[i]].strip('\n')
            print rowinds[i]

            break
        # neighbor_names = map(lambda x: rownames[x].strip('\n'),neighbor_inds)
    
        # np.savetxt(WORK_DIR + 'results/recount/clustering/vptree_test.csv',nearest_neighbors,delimiter=',')

        # print neighbor_inds
        
        # print len(rowinds)
        # print len(rownames)

        # print neighbor_names


        entry = [rownames[rowinds[i]].strip('\n')]
        entry.extend(neighbor_inds)
        results.append(entry)

        # results[rownames[rowinds[i]].strip('\n')] = neighbor_inds

    with open(WORK_DIR + 'results/recount/clustering/recount_nn_k100.csv','w') as out_file:
        writer = csv.writer(out_file,delimiter=',')
        # for item in results.keys():
        for item in results:

            # row = [item]
            # row.extend(results[item])
            # writer.writerow(row)
            writer.writerow(item)
            # out_file.write('%s\n' % item)

    print 'done'



