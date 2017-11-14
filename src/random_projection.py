#!/usr/bin/python

from sklearn import random_projection
import numpy as np
from scipy.stats import ortho_group

def main():
    
    trans = random_projection.GaussianRandomProjection(n_components=26931)

    fit = trans.fit(np.random.rand(1,58037))

    # fit_matrix = ortho_group.rvs(58037)

    # fit_matrix = fit_matrix[:,1:36931]

    
    
    np.savetxt('/gpfs/group/su/lhgioia/map/data/recount/random_proj_rotation.csv',fit.components_,fmt='%.4e',delimiter=',')
    

    return

if __name__ == '__main__':
    main()
