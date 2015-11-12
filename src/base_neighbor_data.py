# -*- coding: utf-8 -*-

import numpy as np
import Data as data

from math import pi, sqrt

data.determine_data()
data.define_vectors_of_inverse_lattice()
data.define_coordinates_of_neighbors()
        
kpt         = data.get_kpt()
omega       = data.get_omega()
Chi_imp     = data.get_chi_imp() 
Iwmax       = data.get_Iwmax()
J           = data.get_J()
lambda_r    = data.get_lambda()
neig_crds   = data.get_neig_crds()
(b1, b2)    = data.get_b1_b2()
qstep       = data.get_qstep()
the_num_of_neib = data.get_the_num_of_neib()