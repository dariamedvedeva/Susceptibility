# -*- coding: utf-8 -*-

#import base_neighbor_data as neig_vectors
import Functions as suscept_calc


#determine_data is in base_neighbor_data.py 

def calculate_chi_of_neighbors():
    suscept_calc.chi_first_neighbor()
    suscept_calc.chi_second_neighbor()
    
def calculate_chi_of_period():
    suscept_calc.calculate_chi_along_x()
    

calculate_chi_of_neighbors()
calculate_chi_of_period()

print("END")

