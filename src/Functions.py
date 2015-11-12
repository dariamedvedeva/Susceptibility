# -*- coding: utf-8 -*-
import numpy as np
import Data as data

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

def chi_first_neighbor():
    order = 1
    chi_neighbor(order)
    print_chi_neighbor(order)

def chi_second_neighbor():
    order = 2
    chi_neighbor(order)
    print_chi_neighbor(order)    
     
def chi_neighbor(order):
    order = order - 1
    
    global Chi_p1, Chi_p2, Chi_loc
    
    Chi_q_local = np.zeros(Iwmax, dtype = np.complex128)
    Chi_p_q1   = np.zeros(Iwmax, dtype = np.complex128)
    Chi_p_q2   = np.zeros(Iwmax, dtype = np.complex128)
    
    #for j in range(the_num_of_neib):
    #    print(neig_crds[order][j])
            
    for qx in range(kpt):
        for qy in range(kpt):
            q = np.array([qx*qstep, qy*qstep])
            # Радиус вектор точки в обратном пространстве (r_k = (q_x b1_x + q_y b2_x, q_x b1_y + q_y b2_y) )
            r = np.array([q[0]*b1[0]+q[1]*b2[0],q[0]*b1[1]+q[1]*b2[1]])
            sum_of_exp = 0
            for j in range(the_num_of_neib):
                sum_of_exp = sum_of_exp + np.exp(np.dot(r, neig_crds[order][j])*1.0j)
                
            J_q = J * sum_of_exp
            
            Chi_q_local = Chi_q_local + 1.0/((1/Chi_imp) + lambda_r - J_q)
                       
            Chi_p_q1   = Chi_p_q1 + (np.exp(-np.dot(r,neig_crds[order][0])*1.0j))/((1/Chi_imp) + lambda_r - J_q)
                        
            Chi_p_q2   = Chi_p_q2 + (np.exp(-np.dot(r,neig_crds[order][2])*1.0j))/((1/Chi_imp) + lambda_r - J_q)
    
    Chi_loc    = Chi_q_local/(kpt*kpt)
    Chi_p1     = Chi_p_q1/(kpt*kpt)
    Chi_p2     = Chi_p_q1/(kpt*kpt)
    

def print_chi_neighbor(order):
    fout = open('chi_' + repr(order) +'_neib.dat', 'w')
    for ch in range(Iwmax):
        Somega = repr(omega[ch].imag) +"\t"
        #Sloc = repr(Chi_loc[ch].real) +"\t"+ repr(Chi_loc[ch].imag) +"\t"
        #S1 = repr(Chi_p1[ch].real) + "\t" + repr(Chi_p1[ch].imag) + "\t"
        S2 = repr(Chi_p2[ch].real) + "\t" + repr(Chi_p2[ch].imag) + "\n"
        outstring = Somega + S2 
        fout.write(outstring)
    fout.close()
    
Chi_qx_eq_0     = np.zeros(Iwmax, dtype = np.complex128)
Chi_qx_eq_pi4   = np.zeros(Iwmax, dtype = np.complex128)
Chi_qx_eq_pi2   = np.zeros(Iwmax, dtype = np.complex128)
Chi_qx_eq_3pi4  = np.zeros(Iwmax, dtype = np.complex128)
Chi_qx_eq_kpt   = np.zeros(Iwmax, dtype = np.complex128)
        

def calculate_chi_period_point(qx_in, qy_in):
    order = 0 # only with 1st neighbor    
    
    Chi_p   = np.zeros(Iwmax, dtype = np.complex128)        
    qx = qx_in
    qy = qy_in
    q = np.array([qx*qstep, qy*qstep])
    # Радиус вектор точки в обратном пространстве (r_k = (q_x b1_x + q_y b2_x, q_x b1_y + q_y b2_y) )
    r = np.array([q[0]*b1[0]+q[1]*b2[0],q[0]*b1[1]+q[1]*b2[1]])
    
    sum_of_exp = 0
    for j in range(the_num_of_neib):
        sum_of_exp = sum_of_exp + np.exp(np.dot(r, neig_crds[order][j])*1.0j) 
                    
    J_q = J * sum_of_exp       
    Chi_p   = Chi_p + 1.0/((1/Chi_imp) + lambda_r - J_q)
    return Chi_p
 
 
def calculate_chi_along_x(): 
    qy = 0
    
    qx = 0
    Chi_qx_eq_0     = calculate_chi_period_point(qx, qy)
    
    qx = kpt/4
    Chi_qx_eq_pi4   = calculate_chi_period_point(qx, qy)
    
    qx = kpt/2
    Chi_qx_eq_pi2   = calculate_chi_period_point(qx, qy)
    
    qx = 3 * kpt/4
    Chi_qx_eq_3pi4  = calculate_chi_period_point(qx, qy)
    
    qx = kpt
    Chi_qx_eq_kpt   = calculate_chi_period_point(qx, qy)
    
    print_chi_period_along_one_axis(Chi_qx_eq_0, Chi_qx_eq_pi4, Chi_qx_eq_pi2, Chi_qx_eq_3pi4, Chi_qx_eq_kpt, "Ox")
    
    
    
def print_chi_period_along_one_axis(chi_00, chi_pi_4, chi_pi_2, chi_3_pi_4, chi_2_pi, axis):
    filename = "chi_period_along_" + axis + ".dat"

    fout = open(filename, 'w')
    for ch in range(Iwmax):
        Somega = repr(omega[ch].imag) +"\t"
        Sqx0 = repr(chi_00[ch].real) +"\t"+ repr(chi_00[ch].imag) +"\t"
        Sqxkpt4 = repr(chi_pi_4[ch].real) +"\t"+ repr(chi_pi_4[ch].imag) +"\t"
        Sqxkpt2 = repr(chi_pi_2[ch].real) +"\t"+ repr(chi_pi_2[ch].imag) +"\t"
        Sqx3kpt4 = repr(chi_3_pi_4[ch].real) +"\t"+ repr(chi_3_pi_4[ch].imag) +"\t"
        Sqxkpt = repr(chi_2_pi[ch].real) +"\t"+ repr(chi_2_pi[ch].imag) +"\n"
   
    outstring = Somega + Sqx0 + Sqxkpt4 + Sqxkpt2 + Sqx3kpt4 + Sqxkpt
    fout.write(outstring)