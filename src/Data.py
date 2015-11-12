# -*- coding: utf-8 -*-
import numpy as np
from math import pi, sqrt

def determine_data():
    global w0, W, J, lambda_r, kpt, sqrt32, period, qstep, the_num_of_neib
    w0 = 0.3861
    W = 0.0182703
    J = 0.0362
    kpt = 50
    
    read_X_imp()
    
    lambda_r = np.zeros(Iwmax, dtype = np.complex128)    
    for i in range(Iwmax):
        lambda_r[i] = computeLambda(i)

    sqrt32      = sqrt(3.0)/2.0
    #   период интегрирования по ЗБ
    period      = 2.0*pi

    # Шаг по сетке
    qstep       = 1.0/kpt

    #кол-во соседей первого и второго порядков (каждого) равно 6ти
    the_num_of_neib = 6 


def define_vectors_of_inverse_lattice():
    # Векторы обратной решетки
    global b1, b2
    b1 = (period/(sqrt32))*np.array([0.0,1.0], dtype=np.double)
    b2 = (period/(sqrt32))*np.array([sqrt32, 0.5], dtype=np.double)


def define_coordinates_of_neighbors():
    n1 = [None] * the_num_of_neib
    n2 = [None] * the_num_of_neib
    
    n1[0] = np.array([1.0,0.0], dtype=np.double)
    n1[1] = np.array([-1.0,0.0], dtype=np.double)
    n1[2] = np.array([ 0.5,sqrt32], dtype=np.double)
    n1[3] = np.array([-0.5,sqrt32], dtype=np.double)
    n1[4] = np.array([ 0.5,-sqrt32], dtype=np.double)
    n1[5] = np.array([-0.5,-sqrt32], dtype=np.double)
    
    n2[0] = np.array([1.5,sqrt32], dtype=np.double)
    n2[1] = np.array([1.5,-sqrt32], dtype=np.double)
    n2[2] = np.array([0.0, sqrt(3.0)], dtype=np.double)
    n2[3] = np.array([0.0,-sqrt(3.0)], dtype=np.double)
    n2[4] = np.array([-1.5,sqrt32], dtype=np.double)
    n2[5] = np.array([-1.5,-sqrt32], dtype=np.double)
    global neig_crds
    neig_crds = [None] * 2
    
    neig_crds[0] = n1
    neig_crds[1] = n2
    
def read_X_imp():
    global Iwmax, omega, Chi_imp
    f = open("chi_imp_r.dat")
    lines = f.readlines()    
    Iwmax = len(lines)    
    
    omega = np.zeros(Iwmax, dtype = np.complex128)
    Chi_imp = np.zeros(Iwmax, dtype = np.complex128)
    
    j = 0
    for i in lines:
        line = i.split()
        omega[j] = np.float64(line[0]) * np.complex128(1j)
        Chi_imp[j] = np.float64(line[1]) + np.float64(line[2]) * np.complex128(1j)
        j += 1    
    
    f.close()    
    
def for_test():
    print(Chi_imp)
    
def computeLambda(i):
    return 2 * w0 * W * W/(omega[i]*omega[i] - w0*w0)

def get_kpt():
    return kpt

def get_omega():
    return omega

def get_chi_imp():
    return Chi_imp

def get_Iwmax():
    return Iwmax

def get_J():
    return J

def get_lambda():
    return lambda_r
    
def get_b1_b2():
    return(b1, b2)

def get_neig_crds():
    return neig_crds

def get_qstep():
    return qstep

def get_the_num_of_neib():
    return the_num_of_neib
