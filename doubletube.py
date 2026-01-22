# -*- coding: utf-8 -*-
"""
Created on Fri Jul  4 10:20:47 2025

@author: schwarz
"""

from TILMedia import Gas, Liquid
from statistics import mean
from heat_exchanger import Heat_exchanger
import numpy as np
mean = np.mean

class Doppelrohr(Heat_exchanger):
    
    def __init__(self, medium_1, medium_2):
        
        if any([medium_1 is None, medium_2 is None]):
            raise Exception("Medium_1 and medium_2 need to be defined at initialization")
            
        super().__init__(medium_1, medium_2)
        
    def calc_Re_1(self, v, X_char, eta=None, rho=None):
        
        if eta is None:
            self.medium_1.setState_pTxi(self.p_1, mean(self.T_1))
            eta = self.medium_1.eta # Pas'
        if rho is None:
            self.medium_1.setState_pTxi(self.p_1, mean(self.T_1))
            rho = self.medium_1.d
        
        self.Re_1 = rho*v*X_char/eta
        
        return self.Re_1
    
    def calc_Re_2(self, v, X_char, eta=None, rho=None):
        
        if eta is None:
            self.medium_2.setState_pTxi(self.p_2, mean(self.T_2))
            eta = self.medium_2.eta # Pas
        if rho is None:
            self.medium_2.setState_pTxi(self.p_2, mean(self.T_2))
            rho = self.medium_2.d # kg/m^3
        
        self.Re_2 = rho*v*X_char/eta
        
        return self.Re_2
    
    def calc_Nu_1(self, Re, d_i, L, Pr=None, formula='laminar'):
        
        if Pr is None:
            self.set_state_1()
            Pr = self.medium_1.Pr
            
        if formula == 'laminar':
            
            Nu_I = 3.66
            Nu_II = 1.615*(Re*Pr*d_i / L)**(1/3)
            self.Nu_1 =  (Nu_I + 0.7**3 + (Nu_II-0.7)**3)**(1/3) # See VDI-Wärmeatlas
            
        return self.Nu_1
    
    def calc_Pr_1(self, **kwargs):
        if len(kwargs) > 0:
            p = kwargs['p']
            T = kwargs['T']
        else:
            p = self.p_1
            T = mean(self.T_1)
            
        self.medium_1.setState_pTxi(p, T)
        
        self.Pr_1 = self.medium_1.Pr
        
        return self.Pr_1
    
    
    def calc_T_1_x(self, T_w, alpha_1=None, m_flow_1=None, c_p=None, A=None):
        
        if alpha_1 is None: alpha_1 = self.alpha_1
        if m_flow_1 is None: m_flow_1 = self.m_flow_1
        if c_p is None:
            self.medium_1.setState_pTxi(self.p_1, mean(self.T_1))
            c_p = self.medium_1.cp
        if A is None:
            A = self.A
            
        T_1_x = T_w-(T_w-self.T_1[0])*np.exp(-alpha_1 * A/(m_flow_1*c_p))
        
        return T_1_x
    
    
# Test functions

dr = Doppelrohr(Gas("DryAir"), Liquid("Water"))
dr.T_1 = [298.15,298.15]
dr.T_2 = [383.15,383.15]
dr.p_1 = 100000
dr.p_2 = 800000
def calc_T_end(L, T_w=383.15, d_i=0.0044, m_flow=1.9470702094005574e-06): # flow at 100 ml/min air
    
    Re = 2.77287
    Nu = dr.calc_Nu_1(Re, d_i, L)
    alpha = dr.calc_alpha_1(Nu, d_i)
    A = np.pi*d_i*L
    T_end = dr.calc_T_1_x(T_w, alpha, m_flow, A=A)
    return T_end

def find_L_min(T_w=383.15, tol=1, increment=0.001):
    L_min = 0
    T_end = 298.15
    while T_w - T_end > 1:
        L_min += increment
        T_end = calc_T_end(L_min, T_w)
        
    return L_min

def length_test(T_w=383.15):
    print(f'Length test with \nT_w={T_w} K\nd_i=0.0044\nm_flow=1.947e-6')
    for L in range(1,10):
        l = L/100
        T_end = calc_T_end(l, T_w=T_w)
        print(f'{l} m: {T_end}')


if __name__ == '__main__':
    
    length_test(300)