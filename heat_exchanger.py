# -*- coding: utf-8 -*-
"""
Created on Fri Jul  4 10:03:55 2025

@author: schwarz
"""

from TILMedia import Gas, Liquid
import numpy as np
mean = np.mean

class Heat_exchanger:
    
    A = None
    P_1 = None
    P_2 = None
    R_1 = None
    R_2 = None
    NTU_1 = None
    NTU_2 = None
    
    m_flow_1 = None
    m_flow_2 = None
    
    alpha_1 = None
    alpha_2 = None
    k = None
    
    # Inlet temp and outlet temp
    T_1 = [None, None]
    T_2 = [None, None]
    
    p_1 = None # Pa
    p_2 = None # Pa
    
    medium_1 = None
    medium_2 = None
    
    def __init__(self, medium_1, medium_2):
        
        self.medium_1 = medium_1
        self.medium_2 = medium_2
        
    # Setter functions
    
    def set_state_1(self, p=None, T=None):
        
        if self.medium_1 is None:
            return Exception("No medium defined for medium_1")
        
        if p is None and self.p_1 is None: 
            self.p_1 = float(input("Set pressure  [Pa] for medium_1: "))
            p = self.p_1
        elif p is None and self.p_1 is not None:
            p = self.p_1
        
            
        if T is None and self.T_1[0] is None: 
            self.T_1[0] = float(input("Set inlet temperature [K] for medium_1: "))
            self.T_1[1] = float(input("Set outlet temperature [K] for medium_1: "))
            T = mean(self.T_1)
        elif T is None and self.T_1[0] is not None:
            T = mean(self.T_1)
                
        self.medium_1.setState_pTxi(p, T)
        
    def set_state_2(self, p=None, T=None):
        
        if self.medium_2 is None:
            return Exception("No medium defined for medium_2")
        
        if p is None and self.p_2 is None: 
            self.p_2 = float(input("Set pressure  [Pa] for medium_2: "))
            p = self.p_2
        elif p is None and self.p_2 is not None:
            p = self.p_2
            
        if T is None and self.T_2[0] is None: 
            self.T_2[0] = float(input("Set inlet temperature [K] for medium_2: "))
            self.T_2[1] = float(input("Set outlet temperature [K] for medium_2: "))
            T=mean(self.T_2)
        elif T is None and self.T_2[0] is not None:
            T = mean(self.T_2)      

        self.medium_2.setState_pTxi(p, T)
        
    def set_A(self, A):
        if A < 0:
            raise Exception(f'A is cannot be smaller than 0: {A}')
        self.A = A
        
        return A
    
    def set_P_1(self, P_1):
        
        self.P_1 = P_1
        
        return P_1
    
    def set_P_2(self, P_2):
        
        self.P_2 = P_2
        
        return P_2
    
    def set_R_1(self, R_1):
        
        self.R_1 = R_1
        
        return R_1
    
    def set_R_2(self, R_2):
        
        self.R_2 = R_2
        
        return R_2
    
    def set_NTU_1(self, NTU_1):
        
        self.NTU_1 = NTU_1
        
        return NTU_1
    
    def set_NTU_2(self, NTU_2):
        
        self.NTU_2 = NTU_2
        
        return NTU_2    
    
    # Calc functions
    def calc_P_1(self, **kwargs):
        if len(kwargs) > 0:
            T_1 = kwargs['T_1']
            T_2 = kwargs['T_2']
            
        else:
            T_1 = self.T_1
            T_2 = self.T_2
            
        self.P_1 = (T_1[0]-T_2[1])/(T_1[0]-T_1[1])
        
        return self.P_1
    
    def calc_P_2(self, **kwargs):
        if len(kwargs) > 0:
            T_1 = kwargs['T_1']
            T_2 = kwargs['T_2']
            
        else:
            T_1 = self.T_1
            T_2 = self.T_2
            
        self.P_2 = (T_2[1]-T_1[0])/(T_2[0]-T_2[1])
        
        return self.P_2
    
    def calc_NTU_1(self, **kwargs):
        if len(kwargs) > 0:
            k = kwargs['k']
            m_flow = kwargs['m_flow']
            A = kwargs['A']
            c_p = kwargs['c_p']
            
        else:
            self.medium_1.setState_pTxi(self.p_1, mean(self.T_1))
            k = self.k # W/m^2K
            m_flow = self.m_flow_1 # kg/s
            A = self.A # m^2
            c_p = self.medium_1.cp # J/kgK
            
        self.NTU_1 = (k*A)/(m_flow*c_p)
        
        return self.NTU_1
    
    def calc_NTU_2(self, **kwargs):
        if len(kwargs) > 0:
            k = kwargs['k']
            m_flow = kwargs['m_flow']
            A = kwargs['A']
            c_p = kwargs['c_p']
            
        else:
            self.medium_2.setState_pTxi(self.p_2, mean(self.T_2))
            k = self.k # W/m^2K
            m_flow = self.m_flow_2 # kg/s
            A = self.A # m^2
            c_p = self.medium_2.cp # J/kgK
            
        self.set_NTU_1((k*A)/(m_flow*c_p))
        
        return self.NTU_2
    
    def calc_R_1(self, **kwargs):
        if len(kwargs) > 0:
            m_flow_1 = kwargs['m_flow_1']
            m_flow_2 = kwargs['m_flow_2']
            c_p_1 = kwargs['c_p_1']
            c_p_2 = kwargs['c_p_2']
            
        else:
            self.medium_1.setState_pTxi(self.p_1, mean(self.T_1))
            self.medium_2.setState_pTxi(self.p_2, mean(self.T_2))
            m_flow_1 = self.m_flow_1
            m_flow_2 = self.m_flow_2
            c_p_1 = self.medium_1.cp
            c_p_2 = self.medium_2.cp
            
        self.R_1 = (m_flow_1*c_p_1)/(m_flow_2*c_p_2)
        
        return self.R_1
    
    def calc_R_2(self, **kwargs):
        
        self.R_2 = 1/self.calc_R_1(kwargs)
        
        return self.R_2
    
    
    def calc_alpha_1(self, Nu, L, lamb=None):
        
        if lamb is None:
            self.set_state_1()
            lamb = self.medium_1.lamb
            
        self.alpha_1 = Nu*lamb / L
        
        return self.alpha_1
    
    def calc_alpha_2(self, Nu, L, lamb=None):
        
        if lamb is None:
            self.set_state_2()
            lamb = self.medium_2.lamb
            
        self.alpha_2 = Nu*lamb / L
        
        return self.alpha_2
    