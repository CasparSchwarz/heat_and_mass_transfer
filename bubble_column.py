# -*- coding: utf-8 -*-
"""
Created on Thu Jan 15 11:53:00 2026

@author: schwarz
"""

from statistics import mean
from heat_exchanger import Heat_exchanger
import numpy as np
from icecream import ic

class Bubble_Column(Heat_exchanger):
    
    def __init__(self, medium):
        super().__init__(self, medium)
        
     
        
    def calc_a(self, rho_f, D, v_f, e_g, sigma, g=9.81):
        
        a = 1/(3*D) * ((g* D**2 * rho_f) / sigma)**(1/8) * ((g * D**3)/ v_f**2)**0.1 * e_g**1.13
        
     
    def calc_phi_z(self, phi_in, z, v_g, beta=None, a=None, **kwargs):
        '''Calculates the relative humidity of a gas at the height z inside
        a bubble column
        
        phi_in: relative humidity of inlet gas
        z:      height of bubble column in m
        v_g:    velocity of the bubbles in m/s
        beta:   mass transfer coefficient in m/s
        a:      volume specific area in m^2/m^3

        returns:
        phi_z:  relative humidity of the gas at z
        '''
        
        if beta is None:
            beta=self.beta
        if a is None:
            a = self.a
            
        phi_z = 1 - (1 - phi_in)*np.exp(- (beta*a)/v_g * z)
        
        return phi_z
        