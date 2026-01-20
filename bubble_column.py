# -*- coding: utf-8 -*-
"""
Created on Thu Jan 15 11:53:00 2026

@author: schwarz
"""

from TILMedia import Gas, Liquid
from statistics import mean
from heat_exchanger import Heat_exchanger
import numpy as np
from icecream import ic
from scipy.optimize import minimize

class Bubble_Column(Heat_exchanger):
    '''Bubble column class to calculate saturation of a liquid inside a gas
    medium_1: gaseous component
    medium_2: liquid component
    
    '''
    
    def __init__(self, medium_1, medium_2):
        super().__init__(medium_1, medium_2)
        
    
    def calc_e_g(self, D, sigma, v_g, C_1, rho_f=None, eta_f=None, g=None):
        '''Calculates the relative gas content in the fluid
        For air and water typicall between 0.022 - 0.238  [Ind. Eng. Chem. Process Des. Develop., Vol. 12, No. 1, 1973]
        
        D:          diameter of the bubble column in m
        sigma:      surface tension of the fluid in N/m
        v_g:        velocity of the gas in m/s
        C_1:        initial dissolved gas in kg/m^3
        rho_f:      density of the fluid in kg/m^3
        eta_f:      dynamic viscosity of the fluid in Pa*s
        g:          gravitational acceleration in m/s^2
        
        returns:
        e_g:        relative gas content
        '''
        
        A = C_1 * ((g* D**2 * rho_f)/sigma)**(1/8) \
            * ((g * D**3 * rho_f**2)/eta_f**2)**(1/12) \
            * (v_g/np.sqrt(g * D))
            
        f = lambda e_g: (e_g / (1-e_g)**4 - A)**2
        
        e_g_solved = minimize(f, 0.05).x
        
        return e_g_solved
        
    def calc_a(self, D, v_f, e_g, sigma, rho_f=None, g=9.81):
        '''Calculates the volume specific area in m^2/m^3
        
        D:          diameter of the bubble column in m
        v_f:        velocity of the fluid in m/s
        e_g:        relative gas content in the fluid
        sigma:      surface tension
        rho_f:      density of the fluid in kg/m^3
        
        returns:
        a:          volume speicifc area in m^2/m^3
        '''
        
        if rho_f is None:
            self.set_state_2()
        
        a = 1/(3*D) * ((g* D**2 * rho_f) / sigma)**(1/8) * ((g * D**3)/ v_f**2)**0.1 * e_g**1.13
        
        return a
        
     
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
        
    
    
    
    
def test():
    water = Liquid('Water')
    air = Gas('DryAir')
    bc = Bubble_Column(air, water)

    

if __name__ == "__main__":
    test()