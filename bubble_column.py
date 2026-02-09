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
    
    D = None # Diameter of bubble column in m
    z = None # Height of the bubble column in m
    V_flow = None # Volume flow of the gas in m^3/s
    beta = 2.4e-2 # Mass transfer coefficient in m/s
                  # Taken from [Transportvorgänge in der Verfahrenstechnik - 2020 Springer Nature]
    
    
    def __init__(self, medium_1, medium_2, **kwargs):
        super().__init__(medium_1, medium_2, **kwargs)
        
    
    def calc_e_g(self, sigma, C_1=0.2, v_g=None, D=None, rho_f=None, eta_f=None, g=9.81):
        '''Calculates the relative gas content in the fluid
        For air and water typicall between 0.022 - 0.238  [Ind. Eng. Chem. Process Des. Develop., Vol. 12, No. 1, 1973]
        
        D:          diameter of the bubble column in m
        sigma:      surface tension of the fluid in N/m
        v_g:        velocity of the gas in m/s
        C_1:        initial gas concentration in kg/m^3
                    value defaults to 0.2 as seen in the source cited above
        rho_f:      density of the fluid in kg/m^3
        eta_f:      dynamic viscosity of the fluid in Pa*s
        g:          gravitational acceleration in m/s^2
        
        returns:
        e_g:        relative gas content
        '''
        
        if D is None and self.D is not None:
            D = self.D
        else:
            return Exception("D is not defined.")
        
        if self.check_none([rho_f, eta_f]):
            self.set_state_2()
            rho_f = self.medium_2.d
            eta_f = self.medium_2.eta
            
        if v_g is None:
            v_g = self.V_flow / (np.pi*D**2 / 4)
            
        if sigma is None and self.useCoolProp:
            sigma = self.medium_2.sigma
        
        A = C_1 * ((g* D**2 * rho_f)/sigma)**(1/8) \
            * ((g * D**3 * rho_f**2)/eta_f**2)**(1/12) \
            * (v_g/np.sqrt(g * D))
            
        f = lambda e_g: (e_g / (1-e_g)**4 - A)**2
        
        e_g_solved = minimize(f, 0.05).x[0]
        
        self.e_g = e_g_solved
        
        return e_g_solved
        
    def calc_a(self, e_g, sigma, nu_f=None, D=None, rho_f=None, g=9.81):
        '''Calculates the volume specific area in m^2/m^3
        
        D:          diameter of the bubble column in m
        nu_f:       kinematic viscosity of the fluid
        e_g:        relative gas content in the fluid
        sigma:      surface tension
        rho_f:      density of the fluid in kg/m^3
        
        returns:
        a:          volume speicifc area in m^2/m^3
        '''
        
        if D is None and self.D is not None:
            D = self.D
        elif D is None and self.D is None:
            return Exception("D is not defined.")
        
        if self.check_none([nu_f, rho_f]):
            self.set_state_2()
            nu_f = self.medium_2.eta / self.medium_2.d 
            rho_f = self.medium_2.d
            
        if sigma is None and self.useCoolProp:
            sigma = self.medium_2.sigma
        
        a = 1/(3*D) * ((g* D**2 * rho_f) / sigma)**(1/2) * ((g * D**3)/ nu_f**2)**0.1 * e_g**1.13
        
        self.a = a
        
        return a
        
     
    def calc_phi_z(self, phi_in, z, v_g=None, beta=None, a=None, **kwargs):
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
            
        if v_g is None:
            v_g = self.V_flow / (np.pi*self.D**2 / 4)
            
        phi_z = 1 - (1 - phi_in)*np.exp(- (beta*a)/v_g * z)
        
        return phi_z
        
    
#%% Test functions    
    
water = Liquid('Water')
air = Gas('DryAir')
bc = Bubble_Column(air, water)
bc_cp = Bubble_Column("Air", "Water", useCoolProp=True)
bc.T_1 = [298.15, 298.15]
bc.T_2 = [298.15, 298.15]
bc.p_1 = 1e5
bc.p_2 = 1e5  

bc_cp.T_1 = [298.15, 298.15]
bc_cp.T_2 = [298.15, 298.15]
bc_cp.p_1 = 1e5
bc_cp.p_2 = 1e5  
 
def test():
    
    isSuccessful = []
    
    bc.V_flow = 0.05 / 3600 # m^3/s
    bc.D = 0.01
    bc_cp.V_flow = 0.05 / 3600 # m^3/s
    bc_cp.D = 0.01
    
    try:
        e_g = bc.calc_e_g(sigma=0.072) # Sigma taken from [Transportvorgänge in der Verfahrenstechnik - 2020 Springer Nature]
        e_g_cp = bc_cp.calc_e_g(sigma=None)
        print("Calculation of e_g successful")
        ic(e_g)
        ic(e_g_cp)
        isSuccessful.append(True)
    except Exception as e:
        print("Calculation of e_g failed", e)
        e_g = 0.03
        isSuccessful.append(False)
    
    try:
        a = bc.calc_a(e_g, sigma=0.072)
        a_cp = bc_cp.calc_a(e_g, sigma=None)
        print("Calculation of a successful")
        ic(a)
        ic(a_cp)
        isSuccessful.append(True)
    except Exception as e:
        print("Calculation of a failed", e)
        a = 30
        isSuccessful.append(False)
    
    try:
        phi_out = bc.calc_phi_z(0.5, 0.13)
        phi_out_cp = bc_cp.calc_phi_z(0.5, 0.13)
        print("Calculation of phi_z successful")
        ic(phi_out)
        ic(phi_out_cp)
        isSuccessful.append(True)
    except Exception as e:
        print("Calculation of phi_z failed", e)
        isSuccessful.append(False)
    
    return all(isSuccessful)
    
    
def phi_out_standard(V_flow, D, z, phi_in=0, sigma=0.072):
    
    bc.V_flow = V_flow
    bc.D = D
    
    e_g = bc.calc_e_g(sigma=sigma)
    bc.calc_a(e_g, sigma=sigma)
    phi_out = bc.calc_phi_z(phi_in, z)
    
    return phi_out

if __name__ == "__main__":
    print(test())
    
    ic(phi_out_standard(V_flow=100 / 60000, D=0.1, z=0.08))