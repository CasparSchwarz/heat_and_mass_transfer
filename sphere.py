# -*- coding: utf-8 -*-
"""
Created on Thu Aug 21 15:57:47 2025

@author: schwarz
"""

from TILMedia import Gas, Liquid
from statistics import mean
from heat_exchanger import Heat_exchanger
import numpy as np
from icecream import ic

class Sphere(Heat_exchanger):
    
    def __init__(self, medium_1):
        super().__init__(medium_1, 'Dummy')
        
        
    # Calculates Nu for sphere according to Gnielinski (VDI Wärmeatlas)
    def calc_Nu(self, Nu_lam=None, Nu_turb=None, Re=None, Pr=None):
        
        if Nu_lam is None:
            Nu_lam = self.calc_Nu_lam(Re, Pr)
        if Nu_turb is None:
            Nu_turb = self.calc_Nu_turb(Re, Pr)
            
        
        return 2+ np.sqrt(Nu_lam**2+Nu_turb**2)
    
    def calc_Nu_lam(self, Re, Pr=None):
        if Pr is None:
            self.set_state_1()
            Pr = self.medium_1.Pr
        
        return 0.664*Re**(1./2.)*Pr**(1./3.)
    
    def calc_Nu_turb(self, Re, Pr=None):
        if Pr is None:
            self.set_state_1()
            Pr = self.medium_1.Pr
        
        return (0.037*Re**(0.8)*Pr) / (1+2.443*Re**(-0.1)*(Pr**(2./3.)-1))
    
    def calc_Re(self, d_k, v=None, rho_v=None, rho=None, eta=None):
        ''' Calculation of Reynolds number for sphreres
        d_k:    Particle diameter in m
        v:      Fluid velocity in m/s
        rho_v:  Density multiplied with velocity
        rho:    Density in kg/m^3
        eta:    Dynamic viscosity in Pa*s
        '''
        
        if eta is None:
            self.set_state_1()
            eta = self.medium_1.eta # Pas
        if rho is None:
            self.set_state_1()
            rho = self.medium_1.d # kg/m³
        if rho_v is not None:
            rho = 1
            v = rho_v
            
        self.Re = rho*v*d_k/eta
        
        return self.Re
        
    def calc_T_1_x(self, T_w, alpha_1=None, m_flow_1=None, c_p=None, A=None):
        '''Calculate the temperature of medium_1 at position x
        x is predetermined by the area A
        
        T_w:     Wall temperature in K
        alpha_1: Heat transfer coefficient between medium_1 and wall in W/m^2K
        m_flow_1: Mass flow of medium 1 in kg/s
        c_p:      Specific heat capacity of medium_1 in kJ/kgK
        A:        Heat exchanger Area in m^2
        '''
        
        if alpha_1 is None: alpha_1 = self.alpha_1
        if m_flow_1 is None: m_flow_1 = self.m_flow_1
        if c_p is None:
            self.medium_1.setState_pTxi(self.p_1, mean(self.T_1))
            c_p = self.medium_1.cp
        if A is None:
            A = self.A
            
        T_1_x = T_w-(T_w-self.T_1[0])*np.exp(-alpha_1 * A/(m_flow_1*c_p))
        
        return T_1_x
        
        
class Particle_Bed(Sphere):
    
    # psi = V_empty / V_total = 1 - rho_total / rho_sphere
    # Lewatit: rho_total = 630 g/l
    # closest packing: 74 % -> Lewatit rho_sphere = 851 g/l
    psi = 0.26 # psi for closest packing of spheres
    f_a = None
    
    
    def __init__(self, medium):
        super().__init__(medium)
        
    def calc_Nu(self, f_a=None, psi=None, **kwargs): 
        
    
        if psi is None:
            psi=self.psi 
            
        if f_a is None:
            f_a = self.calc_f_a(psi)
        
        return super().calc_Nu(**kwargs) * f_a
    
    
    # Calc scaling factor f_a from porosity psi
    # psi = V_empty / V_total
    def calc_f_a(self, psi):
        
        self.f_a = 1 + 1.5*(1-psi)
        
        return self.f_a
    
    
    def calc_A(self, d_k, A_bed_total, L_bed):
        
        V_k = 4/3*(d_k/2)**3*np.pi # Volume of particle
        A_k = 4*(d_k/2)**2*np.pi # Surface area of particle
        
        V_bed_total = A_bed_total*L_bed # Volume of particle bed
        n_k = V_bed_total / V_k # Calculate number of particles
        
        self.A = n_k*A_k # Heat exchange area
        
        return self.A # m^2
    
    # Calculate the temperature of medium_1 at the end of the heat excahgner
    def calc_T_end(self, L, psi, T_w=383.15, d_k=100e-6, dBed=6e-3, m_flow=1.9470702094005574e-06): # flow at 100 ml/min air
        
        A_bed_total = np.pi*dBed**2 /4
        rho_v = m_flow/A_bed_total
        Re = self.calc_Re(d_k, rho_v=rho_v)
        Nu = self.calc_Nu(psi=psi, Re=Re)
        alpha = self.calc_alpha_1(Nu, d_k)
        ic(alpha)
        A = self.calc_A(d_k, A_bed_total, L)
        T_end = self.calc_T_1_x(T_w, alpha, m_flow, A=A)
        return T_end
    
    # Pressure drop calculations
    
    # Euler number
    def calc_Eu_generic(self, dp, d_p, psi, dL, v, rho=None):
        '''
        Euler number describes the relation of a pressure drop to the kinetic
        energy of the fluid experiencing the pressure drop.

        Parameters
        ----------
        dp : float
            pressure drop in Pa
        d_p : floa
            particle diameter in m
        psi : float
            porosity
        dL : float
            length across which the pressure drop occurs in m
        v : float
            velocity of the fluid in m/2
        rho : float, optional
            density in kg/m^3, will be calculated if None is given

        Returns
        -------
        flaot
            Euler number.

        '''
        
        if rho is None:
            self.set_state_1()
            rho = self.medium_1.d
            
        self.Eu = 4/3 * dp/(rho*v**2) * d_p/dL * psi**2/(1-psi)
        
        return self.Eu
    
    def calc_Re(self, *args, **kwargs):
        Re = super().calc_Re(*args, **kwargs) / self.psi # Correct Reynolds number for particle bed
        
        return Re
        
    
    def calc_Eu(self, psi, Re=None, **kwargs):
        '''Correlation for Euler number from VDI Wärmeatlas
        '''
        
        if Re is None:
            Re = self.calc_Re(kwargs)
        
        # Calculate constant based on porosity
        r_delta = 1 / (0.95 / (1-psi)**(1./3.) - 1)
        
        # Calculate parts of final equation
        A_var = 24/Re * (1 + 0.692 * (r_delta + 0.5* r_delta**2) )
        B_var = 4 / np.sqrt(Re) * (1 + 0.12 * r_delta**(1.5) )
        C_var = (0.4 + 0.891 * r_delta * Re**(-0.1) )
        
        self.Eu = A_var + B_var + C_var
        
        return self.Eu
    
    
    def calc_dp(self, d_p, dL, psi, v, Eu=None, Re=None, rho=None):
        '''
        Calculates the pressure drop with Euler number
        '''
        
        if rho is None:
            self.set_state_1()
            rho = self.medium_1.d
            
        if psi is None:
            psi = self.psi
            
        if Eu is None:
            Eu = self.calc_Eu(psi, Re)
            
        dp = 3/4 * Eu * (rho*v**2) * dL/d_p * (1-psi)/psi**2
        
        return dp
    
    
    # Top level method to calculate the pressure drop across the heat exchanger
    def calc_pressure_drop(self, L, d_k=100e-6, dBed=6e-3, m_flow=1.9470702094005574e-06, psi=None):
        
        A_bed_total = np.pi*dBed**2 / 4
        rho_v = m_flow/A_bed_total
        self.set_state_1()
        v = rho_v / self.medium_1.d
        Re = self.calc_Re(d_k, rho_v = rho_v)
        ic(Re)
        
        dp = self.calc_dp(d_k, L, psi, v, Re=Re)
        
        return dp
    

    
#%% ##### TESTER FUNCTIONS #######################################################
pb = Particle_Bed(Gas("DryAir"))
pb.T_1 = [298.15,298.15]
pb.T_2 = [383.15,383.15]
pb.p_1 = 100000   
def temperature_and_pressureDrop_test(): 
    '''Test function for temperature and pressureDrop across packed bed

    '''
    
    passes = []

    try:
        ic(pb.calc_T_end(0.01, psi=0.5, T_w=383.15, d_k=100e-6, m_flow=1.9470702094005574e-06)) # flow at 100 ml/min air
        passes.append(True)
    except Exception as e:
        ic(f'Calculation of T_end unsuccessful: {e}')
        passes.append(False)
        
    try:
        ic(pb.calc_pressure_drop(0.004, psi=None, d_k=70e-6, dBed=12.7e-3, m_flow=1.9470702094005574e-06)) # flow at 100 ml/min air
        passes.append(True)
    except Exception as e:
        ic(f'Calculation of pressure drop unsuccessful: {e}')
        passes.append(False)
        
    if all(passes): 
        ic('test passed')
    else:
        ic('test unsuccessful')
    ic(passes)
        

if __name__ == '__main__':
    
    temperature_and_pressureDrop_test()
    