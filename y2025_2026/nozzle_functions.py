# Functions: nozzle fluid dynamics calculations

import numpy as np

#####################################################################################
# fluids equations
#####################################################################################

def exit_mach_no(gamma,P0,Pe):
    Me = np.sqrt((2/(gamma-1))*((P0/Pe)**((gamma-1)/gamma)-1))
    return Me

def exit_temp(T0, gamma, Me):
    Te = T0/(1+(gamma-1)*(Me**2)/2)
    return Te

def exit_velocity(Me,gamma, R, Te):
    Ve = Me*np.sqrt(gamma*R*Te)
    return Ve

def mass_flow_rate(F, Ve):
    m_dot = F/Ve
    return m_dot

def throat_area(m_dot, P0, T0, gamma,R ):
    At = m_dot/(P0*(np.sqrt(gamma/(R*T0)))*((2/(gamma+1))**((gamma+1)/(2*(gamma-1)))))
    return At

def exit_area(Me, At, gamma):
    Ae = At*(2*(1+((gamma-1)/2)*Me**2)/(gamma+1))**((gamma+1)/(2*(gamma-1)))/Me
    return Ae

def chamber_sizing(char_L, At, CR):
    Vc = At*char_L
    A0 = CR*At
    Lc = Vc/A0
    return (A0, Lc)

def area_to_radius(A):
    r = np.sqrt(A/np.pi)
    return r

#####################################################################################
# Nozzle geometry
#####################################################################################

def chamber_nozzle_sizing(Pe,P0,F,CR,Lstar, T0,gamma,R):
    Me = exit_mach_no(gamma,P0,Pe)
    Te = exit_temp(T0, gamma, Me)
    Ve = exit_velocity(Me,gamma, R, Te)
    m_dot = mass_flow_rate(F, Ve)
    At = throat_area(m_dot, P0, T0, gamma,R )
    Ae = exit_area(Me, At, gamma)
    (A0, Lc) = chamber_sizing(Lstar, At, CR)
    return (A0,At,Ae,Lc,m_dot)
