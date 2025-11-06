import numpy as np
import math

import CoolProp.CoolProp as CP

#####################################################################################
# Channel geometry
#####################################################################################


###Values moved to global variables
#wall_thick = 0.0003
#inregen_wall_thick = 0.001

def channel_width_as_function_of_distance(Nc, r_values, wall_thick, inregen_wall_thick):
    Wc_values = np.zeros_like(r_values)
    Wc_values = 2*np.pi*(r_values + wall_thick)/Nc - inregen_wall_thick
    return Wc_values

def channel_area_as_function_of_distance(hc, Wc_values):
    Ac_values = np.zeros_like(Wc_values)
    Ac_values = hc*Wc_values
    return Ac_values

def calc_hydraulic_diameter(A, w):
    P = 2*w + 2*A/w
    Dh = 4*A/P
    return Dh


#####################################################################################
# Coolant analysis
#####################################################################################

def Re_no_as_function_of_distance(Wc_values, Ac_values, Tl_values, Pl_values, m_dot):
    #P = 50*101325
    #mu = 1036*10**(-6)

    Re = np.zeros_like(Ac_values)

    for i in range(len(Re)):
        # Conditions: Temperature (K), Pressure (Pa)
        mu = CP.PropsSI('VISCOSITY', 'T', Tl_values[i], 'P', Pl_values[i], 'Ethanol')
        print(mu)
        D = Wc_values[i]
        Re[i] = m_dot*D/(mu*Ac_values[i])


    #D = Wc_values 
    #Re = m_dot*D/(mu*Ac_values)
    return Re


def Solve_for_Friction_Factor(epsilon, D, Re, tol=1e-6, max_iter=100):
    """
    Solve the Colebrook-White equation iteratively for the Darcy friction factor, f.

    Parameters:
        epsilon (float): Roughness of the pipe (m).
        D (float): Diameter of the pipe (m).
        Re (float): Reynolds number (dimensionless).
        tol (float): Convergence tolerance (default: 1e-6).
        max_iter (int): Maximum number of iterations (default: 100).

    Returns:
        float: The Darcy friction factor, f.

    Raises:
        ValueError: If the iteration does not converge.
    """
    # Initial guess for f
    f = 0.02

    for i in range(max_iter):
        # Compute the RHS of the Colebrook-White equation
        rhs = -2 * math.log10((epsilon / (3.7 * D)) + (2.51 / (Re * math.sqrt(f))))
        
        # Update f using the LHS = RHS relationship
        f_new = 1 / (rhs ** 2)

        # Check for convergence
        if abs(f_new - f) < tol:
            return f_new

        f = f_new

    # If we reach this point, the iteration did not converge
    raise ValueError("Colebrook-White equation did not converge after {} iterations".format(max_iter))

# inputs: (distance along pipe, pipe area, hyd diamter, Liquid temp, Liquid pressure, mass flow)
def calc_Pressure_Loss(L, Ac, Dh, Tl, Pl, m_dot):
    epsilon = 5 * 10**(-5)  # standard for AM parts

    rho = CP.PropsSI('D', 'T', Tl, 'P', Pl, 'Ethanol')
    mu = CP.PropsSI('VISCOSITY', 'T', Tl, 'P', Pl, 'Ethanol')
    U = m_dot/(rho*Ac)

    Re = rho*Dh*U/mu

    f = Solve_for_Friction_Factor(epsilon, Dh, Re, tol=1e-6, max_iter=100)

    P_drop = f* 0.5*rho*U**(2)*L/Dh
    return P_drop

def get_rho_and_velocity_values(m_dot, Ac_values,Tl_values, Pl_values):
    rho_values = np.zeros_like(Tl_values)
    U_values = np.zeros_like(Tl_values)

    for i in range(len(rho_values)):
        rho_values[i] = CP.PropsSI('D', 'T', Tl_values[i], 'P', Pl_values[i], 'Ethanol')
        U_values[i] = m_dot/(rho_values[i]*Ac_values[i])

    return rho_values, U_values




