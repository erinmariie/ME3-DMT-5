# Functions: Heat Transfer Functions

import numpy as np

from nozzle_functions import conical_area_function,parabolic_area_function, mach_no_as_function_of_distance, temperature_as_function_of_distance
from regen_channel_functions import calc_hydraulic_diameter, calc_Pressure_Loss

import time
import math
from scipy.interpolate import interp1d

import CoolProp.CoolProp as CP

#####################################################################################
# Bertz' equation (hg, gass convective ht tran coef as function of distance)
#####################################################################################
"""
def Prandlt_estimate(gamma): # OG
    return (4*gamma/(9*gamma-5))  # check this !
"""
###Variables originally in this file for 2024
#emissivity = 0.06 # dimensionless,  approx average value obtained from calculator (spreadsheet): https://www.sciencedirect.com/science/article/pii/S0022407318304618
#rc = 0.003 #throat radius of curvature ###MOVE TO VARIABLES
#mu_0 = 6.848*10**(-5)
#Cp_0 = 1.9
#Pr_0 = 0.5758 #mu_0*Cp_0/k_0

#NEED TO UPDATE THE PRANDTL ESTIMATE
def Prandtl_estimate(gamma): # 2
    Pr = 0.5837 # from NASA CEA
    return Pr


def hg_as_function_of_distance(P0,A0, At, Ae, Lc, Lt, L, b, R2, R1, Rn, Tn, Te, gamma,R,Cp,T0,Tw,rc):
    Pr = Prandtl_estimate(gamma)
    #mu = 7.02257 *10**(-5)  ## Temporary## ## need to replace with function for calculation ############

    # c*1
    #c_star = np.sqrt(R*T0)/(np.sqrt(2*gamma/(gamma-1))*(2/(gamma+1))**((gamma+1)/(2*(gamma-1)))) ## may require update ##### based on research ############
    
    # c*2
    c_star = np.sqrt(gamma*R*T0)/(gamma*np.sqrt((2/(gamma+1))**((gamma+1)/(gamma-1)))) ##check !
    
    # c*3
    #c_star = np.sqrt(gamma*R*T0)/(gamma*(2/(gamma+1))**((gamma+1)/(2*gamma-2)))
    
    

    Dt = np.sqrt((4*At/np.pi))

    (x_values, Area_values, throat_index, chamber_end_index) = parabolic_area_function(A0, At, Ae, Lc, Lt, L, b, R2, R1, Rn, Tn, Te)
    (x_values, M_values) = mach_no_as_function_of_distance(A0, At, Ae, Lc, Lt, L, b, R2, R1, Rn, Tn, Te, gamma)
    (x_values, T_values) = temperature_as_function_of_distance(A0, At, Ae, Lc, Lt, L, b, R2, R1, Rn, Tn, Te, gamma,T0)

    mu_values = np.zeros_like(x_values)
    for i in range(len(mu_values)):
        mu_values[i] = calculate_viscosity((0,chamber_end_index,throat_index,-1), i, T_values[i])

    sigma = 1/((0.5*(Tw/T0)*(1+(gamma-1)*(M_values**2)/2)+0.5)**(0.68)*(1+(gamma-1)*(M_values**2)/2)**(0.12))  #correction factor
    hg_values = ((0.026/Dt**0.2)*((mu_values**0.2)*Cp/(Pr**0.6))*(P0/c_star)**(0.8)*(Dt/rc)**(0.1))*(At/Area_values)**(0.9)*sigma
    return (x_values,hg_values)



def calculate_hg(T, i, P0, At, gamma,R,Cp,T0,Tw, Area, M, chamber_index, throat_index, m_dot_total,rc):
    Pr = Prandtl_estimate(gamma)
    #Pr = 0.61
    #mu = 7.02257 *10**(-5)  ## Temporary## ## need to replace with function for calculation ############
    #mu = 6.5*10**(-5)
    mu = calculate_viscosity((0,chamber_index,throat_index,-1), i, T)

    # c*1
    #c_star = np.sqrt(R*T0)/(np.sqrt(2*gamma/(gamma-1))*(2/(gamma+1))**((gamma+1)/(2*(gamma-1)))) ## may require update ##### based on research ############
    
    # c*2
    #c_star = np.sqrt(gamma*R*T0)/(gamma*np.sqrt((2/(gamma+1))**((gamma+1)/(gamma-1))))
    
    # c*3 # effectivly the same as c*2
    #c_star = np.sqrt(gamma*R*T0)/(gamma*(2/(gamma+1))**((gamma+1)/(2*gamma-2)))

    # c*4 # little change
    #c_star = np.sqrt(gamma*R*T)/(gamma*(2/(gamma+1))**((gamma+1)/(2*gamma-2)))

    # c*5 # litle change
    #c_star = P0*At/m_dot

    # c*6
    c_star = np.sqrt((1/gamma)*((gamma+1)/2)**((gamma+1)/(gamma-1))*(R*T0/20)) # molecular mass of ethanol= 46, weigthed avg, M = 20
    
    
    r = np.sqrt(Area/np.pi)
    Dt = np.sqrt((4*At/np.pi))

    sigma = 1/((0.5*(Tw/T0)*(1+(gamma-1)*(M**2)/2)+0.5)**(0.68)*(1+(gamma-1)*(M**2)/2)**(0.12))  #correction factor
    hg = ((0.0026/(Dt**0.2))*((mu**0.2)*Cp/(Pr**0.6))*(P0/c_star)**(0.8)*(Dt/rc)**(0.1))*(At/Area)**(0.9)*sigma
    return (hg)

def calculate_hg2(T, P, rho, i, At,Tw, M, vel, chamber_index, throat_index,mu_0,Cp_0,Pr_0):
    ### may or may not be Dt (could be local)
    #T0 = T(1+((gamma-1)*M**2)/2)
    Tm = (T+Tw)/2

    

    mu_m = calculate_viscosity((0,chamber_index,throat_index,-1), i, Tm)
    rho_m = calculate_density((0,chamber_index,throat_index,-1), i, Tm, P)
      
    Dt = np.sqrt((4*At/np.pi))
    hg = 0.026/(Dt**0.2)*(mu_0**0.2*Cp_0/(Pr_0**(0.6)))*(rho*vel)**(0.8)*(rho_m/rho)**(0.8)*(mu_m/mu_0)**(0.2)
    return (hg)

#####################################################################################
# Calc Comb, gas viscosity + desity from mass fraction values
#####################################################################################

def calculate_viscosity(indices, target_index, target_temperature):
    # Constants for viscosity calculation
    con1000 = [
        [0.65060585, 28.517449, -16690.236, 1.5223271],  # N2
        [0.50714993, -689.66913, 87454.75, 3.0285155],  # H2O
        [0.65060585, 28.517449, -16690.236, 1.5223271],  # CO
        [0.65318879, 51.738759, -62834.882, 1.5227045]   # CO2
    ]
    
    con300 = [
        [0.60443938, -43.632704, -884.41949, 1.897215],  # N2
        [0.7838778, -382.60408, 49040.158, 0.85222785],  # H2O
        [0.60443938, -43.632704, -884.41949, 1.897215],  # CO
        [-0.54330318, -188.23898, 8872.6567, 2.4499362]  # CO2
    ]
    
    # Mass fractions at given indices
    mass_fraction = [
        [0.4628923, 0.4628923, 0.4629015, 0.4629033],  # N2
        [0.1590833, 0.1590833, 0.1557418, 0.1231363],  # H2O
        [0.2821263, 0.2821263, 0.2767468, 0.2259908],  # CO
        [0.0777402, 0.0777402, 0.0861982, 0.1659354]   # CO2
    ]
    
    # Interpolate mass fractions for the given target index
    mass_frac_interp = np.zeros(4)
    for i in range(4):
        mass_frac_interp[i] = np.interp(target_index, indices, mass_fraction[i])
    
    # Calculate viscosity
    viscosity = 0
    mass_sum = 0
    for j in range(4):
        if target_temperature < 1000:
            viscosity += mass_frac_interp[j] * math.exp(
                con300[j][0] * math.log(target_temperature) +
                con300[j][1] / target_temperature +
                con300[j][2] / target_temperature ** 2 +
                con300[j][3]
            )
        else:
            viscosity += mass_frac_interp[j] * math.exp(
                con1000[j][0] * math.log(target_temperature) +
                con1000[j][1] / target_temperature +
                con1000[j][2] / target_temperature ** 2 +
                con1000[j][3]
            )
        mass_sum += mass_frac_interp[j]
    
    viscosity_micropoise = viscosity / mass_sum if mass_sum != 0 else 0  
    viscosity_Pa_s = viscosity_micropoise*10**(-7) # 1microP = 10^-7 Pas

    return viscosity_Pa_s # in Pascal seconds


def calculate_density(indices, target_index, target_temperature, pressure_pa):
    # Molar masses in kg/mol: N2, H2O, CO, CO2
    molar_masses = [0.0280134, 0.0180153, 0.0280101, 0.0440095]

    # Mass fractions as in your code...
    mass_fraction = [
        [0.4628923, 0.4628923, 0.4629015, 0.4629033],  # N2
        [0.1590833, 0.1590833, 0.1557418, 0.1231363],  # H2O
        [0.2821263, 0.2821263, 0.2767468, 0.2259908],  # CO
        [0.0777402, 0.0777402, 0.0861982, 0.1659354]   # CO2
    ]

    # Interpolate mass fractions
    mass_frac_interp = np.zeros(4)
    for i in range(4):
        mass_frac_interp[i] = np.interp(target_index, indices, mass_fraction[i])

    # Convert to mole fractions
    mole_amounts = [mass_frac_interp[i] / molar_masses[i] for i in range(4)]
    total_moles = sum(mole_amounts)
    mole_fractions = [n / total_moles for n in mole_amounts]

    # Calculate average molar mass
    avg_molar_mass = sum([mole_fractions[i] * molar_masses[i] for i in range(4)])

    # Ideal gas law: ρ = P * M / (R * T)
    R = 8.3145  # J/(mol·K)
    density = pressure_pa * avg_molar_mass / (R * target_temperature)  # in kg/m³

    return density




#####################################################################################
# Cooling convection (hl, convective ht tran coef calc)
#####################################################################################

def constant_hl_value(Cp, mu, k, m_dot_l, area, D):
    hl = 0.023*Cp*(m_dot_l/area)*(D*m_dot_l/(area*mu))**(-0.2)*(mu*Cp/k)**(-2/3)
    return hl

# inputs: (single channel m_dot, A, Dh, Cp, mu, Kl)
def calc_hl_value(m_dot_l, A, Dh, Cp, mu, Kl):
    hl = 0.023*(Cp*m_dot_l/A)*(Dh*m_dot_l/A*mu)**(-0.2)*(mu*Cp/Kl)**(-2/3)
    return hl

# improved
# inputs: (single channel m_dot_l1, A, Dh, Cp, mu, Kl)
def calc_hl_value2(m_dot_l1, A, Dh, Cp, mu, Kl, Twl, Tl):
    Re = m_dot_l1*Dh/(A*mu)
    Pr = mu*Cp/Kl
    hl = 0.023*(Kl/Dh)*Re**(0.8)*Pr**(0.4)*(Twl/Tl)**(-0.3)
    return hl

#####################################################################################
# Sequential Heat Transfer Analysis steps
#####################################################################################

# inputs: (Wall temp guess, Gas temp, Gas Mach, Gas gamma, hg from Bartz)
def calc_heat_flux_gas_side(Twg, Tg, M, gamma, hg ):
    Pr = Prandtl_estimate(gamma)
    r = Pr**(1/3)
    Tr = Tg*(1+((gamma-1)*r*M**2)/2)

    q_dot = hg *(Tr-Twg)
    return q_dot

# inputs: (Wall temp guess, Gas temp, Gas Mach, Gas gamma, hg from Bartz)
def calc_radiation_heat_flux_gas_side(Tg,emissivity):
    Stefan_Boltzman_Const = 5.67*10**(-8)  # W/m**2K**4
    q_dot = emissivity*Stefan_Boltzman_Const*(Tg)**(4)
    return q_dot

# inputs: (Wall temp guess, heat flux calculated gas side, Wall conductivity, Wall thickness)
def calc_Twl_from_Twg(Twg, q_dot_g, K, tw):
    Twl = Twg - q_dot_g*tw/K
    return Twl

# inputs: (Wall temp fluid side, coolant temp in previous step, Conv. heat tran. coef. of coolant)
def calc_heat_flux_coolant_side(Twl, Tl_prev, hl):
    q_dot = hl*(Twl-Tl_prev)
    return q_dot

# inputs: (correct (converged) heat flux, coolant temp in previous step, coolant mass flow, coolant Cp, chanel width at point, x step at point)
def calc_new_coolant_temp(q_dot, Tl_prev, ml_dot, Cpl, Wc, x_step):
    Integral = q_dot*Wc*x_step
    Tl = Tl_prev + Integral/(ml_dot*Cpl)
    return Tl

# inputs: (heat flux calc coolant side, gas temp, Mach, gamma, hg)
def calc_Twg_from_heat_flux(q_dot_c, Tg, M, gamma, hg): #used in finding a more accurate Twg guess
    Pr = Prandtl_estimate(gamma)
    r = Pr**(1/3)  #Recovery
    Tr = Tg*(1+((gamma-1)*r*M**2)/2)

    Twg_c = Tr - q_dot_c/hg
    return Twg_c


#####################################################################################
# Sequential Heat Transfer Analysis (basic)
#####################################################################################


#####################################################################################
# Sequential Heat Transfer Analysis extended to radiation + presure drop
#####################################################################################

# inputs: (x values, r values, Mach values, gas temp values, hg values, gas gamma, channel width values, channel area values, inner wall thickness, wall conductivity, input fuel temp.)
def Sequential_Heat_Transfer_Analysis_Extended3(x_values, r_values, M_values,Vg_values, Tg_values, Pg_values,rhog_values, hg_values, gamma, Wc_values, Ac_values, tw, K, Tl_in, Pl_in, m_dot_l, A_values, chamber_index, throat_index, R, P0, Cpg, m_dot):
    max_error = 0.01 #assign accuracy required
    
    Twg_values = np.zeros_like(r_values)
    Twl_values = np.zeros_like(r_values)
    Tl_values = np.zeros_like(r_values)
    
    q_dot_g_conv_values = np.zeros_like(r_values)
    q_dot_g_rad_values = np.zeros_like(r_values)
    q_dot_values = np.zeros_like(r_values)

    Pl_values = np.zeros_like(r_values)
  

    for i in reversed(range(len(Twg_values))):
        q_dot_g_rad = calc_radiation_heat_flux_gas_side(Tg_values[i])
        q_dot_found = False
        #print(i)
        Twg_test_values = np.linspace(Tl_in, np.max(Tg_values), 10) #test 10 values from minimum to maximum possible temps
        #print("Test Twg_test_values:", Twg_test_values)
        #best_guess = (None, None, None, None, None) # (current best guess Twg, q dot error, Twg index, Twl, q_dot_g)
        #error_prev = None
        while (not q_dot_found):
            best_guess = (None, None, None, None, None) # reset best guess so answere doesnt stall
            for I in range(len(Twg_test_values)):

                #find hg using Twg_test_value
                hg = calculate_hg2(Tg_values[i], Pg_values[i], rhog_values[i], i, A_values[throat_index],Twg_test_values[I], M_values[i], Vg_values[i], chamber_index, throat_index)
                print("######### hg ##########", hg)

                print("Twg guess:", Twg_values[i])
                q_dot_g_conv = calc_heat_flux_gas_side(Twg_test_values[I], Tg_values[i], M_values[i], gamma, hg)
                ############
                # function to find q_dot_g_rad
                ############
                q_dot_g = q_dot_g_conv + q_dot_g_rad

                #print("q_dot_g:", q_dot_g)
                Twl = calc_Twl_from_Twg(Twg_test_values[I], q_dot_g, K, tw)
                ######
                if Twl <= 0 or np.isnan(Twl) or np.isinf(Twl):
                    continue
                ######
                #print("Twl guess:", Twl)
                if (i==(len(Twg_values)-1)):
                    Tl_prev = Tl_in
                else:
                    Tl_prev = Tl_values[i+1]

                if (i==(len(Twg_values)-1)):
                    Pl_prev = Pl_in
                else:
                    Pl_prev = Pl_values[i+1]

                Cp = CP.PropsSI('Cpmass', 'T', Tl_prev, 'P', Pl_prev, 'Ethanol') #Specific Heat at Constant Pressure (Cp) in J/(kg·K)
                mu = CP.PropsSI('VISCOSITY', 'T', Tl_prev, 'P', Pl_prev, 'Ethanol') #Dynamic Viscosity (μ) in Pa·s
                Kl = CP.PropsSI('CONDUCTIVITY', 'T', Tl_prev, 'P', Pl_prev, 'Ethanol') #Thermal Conductivity (k) in W/(m·K)
                Dh = calc_hydraulic_diameter(Ac_values[i], Wc_values[i])
                
                print("##############",m_dot_l,"#",Ac_values[i],"#",Dh,"#",Cp,"#",mu,"#",Kl,"#",Twl,"#",Tl_prev)
                #hl = calc_hl_value(m_dot_l, Ac_values[i], Dh, Cp, mu, Kl)
                hl = calc_hl_value2(m_dot_l, Ac_values[i], Dh, Cp, mu, Kl, Twl, Tl_prev)
                #print("######### hl ##########", hl)
                #print("hl:", hl)
                q_dot_c = calc_heat_flux_coolant_side(Twl, Tl_prev, hl)
                #print("q_dot_c:", q_dot_c)

                error = abs(q_dot_g - q_dot_c) / max(abs(q_dot_g), abs(q_dot_c))

                if (best_guess[0] == None):
                    best_guess = (Twg_test_values[I], error, I, Twl, q_dot_g)
                if (error < best_guess[1]):
                    best_guess = (Twg_test_values[I], error, I, Twl, q_dot_g)


            #print(Twg_test_values)
            #print("******best guess:",best_guess)
            
            #error_convergence_threshold = 0.001  # Threshold for convergence

            if (best_guess[1] <= max_error):
                q_dot_found = True


            else:
                # the best test value is identified
                
                #Twg_test_values = np.linspace(Twg_test_values[best_guess[2]-1], Twg_test_values[best_guess[2]+1], 10)
                if best_guess[2] == 0:
                    # If best_guess[2] is the first index, use the first two elements for the range
                    Twg_test_values = np.linspace(Twg_test_values[best_guess[2]], Twg_test_values[best_guess[2] + 1], 10)
                elif best_guess[2] == len(Twg_test_values) - 1:
                    # If best_guess[2] is the last index, use the last two elements for the range
                    Twg_test_values = np.linspace(Twg_test_values[best_guess[2] - 1], Twg_test_values[best_guess[2]], 10)
                else:
                    # Normal case: use the previous and next indices
                    Twg_test_values = np.linspace(Twg_test_values[best_guess[2] - 1], Twg_test_values[best_guess[2] + 1], 10)

                #error_prev = best_guess[1]
                #print(Twg_test_values)
                #time.sleep(1)  # Delay seconds
        
        Twg_values[i] = best_guess[0]
        Twl_values[i] = best_guess[3]

        q_dot_g_rad_values[i] = q_dot_g_rad
        q_dot_g_conv_values[i] = best_guess[4] - q_dot_g_rad
        q_dot_values[i] = best_guess[4]


        if (i==(len(Twg_values)-1)):
            x_step = x_values[i] - x_values[i-1]
        else:
            x_step = x_values[i+1] - x_values[i]

        Tl_values[i] = calc_new_coolant_temp(q_dot_values[i], Tl_prev, m_dot_l, Cp, Wc_values[i], x_step)

        P_drop = calc_Pressure_Loss(x_step, Ac_values[i], Dh, Tl_values[i], Pl_prev, m_dot_l)
        if (i==(len(Twg_values)-1)):
                    Pl_values[i] = Pl_in
        else:
            Pl_values[i] = Pl_values[i+1] - P_drop
        print("hl:", hl)
        print("Twg:",Twg_values[i])
        print("Twl:",Twl_values[i])
        print("q_dot_conv:",q_dot_g_conv_values[i])
        print("q_dot_rad:",q_dot_g_rad_values[i])
        print("q_dot:",q_dot_values[i])
        print("Tl:",Tl_values[i])
        print("Pl:",Pl_values[i])
        print("---------------------------------------------------------")

    return (Twg_values, Twl_values, Tl_values,Pl_values, q_dot_values, q_dot_g_conv_values, q_dot_g_rad_values)



#####################################################################################
# Account for film cooling
#####################################################################################

def adjust_for_film_cooling(x_values, values, x_vals_no, vals_no, x_vals_with, vals_with, kind='linear'):
    """
    Interpolates vals_no and vals_with to x_values and subtracts their difference from values.

    Parameters:
    - x_values: high-resolution x-grid corresponding to `values`
    - values: high-resolution list of values (e.g., Twg_values)
    - x_vals_no: x-points for vals_no (e.g., x_vals_no)
    - vals_no: lower-res list without film cooling (e.g., Twg_vals_no)
    - x_vals_with: x-points for vals_with (e.g., x_vals_with)
    - vals_with: lower-res list with film cooling (e.g., Twg_vals_with)
    - kind: interpolation type (e.g., 'linear', 'cubic')

    Returns:
    - values_film: np.ndarray of corrected values
    """
    
    x_values = np.array(x_values)
    values = np.array(values)

    interp_no = interp1d(x_vals_no, vals_no, kind=kind, fill_value="extrapolate")
    interp_with = interp1d(x_vals_with, vals_with, kind=kind, fill_value="extrapolate")

    vals_no_interp = interp_no(x_values)
    vals_with_interp = interp_with(x_values)

    values_film = values - (vals_no_interp - vals_with_interp)
    return values_film

