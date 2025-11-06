### DETAILED ANALYSIS
#prepares analysis for main analysis
from nozzle_functions import chamber_nozzle_sizing 


def get_nozzle_functions():
    from variables import Pe, P0, F, CR, Lstar, T0, gamma, R
    (A0, At, Ae, Lc, m_dot) = chamber_nozzle_sizing(Pe, P0, F, CR, Lstar, T0, gamma, R)
    print((A0, At, Ae, Lc, m_dot))
    


