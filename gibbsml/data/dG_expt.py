"""
Experimental free-energies for different metal oxides.  Taken from
the NIST-JANAF [1] and Cambridge Dissemination of IT for the Promotion
of Materials Science (DoITPoMS) [2] databases.

[1] Malcolm W. Chase, Jr. NIST-JANAF Thermochemical Tables.
    Washington, DC : New York :
    American Chemical Society; American Institute of Physics
    for the National Institute of Standards and Technology, 1998.

[2] https://www.doitpoms.ac.uk/tlplib/ellingham_diagrams/index.php

"""

__author__ = "Jose A. Garrido Torres, Alexander Urban"
__email__ = "aurban@atomistic.net"
__date__ = "2021-03-18"
__version__ = "1.0"


def get_dG0_expt(T):

    # MO2:
    dG0_TiO2_rut = -941 + 0.18 * T  # Rutile.
    dG0_SiO2_qu  = -907 + 0.18 * T  # Quartz.
    dG0_SnO2     = -575 + 0.2  * T
    dG0_MnO2     = -519 + 0.18 * T
    dG0_ZrO2     = -1092 + 0.18 * T
    dG0_MoO2     = -578 + 0.17 * T
    dG0_WO2      = -581 + 0.17 * T
    dG0_KO2      = -282 + 0.15 * T
    dG0_NbO2     = -784 + 0.17 * T

    # MO:
    dG0_MgO   = -1202 + 0.22 * T
    dG0_CaO   = -1280 + 0.22 * T
    dG0_VO    = -849 + 0.16 * T
    dG0_MnO   = -778 + 0.15  * T
    dG0_FeO   = -538 + 0.14 * T   # Warning Fe --> Fe0.947O.
    dG0_CoO   = -491 + 0.16 * T
    dG0_NiO   = -471 + 0.17  * T
    dG0_CuO   = -305 + 0.17 * T
    dG0_TiO   = -1029 + 0.15 * T
    dG0_BeO   = -1216 + 0.2 * T
    dG0_BaO   = -1114 + 0.21 * T
    dG0_SrO   = -1194 + 0.2 * T
    dG0_NbO   =  -825 + 0.17 * T

    # M2O3:
    dG0_Al2O3 = -1125 + 0.22 * T
    dG0_Cr2O3 = -740 + 0.16  * T
    dG0_V2O3  = -802 + 0.16  * T
    dG0_Fe2O3 = -543 + 0.17  * T
    dG0_Mn2O3 = -636 + 0.17  * T
    dG0_Ti2O3 =  -1001 + 0.17 * T
    dG0_Sc2O3 = -1268 + 0.19 * T

    # M2O4:
    dG0_VO2  = -706 + 0.16 * T  # Warning V2O4.

    # M2O:
    dG0_Cu2O  = -337  + 0.14 * T
    dG0_Li2O  = -1205 + 0.27 * T
    dG0_Na2O  = -843  + 0.28 * T
    dG0_K2O   = -724  + 0.28 * T

    # M2O5:
    dG0_V2O5  = -579 + 0.13 * T   # Warning (liquid phase).
    dG0_Nb2O5 = -755 + 0.17 * T

    # M3O5:
    dG0_Ti3O5 = -974 + 0.17 * T

    # M3O4:
    dG0_Mn3O4 = -692 + 0.17 * T
    dG0_Fe3O4 = -551 + 0.15 * T

    # MO3:
    dG0_MoO3   = -493 + 0.16 * T
    dG0_WO3    = -556 + 0.16 * T

    # Mixed:
    dG0_LiAlO2  = -1.19217958e+03 + 2.28562409e-01 * T
    dG0_MgSiO3  = -1.03001604e+03 + 1.90275343e-01 * T
    dG0_Mg2SiO4 = -1.09441571e+03 + 2.23360614e-01 * T
    dG0_MgTiO3  = -1.05452033e+03 + 2.17230059e-01 * T
    dG0_MgTi2O5 = -1.00604335e+03 + 1.97454003e-01 * T
    dG0_Al2MgO4 = -1.15485171e+03 + 2.26578897e-01 * T
    dG0_NaAlO2  = -1.13500104e+03 + 2.21646594e-01 * T
    dG0_Al2SiO5 = -1.03270670e+03 + 1.87057775e-01 * T
    dG0_Li2SiO3 = -1.09806841e+03 + 2.06274953e-01 * T
    dG0_Li2Si2O5 = -1.02530846e+03 + 1.96245135e-01  * T
    dG0_Li2TiO3 = -1.11159218e+03 + 2.03484859e-01 * T
    dG0_Na2SiO3 = -1.03879943e+03 + 2.10765627e-01 * T
    dG0_Na2Si2O5 = -9.84712696e+02 + 1.87778973e-01 * T
    dG0_K2SiO3  = -1.03000574e+03 + 2.05438329e-01 * T
    dG0_LiFeO2  =  -708 + 0.1656  * T
    dG0_LiFe5O8 = (-2341 + 0.6764 * T) * 2/8
    dG0_Na2ZrO3 = (-1676.27 + 0.348 * T) * 2/3
    dG0_LiCoO2  = -615.1 + 211.8e-03 * T

    # Gas phase molecules:
    dG0_CO = -2.23623129e+02 -1.75809004e-01 * T
    dG0_CO2 = -396 + 0 * T

    results = {
               'temperature' : T,
               'TiO2_rut'    : dG0_TiO2_rut,
               'TiO2'        : dG0_TiO2_rut,
               'SiO2_qu'     : dG0_SiO2_qu,
               'SiO2'        : dG0_SiO2_qu,
               'SnO2'        : dG0_SnO2,
               'MnO2'        : dG0_MnO2,
               'ZrO2'        : dG0_ZrO2,
               'MoO2'        : dG0_MoO2,
               'WO2'         : dG0_WO2,
               'KO2'         : dG0_KO2,
               'NbO2'        : dG0_NbO2,
               'CaO'         : dG0_CaO,
               'MnO'         : dG0_MnO,
               'NiO'         : dG0_NiO,
               'FeO'         : dG0_FeO,
               'CoO'         : dG0_CoO,
               'VO'          : dG0_VO,
               'CuO'         : dG0_CuO,
               'TiO'         : dG0_TiO,
               'BaO'         : dG0_BaO,
               'MgO'         : dG0_MgO,
               'BeO'         : dG0_BeO,
               'SrO'         : dG0_SrO,
               'NbO'         : dG0_NbO,
               'Al2O3'       : dG0_Al2O3,
               'Cr2O3'       : dG0_Cr2O3,
               'V2O3'        : dG0_V2O3,
               'Fe2O3'       : dG0_Fe2O3,
               'Mn2O3'       : dG0_Mn2O3,
               'Ti2O3'       : dG0_Ti2O3,
               'Sc2O3'       : dG0_Sc2O3,
               'VO2'         : dG0_VO2,
               'Cu2O'        : dG0_Cu2O,
               'Li2O'        : dG0_Li2O,
               'Na2O'        : dG0_Na2O,
               'V2O5'        : dG0_V2O5,
               'Nb2O5'       : dG0_Nb2O5,
               'Ti3O5'       : dG0_Ti3O5,
               'Mn3O4'       : dG0_Mn3O4,
               'Fe3O4'       : dG0_Fe3O4,
               'K2O'         : dG0_K2O,
               'MoO3'        : dG0_MoO3,
               'WO3'         : dG0_WO3,
               'LiAlO2'      : dG0_LiAlO2,
               'Li2SiO3'     : dG0_Li2SiO3,
               'Li2Si2O5'    : dG0_Li2Si2O5,
               'LiCoO2'      : dG0_LiCoO2,
               'Li2TiO3'     : dG0_Li2TiO3,
               'LiFeO2'      : dG0_LiFeO2,
               'LiFe5O8'     : dG0_LiFe5O8,
               'MgSiO3'      : dG0_MgSiO3,
               'Mg2SiO4'     : dG0_Mg2SiO4,
               'MgTiO3'      : dG0_MgTiO3,
               'MgTi2O5'     : dG0_MgTi2O5,
               'Al2MgO4'     : dG0_Al2MgO4,
               'NaAlO2'      : dG0_NaAlO2,
               'Na2SiO3'     : dG0_Na2SiO3,
               'Na2Si2O5'    : dG0_Na2Si2O5,
               'Al2SiO5'     : dG0_Al2SiO5,
               'K2SiO3'      : dG0_K2SiO3,
               'Na2ZrO3'     : dG0_Na2ZrO3,
               'CO'          : dG0_CO,
               'CO2'         : dG0_CO2,
               }
    return results
