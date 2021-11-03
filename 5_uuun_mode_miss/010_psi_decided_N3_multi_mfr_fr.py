"""
デトネーション波の枚数 N を振るコード
全角度 360 deg. -> 360/N deg.       Noting that the code outputs only parameters in a range of 0 - 360/N deg. 
流量 1/N -> fillheight 1/N
ΔA/Δtheta N 依存
"""

# Suppress the 'divide by zero' warning encountered for some data points in the trig calls
# so that we can keep the code as similar as possible to the MATLAB equivalent
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
print('low order parametric analysis?')

import time
start = time.time()

import cantera as ct
import math
from sdtoolbox.thermo import soundspeed_fr, soundspeed_eq
from sdtoolbox.postshock import CJspeed, PostShock_eq, PostShock_fr
from scipy.interpolate import pchip # creates a PCHIP interpolation function that can then be given a query point
from scipy.optimize import fminbound
from sys import exit
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os

start = time.time()

#===== デトネーション波の枚数 N
N_deto = 3.
#===== 管摩擦係数
Cf = 0.008
#===== psi_sl
psi_sl = 30. /360. * 2. * np.pi
#=====plenum temperture
T_exp = 300.0 # Temperature [K]
#=====calculation for gas
ER = 1.
coef_C2H4 = 1
coef_O2 = 3
mol_ratio = coef_C2H4/coef_O2*ER
fuel = 'C2H4'
oxygen = 'O2'
mech = 'gri30_highT.cti'
# ===== folder
path = "Test_sample_Ndeto_" + str(int(N_deto)) + "_fr"
os.makedirs(path, exist_ok=True)
os.chdir(path)
path1 = "Cf" + str(Cf) + "_psi_sl" + str(psi_sl/2./np.pi*360.) + "deg_" + "Tple" + str(T_exp) + "K"
os.makedirs(path1, exist_ok=True)
os.chdir("../")
#===========================================================================
#======================= sawada example====================================
#===========================================================================
#=====geometory
rm = 35. * 10. ** (-3.) # [m]
delta_a = 8. * 10. ** (-3.) # [m]
Acom = 2. * np.pi * rm * delta_a # combustor cross section area [m2]
Ainj = 120. * 2. * (0.0005**2.*np.pi) # incjection area [m2]
r_in = rm - delta_a / 2.
r_out = rm + delta_a / 2.

def func_atm(atm):
    #=====calculation for gas
    gas = ct.Solution(mech)
    q = fuel + ':' + str(mol_ratio) + ', ' + oxygen + ':1.0';
    # ===== error
    eps = 10 ** (-6.)
    num_angle = 3600 # delta_theta = 360. / num_angle
    #===========================================================================
    #=========================injection (state1) ===============================
    #===========================================================================

    #===== from Mizener
    T_ple = T_exp
    P_ple = atm * 101300

    # ===========================================
    # Plenum:
    # ===========================================
    ###### Set plenum condition
    gas.TPX = T_ple,P_ple,q
    s0 = gas.entropy_mass
    h0 = gas.enthalpy_mass
    rho_ple = gas.density
    R_ple = ct.gas_constant / gas.mean_molecular_weight
    a_fr_ple = soundspeed_fr(gas)
    gamma_ple = a_fr_ple**2.*rho_ple/P_ple
    print('plenum section')
    print('     gamma_ple', gamma_ple, '-')
    print('     density', gas.density, 'kg/m3')
    print('     pressure', P_ple, 'Pa','////', gas.P/1000000, 'MPa')
    print('     temperature', T_ple, 'K')
    print('#=====================================================================================================================')
    print('#=====================================================================================================================')

    ###### ==============================================================================================================
    ###### critical values through injector holes
    ###### ==============================================================================================================
    Pcr = P_ple * (2./(gamma_ple + 1.)) ** ((gamma_ple)/(gamma_ple-1.))                  # Eq.(3)
    Pcom = Pcr
    Tcom = T_ple * (Pcom/P_ple) ** ((gamma_ple-1.)/gamma_ple)                                                # Eq.(4)
    Vinj = ((2. * gamma_ple)/(gamma_ple - 1.) * R_ple * T_ple * (1.-(Pcom/P_ple)**((gamma_ple-1.)/gamma_ple))) ** (1./2.)

    ############### set gas (Tcom, Pcom, q)
    gas.TPX = Tcom, Pcom, q
    h1 = gas.enthalpy_mass
    s1 = gas.entropy_mass
    rho1 = gas.density
    w1 = gas.mean_molecular_weight
    R1 = ct.gas_constant/w1
    a1_fr = soundspeed_fr(gas)
    gamma1_fr =  a1_fr**2.*rho1/Pcom
    mfr_inj = rho1 * Ainj * Vinj
    print('Layer detonation computation for '+mech+' with composition '+q)
    print('State 1 - Initial state of reacting layer')
    print('   Pressure '+str(gas.P)+' (Pa)')
    print('   Temperature '+str(gas.T)+' (K)')
    print('   Density '+str(rho1)+' (kg/m3)')
    print('   Sound speed (frozen) '+str(a1_fr)+' (m/s)')
    print('   Enthalpy '+str(h1)+' (J/kg)')
    print('   Entropy '+str(s1)+' (J/kg K)')
    print('   gamma (frozen) '+str(gamma1_fr)+' ')
    print('   mfr_inj ', mfr_inj)

    ##
    # Find CJ speed
    U_CJ = CJspeed(Pcom, Tcom, q, mech)
    # Evaluate CJ gas state
    gas = PostShock_eq(U_CJ,Pcom, Tcom, q, mech)
    x2 = gas.X
    P2 = gas.P
    T2 = gas.T
    rho2 = gas.density
    a2_eq = soundspeed_eq(gas)
    s2 = gas.entropy_mass
    h2 = gas.enthalpy_mass
    w2 = rho1*U_CJ/rho2
    u2 = U_CJ-w2
    gamma2 = a2_eq**2.*rho2/P2

    print('State 2 - CJ ')
    print('   CJ speed '+str(U_CJ)+' (m/s)')
    print('   Pressure '+str(P2)+' (Pa)')
    print('   Temperature '+str(T2)+' (K)')
    print('   Density '+str(rho2)+' (kg/m3)')
    print('   Enthalpy '+str(h2)+' (J/kg)')
    print('   Entropy '+str(s2)+' (J/kg K)')
    print('   w2 (wave frame) '+str(w2)+' (m/s)')
    print('   u2 (lab frame) '+str(u2)+' (m/s)')
    print('   a2 (equilibrium) '+str(a2_eq)+' (m/s)')
    U2 = w2
    ##
    gas.TPX = T2,P2,x2
    print("VinjVinjVinjVinjVinjVinjVinjVinjVinjVinjVinjVinj", Vinj)
    h3_U2 = (h2 + U2**2./2.) # SDT
    print("h3_U2h3_U2h3_U2h3_U2h3_U2", h3_U2)
    s3 = gas.entropy_mass

    ##### ================================================================================================================ waveheight : 1/N
    ##### waveheight3 : from pre_det
    dt_com = 2 * np.pi * rm / U_CJ
    waveheight = Vinj * dt_com
    waveheight = waveheight / N_deto
    print("waveheight ->", waveheight)

    ##### calculate maximum value of tau
    tau0 = 1./2. * Cf * rho2 * abs(U2 - U_CJ) * (U2 - U_CJ)
    tau_psi_sl0 = 1./2. * Cf * rho2 * abs(U2 * math.cos(psi_sl) - U_CJ) * (U2* math.cos(psi_sl) - U_CJ)
    print("tau0tau0tau0tau0tau0tau0tau0tau0tau0", tau0)

    #===========================================================================
    #======================= angle_lists====================================
    #===========================================================================
    dtheta = 2. * np.pi / num_angle
    theta_lists_origin = [0] * (num_angle)
    for i in range(num_angle):
        theta_lists_origin[i] = dtheta * i
    theta_lists = theta_lists_origin

    #===========================================================================
    #============================== wave height ================================
    #===========================================================================
    A_sonic = waveheight * delta_a
    def func_dMwc(Mwc,p):
        gas.SPX = s3, p, x2
        v2 = 2.0*(h3_U2 - gas.enthalpy_mass)
        v = math.sqrt(v2)
        dMwc = abs(Mwc - v/soundspeed_fr(gas))
        return dMwc

    def func_delta_Eq16(A, Mwc):
        Pwc_a = 0.001 * 101300
        Pwc_b = 0.002 * 101300
        dMwc_b = func_dMwc(Mwc, Pwc_b)
        dMwc_a = func_dMwc(Mwc, Pwc_a)
        while abs(dMwc_b) > eps:
            Pwc_s = (Pwc_a * dMwc_b - Pwc_b * dMwc_a)/(dMwc_b - dMwc_a)
            Pwc_a, Pwc_b = Pwc_b, Pwc_s
            dMwc_a = dMwc_b
            dMwc_b = func_dMwc(Mwc, Pwc_b)
        gas.SPX = s3, Pwc_b,x2
        Twc = gas.T
        Rwc = ct.gas_constant / gas.mean_molecular_weight
        rhowc = gas.density
        hwc = gas.enthalpy_mass
        a_frwc = soundspeed_fr(gas)
        Wwc = Mwc * a_frwc
        hwc_total = hwc + Wwc ** 2. / 2.
        gammawc = a_frwc**2.*rhowc/Pwc_b
        delta_Eq16 = (A/A_sonic) ** 2. - (1./Mwc) ** 2. * ((2. + (gammawc-1.)*Mwc**2.)/(gammawc+1.)) ** ((gammawc+1.)/(gammawc-1.))
        return delta_Eq16, Pwc_b, gammawc, Twc, Rwc, a_frwc, Wwc, rhowc, hwc, hwc_total

    ######## C section
    Mwc_lists = []
    Pwc_lists = []
    Twc_lists = []
    rhowc_lists = []
    a_frwc_lists = []
    gammawc_lists = []
    # Rwc_lists = []
    Vwc_lists = []
    Vwc_lists_psi_sl = []
    dKE_lists = []

    ######## add zeta
    zeta_lists = []

    dA_lists = []

    dTz_in_lists = []
    dTz_out_lists = []

    ######## all section
    ### tau
    tau_lists = []
    tau_lists_psi_sl = []
    tau_lists_half_psi_sl = []
    ### dftheta
    dftheta_in_lists = []
    dftheta_out_lists = []
    dftheta_in_lists1 = []
    dftheta_out_lists1 = []
    ### dTz
    dTz_in_lists = []
    dTz_out_lists = []
    dTz_in1_lists = []
    dTz_out1_lists = []

    ######### negative section
    ### tau_N
    tau_lists_N = []
    tau_lists_psi_sl_N = []
    tau_lists_half_psi_sl_N = []
    ### dftheta_N
    dftheta_in_lists_N = []
    dftheta_out_lists_N = []
    dftheta_in_lists1_N = []
    dftheta_out_lists1_N = []
    ### dTz_N
    dTz_in_lists_N = []
    dTz_out_lists_N = []
    dTz_in1_lists_N = []
    dTz_out1_lists_N = []

    ######## inj section
    Pinj_lists = []
    Tinj_lists = []
    Vinj_lists = []
    rhoinj_lists = []
    a_frinj_lists = []
    gammainj_lists = []
    # Rinj_lists = []

    ######## add zeta
    zeta_lists = []

    Pwc = P_ple*1.0001
    theta_num = 0

    while Pwc >= P_ple:
        if theta_num < (num_angle/N_deto):
            theta = theta_lists[theta_num]
            zeta = waveheight + rm * theta * math.tan(psi_sl)
            A = delta_a * zeta
            Mwc_a = 2.
            Mwc_b = 2.4
            delta_Eq16_b = func_delta_Eq16(A, Mwc_b)[0]
            delta_Eq16_a = func_delta_Eq16(A, Mwc_a)[0]
            while abs(delta_Eq16_b) > eps:
                Mwc_s = (Mwc_a*delta_Eq16_b-Mwc_b*delta_Eq16_a)/(delta_Eq16_b-delta_Eq16_a)
                Mwc_a,Mwc_b = Mwc_b, Mwc_s
                delta_Eq16_a = delta_Eq16_b
                delta_Eq16_b, Pwc_b, gammawc_b, Twc_b, Rwc_b, awc_b, Wwc_b, rhowc_b, hwc_b, hwc_total_b = func_delta_Eq16(A, Mwc_b)
            ######## C section
            gas.SPX = s3, Pwc_b, x2
            Pwc = Pwc_b
            Mwc = Mwc_b
            rhowc = gas.density
            Mwc_lists.append(Mwc_b)
            Pwc_lists.append(Pwc_b)
            Twc_lists.append(gas.T)
            rhowc_lists.append(rhowc)
            # Rwc_lists.append(ct.gas_constant / gas.mean_molecular_weight)
            a_frwc = soundspeed_fr(gas)
            a_frwc_lists.append(a_frwc)
            gammawc_lists.append(a_frwc**2*rhowc/Pwc_b)
            Vwc = Mwc * a_frwc - U_CJ
            Vwc_psi_sl = Mwc * a_frwc * math.cos(psi_sl) - U_CJ
            Vwc_half_psi_sl = Mwc * a_frwc * math.cos(psi_sl/2.) - U_CJ
            Vwc_lists_psi_sl.append(Vwc_psi_sl)
            Vwc_lists.append(Vwc)
            
            # ===== friction : all of kinetic energy contribute to the azimuthal component
            ### 流路断面積は実効断面積がというだけの話であって濡れ面積と関係ない
            dA_sum = sum(dA_lists)
            # dA_bottom = delta_a * dtheta * rm * phai
            # dA_fresh = 1./2. * rm * (theta - theta_in) * waveheight * (theta-theta_in)/(2.*np.pi-theta_in)
            dA = waveheight * rm * theta + 1./2. * rm*theta * (rm*theta) * math.tan(psi_sl) - dA_sum    # + dA_bottom
            dA_lists.append(dA)
            dA_in = r_in / rm * dA
            dA_out = rm / r_out * dA

            # ===== dKE
            dKE = rhowc * abs(Vwc) * Vwc / 2.
            dKE_lists.append(dKE)

            # ===== tau
            tau = Cf * rhowc * abs(Vwc) * Vwc / 2.
            tau_psi_sl = Cf * rhowc * abs(Vwc_psi_sl) * Vwc_psi_sl / 2.

            # ===== dftheta
            dftheta_in = dA_in * tau
            dftheta_out = dA_out * tau
            dftheta_in1 = dA_in * tau_psi_sl
            dftheta_out1 = dA_out * tau_psi_sl
            
            # ===== dTz
            dTz_in = r_in * dftheta_in
            dTz_out = r_out * dftheta_out
            dTz_in1 = r_in * dftheta_in1
            dTz_out1 = r_out * dftheta_out1
            
            ##### all section
            ### tau
            tau_lists.append(tau)
            tau_lists_psi_sl.append(tau_psi_sl)
            ### dftheta
            dftheta_in_lists.append(dftheta_in)
            dftheta_out_lists.append(dftheta_out)
            dftheta_in_lists1.append(dftheta_in1)
            dftheta_out_lists1.append(dftheta_out1)
            ### dTz
            dTz_in_lists.append(dTz_in)
            dTz_out_lists.append(dTz_out)
            dTz_in1_lists.append(dTz_in1)
            dTz_out1_lists.append(dTz_out1)

            if tau < 0.2 * tau0:
                ##### negative Vwc section
                ### tau
                tau_lists_N.append(tau)
                ### dftheta (no use)
                ### dTz
                dTz_in_lists_N.append(dTz_in)
                dTz_out_lists_N.append(dTz_out)
            if tau_psi_sl < 0.2 * tau_psi_sl0:
                ##### negative Vwc section
                ### tau
                tau_lists_psi_sl_N.append(tau_psi_sl)
                ### dftheta (no use)
                ### dTz
                dTz_in1_lists_N.append(dTz_in1)
                dTz_out1_lists_N.append(dTz_out1)

            ######## add zeta
            zeta_lists.append(zeta)
            ######## next theta_num
            theta_num = theta_num + 1

        else:
            print("Pc at all points is under P0inj")
            theta_in = 0
            theta_ch = 0
            break

    # ===== theta_in
    theta_in = theta_lists[theta_num-1]
    theta_in_num = theta_num - 1
    print("/////////////////////////////////////////////////////////////////////////////////// theta_in -> ", theta_in * 360 / 2/ np.pi)

    while Pwc >= Pcr:
        if theta_num < (num_angle/N_deto):
            theta = theta_lists[theta_num]
            zeta = waveheight + rm * theta * math.tan(psi_sl) - waveheight * ((theta-theta_in)/(2*np.pi-theta_in))
            A = delta_a * zeta
            Mwc_a = 3.1
            Mwc_b = 3.2
            delta_Eq16_b = func_delta_Eq16(A, Mwc_b)[0]
            delta_Eq16_a = func_delta_Eq16(A, Mwc_a)[0]
            while abs(delta_Eq16_b) > eps:
                Mwc_s = (Mwc_a*delta_Eq16_b-Mwc_b*delta_Eq16_a)/(delta_Eq16_b-delta_Eq16_a)
                Mwc_a,Mwc_b = Mwc_b, Mwc_s
                delta_Eq16_a = delta_Eq16_b
                delta_Eq16_b, Pwc_b, gammawc_b, Twc_b, Rwc_b, awc_b, Wwc_b, rhowc_b, hwc_b, hwc_total_b = func_delta_Eq16(A, Mwc_b)
            ######## C section
            gas.SPX = s3, Pwc_b, x2
            Pwc = Pwc_b
            Mwc = Mwc_b
            rhowc = gas.density
            Mwc_lists.append(Mwc_b)
            Pwc_lists.append(Pwc_b)
            Twc_lists.append(gas.T)
            rhowc_lists.append(rhowc)
            # Rwc_lists.append(ct.gas_constant / gas.mean_molecular_weight)
            a_frwc = soundspeed_fr(gas)
            a_frwc_lists.append(a_frwc)
            gammawc_lists.append(a_frwc**2*rhowc/Pwc_b)
            Vwc = Mwc * a_frwc - U_CJ
            Vwc_psi_sl = Mwc * a_frwc * math.cos(psi_sl) - U_CJ
            Vwc_half_psi_sl = Mwc * a_frwc * math.cos(psi_sl/2.) - U_CJ
            Vwc_lists_psi_sl.append(Vwc_psi_sl)
            Vwc_lists.append(Vwc)

            # ===== friction : all of kinetic energy contribute to the azimuthal component
            ### 流路断面積は実効断面積がというだけの話であって濡れ面積と関係ない
            dA_sum = sum(dA_lists)
            # dA_bottom = delta_a * dtheta * rm * phai
            dA_fresh = 1./2. * rm * (theta - theta_in) * waveheight * (theta-theta_in)/(2.*np.pi-theta_in)
            dA = waveheight * rm * theta + 1./2. * rm*theta * (rm*theta) * math.tan(psi_sl) - dA_sum - dA_fresh       # + dA_bottom
            dA_lists.append(dA)
            dA_in = r_in / rm * dA
            dA_out = rm / r_out * dA

            # ===== dKE
            dKE = rhowc * abs(Vwc) * Vwc / 2.
            dKE_lists.append(dKE)

            # ===== tau
            tau = Cf * rhowc * abs(Vwc) * Vwc / 2.
            tau_psi_sl = Cf * rhowc * abs(Vwc_psi_sl) * Vwc_psi_sl / 2.

            # ===== dftheta
            dftheta_in = dA_in * tau
            dftheta_out = dA_out * tau
            dftheta_in1 = dA_in * tau_psi_sl
            dftheta_out1 = dA_out * tau_psi_sl
            
            # ===== dTz
            dTz_in = r_in * dftheta_in
            dTz_out = r_out * dftheta_out
            dTz_in1 = r_in * dftheta_in1
            dTz_out1 = r_out * dftheta_out1
            
            ##### all section
            ### tau
            tau_lists.append(tau)
            tau_lists_psi_sl.append(tau_psi_sl)
            ### dftheta
            dftheta_in_lists.append(dftheta_in)
            dftheta_out_lists.append(dftheta_out)
            dftheta_in_lists1.append(dftheta_in1)
            dftheta_out_lists1.append(dftheta_out1)
            ### dTz
            dTz_in_lists.append(dTz_in)
            dTz_out_lists.append(dTz_out)
            dTz_in1_lists.append(dTz_in1)
            dTz_out1_lists.append(dTz_out1)

            if tau < 0.2 * tau0:
                ##### negative Vwc section
                ### tau
                tau_lists_N.append(tau)
                ### dftheta (no use)
                ### dTz
                dTz_in_lists_N.append(dTz_in)
                dTz_out_lists_N.append(dTz_out)
            if tau_psi_sl < 0.2 * tau_psi_sl0:
                ##### negative Vwc section
                ### tau
                tau_lists_psi_sl_N.append(tau_psi_sl)
                ### dftheta (no use)
                ### dTz
                dTz_in1_lists_N.append(dTz_in1)
                dTz_out1_lists_N.append(dTz_out1)

            ######## add zeta
            zeta_lists.append(zeta)
            ######## next theta_num
            theta_num = theta_num + 1
        else:
            print("Pc at all points is under Pcr")
            theta_ch = 0
            break

    # ===== theta_ch
    theta_ch = theta_lists[theta_num-1]
    theta_ch_num = theta_num - 1
    print("/////////////////////////////////////////////////////////////////////////////////// theta_ch -> ", theta_ch * 360 / 2/ np.pi)

    while theta_num < (num_angle/N_deto):
        theta = theta_lists[theta_num]
        zeta = waveheight + rm * theta * math.tan(psi_sl) - waveheight * ((theta-theta_in)/(2.*np.pi-theta_in))
        A = delta_a * zeta
        Mwc_a = 3.6
        Mwc_b = 3.7
        delta_Eq16_b = func_delta_Eq16(A, Mwc_b)[0]
        delta_Eq16_a = func_delta_Eq16(A, Mwc_a)[0]
        while abs(delta_Eq16_b) > eps:
            Mwc_s = (Mwc_a*delta_Eq16_b-Mwc_b*delta_Eq16_a)/(delta_Eq16_b-delta_Eq16_a)
            Mwc_a,Mwc_b = Mwc_b, Mwc_s
            delta_Eq16_a = delta_Eq16_b
            delta_Eq16_b, Pwc_b, gammawc_b, Twc_b, Rwc_b, awc_b, Wwc_b, rhowc_b, hwc, hwc_total = func_delta_Eq16(A, Mwc_b)
        ######## C section
        gas.SPX = s3, Pwc_b, x2
        Pwc = Pwc_b
        Mwc = Mwc_b
        rhowc = gas.density
        Mwc_lists.append(Mwc_b)
        Pwc_lists.append(Pwc_b)
        Twc_lists.append(gas.T)
        rhowc_lists.append(rhowc)
        # Rwc_lists.append(ct.gas_constant / gas.mean_molecular_weight)
        a_frwc = soundspeed_fr(gas)
        a_frwc_lists.append(a_frwc)
        gammawc_lists.append(a_frwc**2.*rhowc/Pwc_b)
        Vwc = Mwc * a_frwc - U_CJ
        Vwc_psi_sl = Mwc * a_frwc * math.cos(psi_sl) - U_CJ
        Vwc_half_psi_sl = Mwc * a_frwc * math.cos(psi_sl/2.) - U_CJ
        Vwc_lists_psi_sl.append(Vwc_psi_sl)
        Vwc_lists.append(Vwc)

        # ===== friction : all of kinetic energy contribute to the azimuthal component
        ### 流路断面積は実効断面積がというだけの話であって濡れ面積と関係ない
        dA_sum = sum(dA_lists)
        # dA_bottom = delta_a * dtheta * rm * phai
        dA_fresh = 1./2. * rm * (theta - theta_in) * waveheight * (theta-theta_in)/(2.*np.pi-theta_in)
        dA = waveheight * rm * theta + 1./2. * rm*theta * (rm*theta) * math.tan(psi_sl) - dA_sum - dA_fresh       # + dA_bottom
        dA_lists.append(dA)
        dA_in = r_in / rm * dA
        dA_out = rm / r_out * dA

        # ===== dKE
        dKE = rhowc * abs(Vwc) * Vwc / 2.
        dKE_lists.append(dKE)

        # ===== tau
        tau = Cf * rhowc * abs(Vwc) * Vwc / 2.
        tau_psi_sl = Cf * rhowc * abs(Vwc_psi_sl) * Vwc_psi_sl / 2.
        tau_half_psi_sl = Cf * rhowc * abs(Vwc_half_psi_sl) * Vwc_half_psi_sl / 2.

        # ===== dftheta
        dftheta_in = dA_in * tau
        dftheta_out = dA_out * tau
        dftheta_in1 = dA_in * tau_psi_sl
        dftheta_out1 = dA_out * tau_psi_sl
        
        # ===== dTz
        dTz_in = r_in * dftheta_in
        dTz_out = r_out * dftheta_out
        dTz_in1 = r_in * dftheta_in1
        dTz_out1 = r_out * dftheta_out1
        
        ##### all section
        ### tau
        tau_lists.append(tau)
        tau_lists_psi_sl.append(tau_psi_sl)
        ### dftheta
        dftheta_in_lists.append(dftheta_in)
        dftheta_out_lists.append(dftheta_out)
        dftheta_in_lists1.append(dftheta_in1)
        dftheta_out_lists1.append(dftheta_out1)
        ### dTz
        dTz_in_lists.append(dTz_in)
        dTz_out_lists.append(dTz_out)
        dTz_in1_lists.append(dTz_in1)
        dTz_out1_lists.append(dTz_out1)

        if tau < 0.2 * tau0:
            ##### negative Vwc section
            ### tau
            tau_lists_N.append(tau)
            ### dftheta (no use)
            ### dTz
            dTz_in_lists_N.append(dTz_in)
            dTz_out_lists_N.append(dTz_out)
        if tau_psi_sl < 0.2 * tau_psi_sl0:
            ##### negative Vwc section
            ### tau
            tau_lists_psi_sl_N.append(tau_psi_sl)
            ### dftheta (no use)
            ### dTz
            dTz_in1_lists_N.append(dTz_in1)
            dTz_out1_lists_N.append(dTz_out1)

        ######## add zeta
        zeta_lists.append(zeta)
        ######## next theta_num
        theta_num = theta_num + 1

    ###### all section
    # ===== mean_tau
    mean_tau = sum(tau_lists)/len(tau_lists)
    print("===================================================================================")
    print("================mean_tau", mean_tau)
    print("===================================================================================")

    # ===== mean_tau_psi_sl
    mean_tau_psi_sl = sum(tau_lists_psi_sl)/len(tau_lists_psi_sl)
    print("===================================================================================")
    print("================mean_tau", mean_tau_psi_sl)
    print("===================================================================================")

    print("")
    print("")
    print("")

    ###### negative section
    # ===== mean_tau_N
    mean_tau_N = sum(tau_lists_N)/len(tau_lists_N)
    print("===================================================================================")
    print("================mean_tau_N", mean_tau_N)
    print("===================================================================================")

    # ===== mean_tau_psi_sl_N
    mean_tau_psi_sl_N = sum(tau_lists_psi_sl_N)/len(tau_lists_psi_sl_N)
    print("===================================================================================")
    print("================mean_tau_psi_sl_N", mean_tau_psi_sl_N)
    print("===================================================================================")

    print("")
    print("")
    print("")


    ###### all section
    # ===== Tz
    Tz = sum(dTz_in_lists) + sum(dTz_out_lists)
    print("===================================================================================")
    print("================Tz", Tz)
    print("===================================================================================")

    # ===== Tz1
    Tz_psi_sl = sum(dTz_in1_lists) + sum(dTz_out1_lists)
    print("===================================================================================")
    print("================Tz_psi_sl", Tz_psi_sl)
    print("===================================================================================")

    print("")
    print("")
    print("")

    ###### negative section
    # ===== Tz1_N
    Tz_N = sum(dTz_in_lists_N) + sum(dTz_out_lists_N)
    print("===================================================================================")
    print("================Tz_N", Tz_N)
    print("===================================================================================")

    # ===== Tz2_N
    Tz_psi_sl_N = sum(dTz_in1_lists_N) + sum(dTz_out1_lists_N)
    print("===================================================================================")
    print("================Tz_psi_sl_N", Tz_psi_sl_N)
    print("===================================================================================")










    ######## C section
    Mwc_lists_output = Mwc_lists * int(N_deto)
    Pwc_lists_output = Pwc_lists * int(N_deto)
    Twc_lists_output = Twc_lists * int(N_deto)
    rhowc_lists_output = rhowc_lists * int(N_deto)
    a_frwc_lists_output = a_frwc_lists * int(N_deto)
    gammawc_lists_output = gammawc_lists * int(N_deto)
    # Rwc_lists_output = Rwc_lists * int(N_deto)
    Vwc_lists_output = Vwc_lists * int(N_deto)
    Vwc_lists_psi_sl_output = Vwc_lists_psi_sl * int(N_deto)
    dKE_lists_output = dKE_lists * int(N_deto)
    ######## add zeta
    zeta_lists_output = zeta_lists * int(N_deto)
    ######## all section
    ### tau
    tau_lists_output = tau_lists * int(N_deto)
    tau_lists_psi_sl_output = tau_lists_psi_sl * int(N_deto)
    ######## inj section
    Pinj_lists_output = Pinj_lists * int(N_deto)
    Tinj_lists_output = Tinj_lists * int(N_deto)
    Vinj_lists_output = Vinj_lists * int(N_deto)
    rhoinj_lists_output = rhoinj_lists * int(N_deto)
    a_frinj_lists_output = a_frinj_lists * int(N_deto)
    gammainj_lists_output = gammainj_lists * int(N_deto)
    # Rinj_lists_output = Rinj_lists * int(N_deto)


    os.chdir(path)
    os.chdir(path1)
    df_lists = pd.DataFrame(index=[],columns=[])
    ###################################################inlet calculation
    df_lists['theta [rad]'] = theta_lists
    df_lists['Mwc [-]'] = Mwc_lists_output
    df_lists['Pwc [Pa]'] = Pwc_lists_output
    df_lists['Twc [K]'] = Twc_lists_output
    df_lists['rhowc [kg/m3]'] = rhowc_lists_output
    df_lists['a_frwc [m/s]'] = a_frwc_lists_output
    df_lists['gammawc [-]'] = gammawc_lists_output
    df_lists['Vwc [m/s]'] = Vwc_lists_output
    df_lists['Vwc_psi_sl [m/s]'] = Vwc_lists_psi_sl_output
    # df_lists['Rwc []'] = Rwc_lists_output
    df_lists['tau [-]'] = tau_lists_output
    df_lists['tau_psi_sl [-]'] = tau_lists_psi_sl_output
    df_lists['zeta_inlet [m]'] = zeta_lists_output
    df_lists['dKE'] = dKE_lists_output

    df_lists.to_csv(str(atm) + "_theta.csv",index=False)

    df_value = pd.DataFrame(index=[],columns=['result'])
    df_value['mass flow rate [g/s]'] = [mfr_inj]
    df_value["ER [-]"] = ER
    df_value['P plenum [Pa]'] = [P_ple]
    df_value['T plenum [K]'] = [T_ple]
    df_value['pressure state1(Pcom) [Pa]'] = [Pcom]
    df_value['temperature state1(Tcom) [K]'] = [Tcom]
    df_value['density state1(rho1) [kg/m3]'] = [rho1]
    # df_value['soud speed state1 [m/s]'] = [a1_fr]
    # df_value['enthalpy state1 [J/kg]'] = [h1]
    # df_value['entropy state1 [J/kg/K]'] = [s1]
    # df_value['gamma state1 [-]'] = [gamma1_fr]
    df_value['CJ speed [m/s]'] = [U_CJ]
    # df_value['pressure state2 [Pa]'] = [P2]
    # df_value['temperature state2 [K]'] = [T2]
    # df_value['density state2 [kg/m3]'] = [rho2]
    # df_value['enthalpy state2 [J/kg]'] = [h2]
    # df_value['entropy state2 [J/kg/K]'] = [s2]
    # df_value['w2(wave frame) [m/s]'] = [w2]
    # df_value['u2(lab frame) [m/s]'] = [u2]
    # df_value['sound speed state2 [m/s]'] = [a2_eq]
    # df_value['wave angle [deg]'] = [angle_det]
    # df_value['downstream flow angle [deg]'] = [theta_det]
    # df_value['triple point speed [m/s]'] = [u_tp]
    # df_value['downstream speed(wave frame)(x) [m/s]'] = [U2]
    # df_value['speed tangential to deto(z) [m/s]'] = [vt_det]
    # df_value['Mc (upstream of oblique shock) [-]'] = [Mc]
    # df_value['M3 [-]'] = [M3]
    # df_value['psi_sh [deg]'] = [psi_sh /2/ np.pi* 360]
    # df_value['psi_sl [deg]'] = [psi_sl/2/np.pi *360]
    # df_value['theta_in [deg]'] = [theta_in/2/np.pi*360]
    # df_value['theta_ch [deg]'] = [theta_ch/2/np.pi*360]
    # df_value['waveheight [m]'] = [waveheight]
    # df_value['fz [N]'] = [fz]
    # df_value['ftheta [N]'] = [ftheta]
    # df_value['force ratio [-]'] = [force_ratio]
    # df_value['Tz [Nm]'] = [Tz]
    df_value['mass flow rate inj [g/s]'] = [mfr_inj]
    # df_value['mass flow rate exit [g/s]'] = [mfr_exit]
    # df_value['Ispz [s]'] = [Ispz]
    # df_value['mean psi esit [deg]'] = [psi_exit_mean/2/np.pi*360]
    # df_value['mass flow rate error [%]'] = [mfr_error * 100]
    df_value['Vinj [m/s]'] = [Vinj]
    # df_value['fz_phai [N]'] = [fz_phai]
    # df_value['ftheta_phai [N]'] = [ftheta_phai]
    # df_value['force ratio_phai [-]'] = [force_ratio_phai]
    # df_value['Tz_phai [Nm]'] = [Tz_phai]
    # df_value['mass flow rate inj_phai [g/s]'] = [mfr_inj_phai]
    # # df_value['mass flow rate exit_phai [g/s]'] = [mfr_exit_phai]
    # df_value['Ispz_phai [s]'] = [Ispz_phai]
    # df_value['phai'] = [phai]
    # df_value['he_mean'] = [he_mean]
    # df_value['he_total_mean'] = [he_total_mean]

    df_value['mean_tau [Pa]'] = [mean_tau]
    df_value['mean_tau_psi_sl [Pa]'] = [mean_tau_psi_sl]

    df_value['mean_tau_N [Pa]'] = [mean_tau_N]
    df_value['mean_tau_psi_sl_N [Pa]'] = [mean_tau_psi_sl_N]

    df_value['Tz [Nm]'] = [Tz]
    df_value['Tz_psi_sl [Nm]'] = [Tz_psi_sl]

    df_value['Tz_N [Nm]'] = [Tz_N]
    df_value['Tz_psi_sl [Nm]'] = [Tz_psi_sl]

    df_value.to_csv(str(atm) + "_value.csv",index=False)
    os.chdir("../")
    os.chdir("../")

atm_max = 3.
atm_min = 0.1
delta = 0.1
total = (atm_max-atm_min)/delta
for atm_num in range(int(total)):
    atm = atm_min + (atm_max-atm_min)/ float(total)*atm_num
    func_atm(atm)

def func_elapsedtime(start):
    elapsed_time = time.time() - start  #処理時間算出
    print(elapsed_time)                 #処理時間表示
func_elapsedtime(start)






























































