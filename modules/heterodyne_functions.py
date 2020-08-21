from .config import c ,d_in,foc,freq,ω_0_in, ω_0_in_list
from .beam_waist_functions import cal_ω_0_out_and_d_out
import numpy as np
import scipy.optimize as opt # data fitting

# λ (wavelength) calc
def wavelen(f,c):
    return c/(f*1E9)
# rayleigh (rayleigh length) calc
def rayleigh(w_0,λ):
    return w_0**2*np.pi/λ
# w_z (beam radius) calc
def w_z(z, w_0, z_0):
    return w_0*np.sqrt(1+(z/z_0)**2)
# R (radius of wavefronts) calc
def R(z_0,z):
    return z*(1+(z_0/z)**2)
# η (guoy phase) calc
def η(z_0,z):
    return np.arctan(z/z_0)
# k (wave number)
def k(λ):
    return 2*np.pi/λ


# Frequency Array
#function for E (for z!=d_in+d_out)
def E(z, ρ, f_list=[95-6.35,95+6.35], A_0_list=[0.5,0.5], d_in=d_in, foc=foc, c=c,ω_0_in_list=ω_0_in_list):
    # Error message
    if len(f_list)!=len(A_0_list):
        raise ValueError("frequency list and amplitude list do not match")
    # save solutions here
    E_results=0
    # calculation
    for i,(f,A_0) in enumerate(zip(f_list,A_0_list)):
        ω_0_in=ω_0_in_list[min(ω_0_in_list.keys(), key=lambda x:abs(x-freq))]
        λ_=wavelen(f,c)
        d_out, w_0_ = cal_ω_0_out_and_d_out(d_in,ω_0_in,foc,λ_)
        # z transform
        z_=z-(d_in+d_out)
        # more params
        z_0_=rayleigh(w_0_,λ_)
        w_z_=w_z(z_, w_0_, z_0_)
        R_=R(z_0_,z_)
        η_=η(z_0_,z_)
        k_=k(λ_)
        E_results+= A_0 * np.exp(- (ρ/w_z_)**2) * np.exp(1j*k_*ρ**2/(2*R_)) * np.exp(1j*(k_*z_-η_))
#results
    return E_results
# function for i
def I_cross(z, ρ, f_list=[95-6.35,95+6.35], A_0_list=[0.5,0.5], d_in=d_in, ω_0_in=ω_0_in, foc=foc, c=c,ω_0_in_list=ω_0_in_list):
    peak = abs(E(z, 0, f_list, A_0_list))**2
    return abs(E(z, ρ, f_list, A_0_list))**2 -peak/np.e**2
# numerically obtain beam radius for combined frequencies
def I_comb(z, f_list=[95-6.35,95+6.35], A_0_list=[0.5,0.5], d_in=d_in, ω_0_in=ω_0_in, foc=foc, c=c,ω_0_in_list=ω_0_in_list):
    return abs(opt.fsolve(lambda x: I_cross(z=z, ρ=x,f_list=f_list, A_0_list=A_0_list), 30))[0]