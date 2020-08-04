from .config import d_in,foc,λ,ω_0_in
import numpy as np


# function for calc. transversile profile
def cal_trans_prof(z,ω_0,λ): ### λ aprrox by λ=c/f
    z_r=np.pi*ω_0**2/λ
    return    ω_0*np.sqrt(1+(z/z_r)**2)


# function for calc. d_out and w_0_out
def cal_ω_0_out_and_d_out(d_in=d_in,ω_0_in=ω_0_in,foc=foc,λ=λ):
    C=-1/foc
    z_c=np.pi*ω_0_in**2/λ
    d_out=-(d_in*(C*d_in+1)+C*z_c**2)/((C*d_in+1)**2+C**2*z_c**2)
    ω_0_out=ω_0_in/np.sqrt((C*d_in+1)**2+C**2*z_c**2)
    return d_out,ω_0_out


# function that creates beam waist plot
def create_plot(X, d_in=d_in, ω_0_in=ω_0_in, foc=foc, λ=λ, d_out=None, ω_0_out=None, ):
    if (("d_out" and "ω_0_out" in locals()) and (d_out == None or ω_0_out == None)) == True:  # measure to save compute
        d_out, ω_0_out = cal_ω_0_out_and_d_out(d_in, ω_0_in, foc, λ)

    if (np.issubdtype(type(X), int) == True or np.issubdtype(type(X), float) == True) == False:
        res = np.copy(X)
        for i, x in enumerate(res):
            if x <= d_in:
                res[i] = cal_trans_prof(x, ω_0_in, λ)
            else:
                res[i] = cal_trans_prof((d_out + d_in) - x, ω_0_out, λ)
    else:
        if X <= d_in:
            res = cal_trans_prof(X, ω_0_in, λ)
        else:
            res = cal_trans_prof((d_out + d_in) - X, ω_0_out, λ)

    return res

