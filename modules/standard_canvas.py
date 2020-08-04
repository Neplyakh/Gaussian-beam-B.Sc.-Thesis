# Beam intensity plotting

from .beam_waist_functions import  cal_ω_0_out_and_d_out
from .config import d_in
import numpy as np

# calculate d_out and ω_out for canvas
d_out,ω_0_out = cal_ω_0_out_and_d_out()


# creating canvas
max_r=100 #in mm
max_z=np.int((np.round(d_in+1.5*d_out))) # in mm
canvas_size=(max_r*2+1,max_z) #uneven for symmetry
canvas=np.zeros(canvas_size)