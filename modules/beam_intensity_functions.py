from .config import d_in,λ,ω_0_in, c, foc, ω_0_in_list
from .standard_canvas import d_out,ω_0_out, canvas
from .beam_waist_functions import cal_trans_prof, cal_ω_0_out_and_d_out, create_plot
from tqdm.notebook import tqdm
import numpy as np
# calculation of beam intensity

# calculation of beam intensity
def calculate_beam(r,ω_0,ω_z):
    return ((ω_0/ω_z)**2*np.exp(-2*r**2/ω_z**2))

# calculation plot
def calculate_one(λ=λ, d_out=d_out, ω_0_out=ω_0_out, d_in=d_in, ω_0_in=ω_0_in, canvas=canvas, i=None, leng=None):
    # calculations
    just = (canvas.shape[0] - 1) / 2  # justify
    # very naive intensity for post reflection beam (for establishing continuity)
    I_cell = (cal_trans_prof(d_in, ω_0_in, λ) / cal_trans_prof(-d_out, ω_0_out, λ) * ω_0_in / ω_0_out) ** 2
    image = np.copy(canvas)
    # progress bar: hide after completion?, labeling?
    if i == None:
        leave_bool = True
        desc_lab = "Progress"
    else:
        leave_bool = False
        desc_lab = "image {} of {}".format(i, leng)
    # the main calculations
    for z in tqdm(range(0, canvas.shape[1]), leave=leave_bool, desc=desc_lab):
        for r in range(0, canvas.shape[0]):
            if z <= d_in:
                ω_z = cal_trans_prof(z, ω_0_in, λ)
                image[r, z] = calculate_beam((r - just), ω_0_in, ω_z)
            else:
                ω_z = cal_trans_prof(z - (d_in + d_out), ω_0_out, λ)
                image[r, z] = calculate_beam((r - just), ω_0_out, ω_z) * I_cell
    return image

# save values for later use:
def calculate_many(frequencies, d_in=d_in, ω_0_in_list=ω_0_in_list, foc=foc, c=c, canvas=canvas):
    saved = []
    beam_saved=[]
    # calculate
    for i,ν in tqdm(enumerate(frequencies),total=len(frequencies),desc="total progress"):
        ω_in_now=ω_0_in_list[min(ω_0_in_list.keys(), key=lambda x:abs(x-ν))]
        λ_now=c/(ν*1E9)
        d_out,ω_0_out = cal_ω_0_out_and_d_out(d_in,ω_in_now,foc,c/(ν*1E9))
        saved.append(calculate_one(λ=λ_now, d_out=d_out, ω_0_out=ω_0_out, d_in=d_in, ω_0_in=ω_in_now, canvas=canvas, i=i+1, leng=len(frequencies)))
        z_vals=np.linspace(0,d_in+2*d_out,1000)
        r_vals=create_plot(X=z_vals, d_in=450, ω_0_in=ω_in_now, foc=foc, λ=λ_now, d_out=d_out, ω_0_out=ω_0_out)
        beam_saved.append([z_vals,r_vals])
    return saved, beam_saved
