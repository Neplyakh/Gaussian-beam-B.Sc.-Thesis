# Define Physical constants
c=2.99792458E11 # in mm/s

# Define general input paramters
ω_0_in_list = {70:7.85,75:7.67,80:7.47,85:7.26,90:7.03,95:7.16,100:7.23,105:7.26,110:7.24}
d_in=449.7 # in mm
freq=95
ω_0_in=ω_0_in_list[min(ω_0_in_list.keys(), key=lambda x:abs(x-freq))]
freq=freq*1E9
# elliptical mirror
R1=1
R2=1
# extracted params==>
foc=339.5
λ=c/freq


