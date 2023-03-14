from plt_function import plotfunc
import numpy as np
import matplotlib.pyplot as plt
wavelengtharr=np.linspace(700,800,101)
E0=1.67#ev
f=0.01
gamma=1*10**(-3)#ev
d1=40#nm
d2=0.65#nm
d3=15#nm
d4=300#nm
result=plotfunc("ag/al2o3",wavelengtharr,E0,f,gamma,d1,d2,d3,d4)
Rarr=result[1]
Rarrs=result[2]
contrastarr=result[3]
chiarr_real=result[4]
chiarr_imag=result[5]

fig,axis=plt.subplots(5)
axis[0].plot(wavelengtharr,Rarr,label="whole stack")
axis[0].legend()
axis[1].plot(wavelengtharr,Rarrs,label="substrate")
axis[1].legend()
axis[2].plot(wavelengtharr,contrastarr,label="contrast")
axis[2].legend()

axis[3].plot(wavelengtharr,chiarr_real,label="chi(real)")
axis[3].legend()

axis[4].plot(wavelengtharr,chiarr_imag,label="chi(imag)")
axis[4].legend()
axis[4].set_xlabel("wavelength(nm)")

plt.show()
