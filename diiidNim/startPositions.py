import numpy as np



r_m = 1.750
z_m = 2.537e-2

a=0.567
kappa=1.829


ppc = 25
cords = 4

npts = ppc * cords


phi = 0.0
f1=open("start_positions.dat",'w',)
print npts
f1.write(str(npts)+'\n')
for it in np.linspace(0,2*np.pi,cords,endpoint=False):
    for ir in np.linspace(1.0,0,ppc,endpoint=False):
       thisZ = z_m + ir * kappa * a * np.sin(it)
       thisR = r_m + ir * a * np.cos(it)
       print thisR,'\b,', thisZ, '\b,', phi
       f1.write( str(thisR)+', '+ str(thisZ)+ ', '+ str(phi)+'\n')
    
f1.close()    