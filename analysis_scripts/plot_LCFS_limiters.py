import matplotlib.pyplot as plt
import numpy as np
pi=np.pi

for filename,colour in zip(['LCFS_DIII-D.dat_1.05','LCFS_DIII-D.dat_1.10','LCFS_DIII-D.dat_1.15','LCFS_DIII-D.dat_1.20','LCFS_DIII-D.dat_1.30','LCFS_DIII-D.dat_1.40','LCFS_DIII-D.dat_1.50',],['r','g','b','y','m','k','r',]):
    file=open(filename)
    lines=file.readlines()
    del(lines[0])
    del(lines[0])
    r=[]
    for line in lines:
        r.append(float(line.split()[0]))
    #angles_new=np.linspace(0.,(3599./3600.)*2.*pi,3600)
    angles_new=np.linspace(pi,(3599./3600.)*3.*pi,3600)
    r=np.array(r)
    plt.polar(angles_new,np.sqrt(r),colour)

plt.show()