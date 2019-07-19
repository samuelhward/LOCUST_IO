import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import physical_constants as const
from scipy.special import erf

# Test particle parameters and energy grid
ma = const["deuteron mass"][0]
qa = const["elementary charge"][0]
Egrid = np.linspace(1, 120e3, 100)

# Plasma species parameters
Nion = 1
anum = np.array([2])
znum = np.array([1])
ne = 5.9567894e19
ni = np.array([5.9567894e19])
Te = 4.1185698e3
Ti = 9.4194805e3
va    = np.sqrt( const["elementary charge"][0]*2*Egrid/ma )

plt.figure()

h1 = plt.subplot(1, 1, 1)

h1.set_xlabel("Energy [keV]")
h1.set_ylabel("Energy drift [J/s]")

h1.set_xlim(0, 120)
h1.set_ylim(-3e-13,0)

# Bacgkround quantities, the last item is the Coulomb logarithm
eback = [Te, -const["elementary charge"][0], const["electron mass"][0],
         17]
iback = [Ti,  const["elementary charge"][0], const["proton mass"][0],
         21]
for back in [eback, iback]:
    mb = back[2];
    qb = back[1];

    clog = back[3];
    cab  = qa*qa*qb*qb * clog \
           /( 8*np.pi* np.power(const["electric constant"][0],2) );

    vb = np.sqrt(2*back[0]*const["elementary charge"][0]/mb);
    x  = va/vb;

    h = 2*erf(x)/x;
    g = ( ( 4*np.exp(-x*x) / ( np.sqrt(np.pi)*x ) - 2*erf(x) / (x*x) ) ) / x;

    nu = (ne*cab/(ma*ma)) * (1 + ma/mb) * (1/np.power(vb,3))*g;
    D  = (ne*cab/(ma*ma)) * (1/vb) * h;

    ve = 2*nu*(0.5*va*va*ma) + ma*D;


    h1.plot(Egrid/1e3, ve, marker='o', fillstyle='none');

plt.show()
