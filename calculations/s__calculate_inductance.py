#%% init
import numpy as np
from f__physical_constants import physical_constants

p = physical_constants()

#%% inputs
#see Van Duzer and Turner pg 110 - 115

#fundamental
lam = 90e-9#penetration depth of Nb
lam1 = lam#penetration depth of inductor
lam2 = lam#penetration depth of ground plane

#geometrical
b1 = 100e-9#thickness of microstrip
b2 = 100e-9#thickness of ground plane
d = 50e-9#thickness of dielectric spacer
K = 1.5#fringe factor

#%% total inductance per square
L = (p['mu0']*d/K)*(1+(lam1/d)*(np.tanh(b1/lam1))**(-1)+(lam2/d)*(np.tanh(b2/lam2))**(-1))
print('\nInductance per square = '+str(L*1e15)+' fH/sq\n')

