import numpy as np

def physical_constants():

    p = dict(h = 6.62606957e-34,#Planck's constant in kg m^2/s
         hBar = 6.62606957e-34/2/np.pi,
         c = 299792458,#speed of light in meters per second
         epsilon0 = 8.854187817e-12,#permittivity of free space in farads per meter
         mu0 = 4*np.pi*1e-7,#permeability of free space in volt seconds per amp meter
         kB = 1.3806e-23,#Boltzmann's constant
         eE = 1.60217657e-19,#electron charge in coulombs
         mE = 9.10938291e-31,#mass of electron in kg
         eV = 1.60217657e-19,#joules per eV
         Ry = 9.10938291e-31*1.60217657e-19**4/(8*8.854187817e-12**2*(6.62606957e-34/2/np.pi)**3*299792458),#13.3*eV;#Rydberg in joules
         a0 = 4*np.pi*8.854187817e-12*(6.62606957e-34/2/np.pi)**2/(9.10938291e-31*1.60217657e-19**2),#estimate of Bohr radius
         Phi0 = 6.62606957e-34/(2*1.60217657e-19)#flux quantum
         )

    return p