#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 11 13:31:02 2019

@author: bryceprimavera
"""

def washer_no_plane(w_wide = 25, d_in = 13, lambda1 = .039):
    file = open('squid_no_plane.inp','w')  
    file.write('*Washer\n') 
    file.write('.units um\n')
    file.write('.default z=0 lambda = %f\n'%(lambda1)) # Penetration Depth
    
    thickness = 0.2
    d_in = 13
    w_wide = 25
    w_slit = 2
    
    # Actual Washer
    
    #Defines the ground plane with 3 corners, fast henry figures out the 4th, "+" means continue the line
    file.write('g1 x1=%f y1=%f z1=%f\n'%(0, 0, 0))
    file.write('+ x2=%f y2=%f z2=%f\n'%(0, d_in + 2*w_wide, 0))
    file.write('+ x3=%f y3=%f z3=%f\n'%(d_in + 2*w_wide, d_in + 2*w_wide, 0))
    
    file.write('+ thick=%f\n'%(thickness))
    
    #how many segments to divide x and y directions into
    file.write('+ seg1=%i seg2=%i\n'%(100, 100))
    
#    file.write('+ nin1 (25,10)\n')
#    file.write('+ nout1 (35,10 )\n')
    
    # Punch out center hole
    file.write('+ hole rect (%f,%f,%f,%f,%f,%f)\n'%(w_wide, w_wide, 0, d_in + w_wide, d_in + w_wide, 0))
    
    #Punch out gap
    file.write('+ hole rect (%f,%f,%f,%f,%f,%f)\n'%(.5*(d_in + 2*w_wide-w_slit), 0, 0, .5*(d_in + 2*w_wide+w_slit), w_wide, 0))
    
    #Define input and output nodes - fast henry has to find the nodes in the plane that are closest to these values
    file.write('+ nin1 (%f,%f,%f)\n'%(.5*(d_in + 2*w_wide-1.5*w_slit) - .1, .5*w_wide, 0))
    file.write('+ nout1 (%f,%f,%f)\n'%(.5*(d_in + 2*w_wide+1.5*w_slit) + .1, .5*w_wide, 0))
    #file.write('+ nplane (%f,%f,%f)\n'%(.5*(d_in + 2*w_wide+1.5*w_slit) + .1, .5*w_wide, -2*thickness))

    # I don't think this is necessary, but it's the way the examples do it
    file.write('.equiv ninput nin1\n')
    file.write('.equiv noutput nout1\n')
    #file.write('.equiv nplane nout1\n')

    # Says which nodes to use as ports in the impedance matrix
    file.write('.external ninput noutput\n')
    
    #Frequencies to test out - doesn't really matter
    file.write('.freq fmin=.1 fmax = 10 ndec = 1\n')

    file.write('.end')
    file.close()
    
    return


def washer_plane(w_wide = 25, d_in = 13, lambda1 = .039):
    file = open('squid_plane.inp','w')  
    file.write('*Washer\n') 
    file.write('.units um\n')
    file.write('.default z=0 lambda = %f\n'%(lambda1))
    
    thickness = 0.2
    d_in = 13
    w_wide = 25
    w_slit = 2
    
    # Big Ground Plane Below Squid
    file.write('g0 x1=%f y1=%f z1=%f\n'%(-(d_in + 2*w_wide), -(d_in + 2*w_wide), -2*thickness))
    file.write('+ x2=%f y2=%f z2=%f\n'%(-(d_in + 2*w_wide), 2*(d_in + 2*w_wide), -2*thickness))
    file.write('+ x3=%f y3=%f z3=%f\n'%(2*(d_in + 2*w_wide), 2*(d_in + 2*w_wide), -2*thickness))
    file.write('+ thick=%f\n'%(thickness))
    file.write('+ seg1=%i seg2=%i\n'%(49, 49))
    
    # Actual Washer
    file.write('g1 x1=%f y1=%f z1=%f\n'%(0, 0, 0))
    file.write('+ x2=%f y2=%f z2=%f\n'%(0, d_in + 2*w_wide, 0))
    file.write('+ x3=%f y3=%f z3=%f\n'%(d_in + 2*w_wide, d_in + 2*w_wide, 0))
    file.write('+ thick=%f\n'%(thickness))
    file.write('+ seg1=%i seg2=%i\n'%(100, 100))
    
#    file.write('+ nin1 (25,10)\n')
#    file.write('+ nout1 (35,10 )\n')
    
    file.write('+ hole rect (%f,%f,%f,%f,%f,%f)\n'%(w_wide, w_wide, 0, d_in + w_wide, d_in + w_wide, 0))
    file.write('+ hole rect (%f,%f,%f,%f,%f,%f)\n'%(.5*(d_in + 2*w_wide-w_slit), 0, 0, .5*(d_in + 2*w_wide+w_slit), w_wide, 0))
        
    file.write('+ nin1 (%f,%f,%f)\n'%(.5*(d_in + 2*w_wide-1.5*w_slit) - .1, .5*w_wide, 0))
    file.write('+ nout1 (%f,%f,%f)\n'%(.5*(d_in + 2*w_wide+1.5*w_slit) + .1, .5*w_wide, 0))
    #file.write('+ nplane (%f,%f,%f)\n'%(.5*(d_in + 2*w_wide+1.5*w_slit) + .1, .5*w_wide, -2*thickness))

    file.write('.equiv ninput nin1\n')
    file.write('.equiv noutput nout1\n')
    #file.write('.equiv nplane nout1\n')

    file.write('.external ninput noutput\n')
    file.write('.freq fmin=.1 fmax = 10 ndec = 1\n')

    file.write('.end')
    file.close()
    
    return

    
