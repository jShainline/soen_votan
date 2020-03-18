#%% initialization
import numpy as np
from pylab import *

# Import instruments
from amcc.instruments.srs_sim970 import SIM970
from amcc.instruments.srs_sim928 import SIM928
from amcc.instruments import Switchino

# Setup instruments
dmm = SIM970('GPIB0::4', sim900port = 7)
vs1 = SIM928('GPIB0::4', sim900port = 3)
vs2 = SIM928('GPIB0::4', sim900port = 4)
vs3 = SIM928('GPIB0::4', sim900port = 5)
switch = Switchino('COM7')

#import functions
from vt__meas_util import iv_sweep_srs__current_bias, flux_purge_srs, sq_voltage_vs_incoil_current


#%% measurement specifics

#squid ports
incoil_v_source_srs = vs1
squid_v_source_srs = vs2
addflux_v_srs = vs3
squid_v_meas_srs = dmm
squid_v_meas_srs_dmm_channel = 4

dmm.set_impedance(gigaohm=False, channel = squid_v_meas_srs_dmm_channel)

#resistors in series with voltage sources
incoil_i_source_res = 1e4
squid_i_source_res = 1e4
addflux_i_source_res = 1e4

incoil_v_source_srs.set_voltage(0)
squid_v_source_srs.set_voltage(0)
addflux_v_srs.set_voltage(0)

#%% purge flux

flux_purge_srs(squid_v_source_srs,squid_i_source_res,1e-3,15)

#%% get squid I-V in the absence of incoil flux bias
device_name = 'vt01_nw_22_sqb_device4'
sq_current_bias_values = np.arange(0,300e-6,1e-6)
V = iv_sweep_srs__current_bias(squid_v_source_srs,squid_v_meas_srs,squid_v_meas_srs_dmm_channel,sq_current_bias_values,squid_i_source_res, delay = 0.75, device_name = device_name)
    
I = sq_current_bias_values
fig, axes = plt.subplots(1,1)
axes.plot(V*1e3,I*1e6,label = '1')
axes.set_xlabel(r'Voltage across SQUID (mV)', fontsize=20)
axes.set_ylabel(r'Current applied to SQUID (uA)', fontsize=20)

#ylim((ymin_plot,ymax_plot))
#xlim((xmin_plot,xmax_plot))

#axes.legend(loc='best')
grid(True,which='both')
plt.show()

#p = np.polyfit(I,V,1)
#print('%0.1f Ohm' % (p[0]))

#title(str(np.around(p[0],decimals = 2))+' ohm')
#
#time.sleep(1)

#%% measure SQUID voltage as a function of incoil current for several values of SQUID current bias
incoil_current_bias_values = np.arange(0,1e-3,1e-6)
sq_Ic = 180e-6
#sq_current_bias_values = np.linspace(0.5,1.5,10)*sq_Ic 
sq_current_bias_values = np.linspace(170e-6,210e-6,5) 
measurement_delay = 0.75

I_sq_bias,I_incoil_current,V_sq_meas = sq_voltage_vs_incoil_current(squid_v_source_srs,squid_i_source_res,squid_v_meas_srs,squid_v_meas_srs_dmm_channel,
                                                         incoil_v_source_srs,incoil_i_source_res,sq_current_bias_values,incoil_current_bias_values,measurement_delay,device_name)

    
fig, axes = plt.subplots(1,1)
for ii in range(len(I_sq_bias)):
    v_vec = V_sq_meas[:,ii]*1e3
    axes.plot(I_incoil_current[:]*1e6, v_vec, 'o-', linewidth = 1, markersize = 3, label = 'I_sq = {0} uA'.format(I_sq_bias[ii]*1e6))
    axes.set_ylabel(r'Voltage across SQUID (mV)', fontsize=20)
    axes.set_xlabel(r'Incoil current (uA)', fontsize=20)

#ylim((ymin_plot,ymax_plot))
#xlim((xmin_plot,xmax_plot))

axes.legend(loc='best')
grid(True,which='both')
plt.show()

time.sleep(1)
