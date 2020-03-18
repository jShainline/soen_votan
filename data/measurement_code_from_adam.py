#%%

# Import instruments
from amcc.instruments.srs_sim970 import SIM970
from amcc.instruments.srs_sim928 import SIM928
from amcc.instruments import Switchino

# Setup instruments
dmm = SIM970('GPIB0::4', sim900port = 7)
vs = SIM928('GPIB0::4', sim900port = 4)
switch = Switchino('COM7')

dmm.set_impedance(gigaohm=False, channel = 2)
dmm.set_impedance(gigaohm=False, channel = 3)

#%%
def run_iv_sweep_srs(voltages, R_series, delay = 0.75):
    vs.reset()
    vs.set_output(True)
    time.sleep(2)
    V = []
    I = []
    for v in voltages:
        vs.set_voltage(v)
        time.sleep(delay)
#        v1 = dmm.read_voltage(channel = 1)
        v1 = v
        v2 = dmm.read_voltage(channel = 2)
        v3 = dmm.read_voltage(channel = 3)
        V.append(v3)
        I.append((v1-v2)/R_series)
    vs.set_voltage(0)
    return np.array(V),np.array(I)


def iv(port_pair = [1,2], voltages = np.linspace(0,1,25),
             R_series = 10e3, delay = 0.75):
    switch.select_ports(port_pair = port_pair)
    V, I = run_iv_sweep_srs(voltages, R_series, delay = delay)
    return V, I


port_pairs = [
        [1,2],
        [5,6],
        ]
for port_pair in port_pairs:
    voltages = np.linspace(0,1,15)
    V, I = iv(port_pair = port_pair, voltages = voltages, R_series = 10e3, delay = 1.5)
    
    figure()
    plot(V*1e3,I*1e6,'.')
    xlabel('Voltage (mV)')
    ylabel('Current (uA)')
    title(str(port_pair))

    p = np.polyfit(I,V,1)
    print(port_pair)
    print('%0.1f Ohm' % (p[0]))
