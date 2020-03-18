%% initialize
clc
clear all
close all

[fontName,fontSize,fontSize_legend,bRGY,scrsz] = f_plotting;
p = f_physicalConstants;

%%
tau_a = 50e-9;
tau_b = 70e-9;

Ic_jj_vec = [40 60 80]*1e-6;
r_jj_vec = 250e-6./Ic_jj_vec;

r_spd_vec = tau_a*r_jj_vec./(tau_b-tau_a);
L_spd_vec = tau_b*r_spd_vec;

%%
for ii = 1:length(Ic_jj_vec)
    fprintf('Ic = %g uA:\n  r_spd = %g Ohm\n  L_spd = %g nH\n  tau_a = %g ns\n  tau_b = %g ns\n\n',Ic_jj_vec(ii)*1e6,r_spd_vec(ii),L_spd_vec(ii)*1e9,1e9*L_spd_vec(ii)/(r_spd_vec(ii)+r_jj_vec(ii)),1e9*L_spd_vec(ii)/r_spd_vec(ii))
end