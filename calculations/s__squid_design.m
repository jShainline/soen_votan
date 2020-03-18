%% initialize
clc
clear all
% close all

[fontName,fontSize,fontSize_legend,bRGY,scrsz] = f_plotting;
p = f_physicalConstants;

%% inputs

Ic_uA = 100;%junction critical current in uA
I_applied_max = 20e-6;%maximum applied current
Jc_kAcm2 = 1;%critical current density in kA/cm^2
beta_c = 0.3;
beta_L = 1;
k = 0.7;%mutual inductor coupling factor

%% junctions

C_per_current = 1.15e-14/1e-6;%junction capacitance per amp of critical current
Jc = Jc_kAcm2*1e3*1e4;%critical current density in A/m^2
Ic = Ic_uA*1e-6;
Aj = Ic/Jc;
dj = (4*Aj/pi)^(1/2);
Cj = C_per_current*Ic;
R_shunt = ( (p.Phi0*beta_c)/(2*pi*Ic*Cj) )^(1/2);

fprintf('JJs:\n')
fprintf('Jc: %g kA/cm^2\n',Jc*1e-7)
fprintf('Ic: %g uA\n',Ic*1e6)
fprintf('Aj: %g um^2\n',Aj*1e12)
fprintf('dj: %g um\n',dj*1e6)
fprintf('Cj: %g fF\n',Cj*1e15)
fprintf('R_shunt: %g ohm\n',R_shunt)

%% inductors

L1 = p.Phi0*Ic/(2*k^2*beta_L*(I_applied_max)^2);
L2 = p.Phi0*beta_L/(2*Ic);
V_phi = Ic*R_shunt/2;

fprintf('\nSQUIDs:\n')
fprintf('Mutual inductor coupling k: %g\n',k)
fprintf('SQUID input inductance L1: %g pH\n',L1*1e12)
fprintf('SQUID self-inductance L2: %g pH\n',L2*1e12)
fprintf('V_phi: %g uV\n',V_phi*1e6)

%% washer 


