%% initialize
clc
clear all
close all

[fontName,fontSize,fontSize_legend,bRGY,scrsz] = f_plotting;
p = f_physicalConstants;

%% time-averaged jj voltage

I_c = 200e-6;
r = 0.5;
I_vec = linspace(I_c,3*I_c,100);
V_time_averaged = r*sqrt(I_vec.^2-I_c^2);

%%
figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
plot(I_vec*1e6,V_time_averaged*1e6,'Color',bRGY(3,:),'LineStyle','-','LineWidth',3)
ylabel('Voltage [uV]','FontSize',fontSize,'FontName','Times')
xlabel('Current [uA]','FontSize',fontSize,'FontName','Times')
set(gca,'FontSize',fontSize,'FontName',fontName)
% ylim([0 1.1])
