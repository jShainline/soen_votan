%% initialize
clc
% clear all
close all

[fontName,fontSize,fontSize_legend,bRGY,scrsz] = f_plotting;
p = f_physicalConstants;

%%
contact_resistance = 30;
jc_vec = [2e7 1.5e7 1e7 0.5e7];
d_vec_jja = sqrt(4*40e-6./(pi*jc_vec));
d_vec_jja = linspace(d_vec_jja(1),d_vec_jja(end),8);
d_vec_jjb = sqrt(4*100e-6./(pi*jc_vec));
d_vec_jjb = linspace(d_vec_jjb(1),d_vec_jjb(end),8);
areas_mask_jja = 1e12*pi*(d_vec_jja/2).^2; 
areas_mask_jjb = 1e12*pi*(d_vec_jjb/2).^2; 

%40 jjs in series
res_jjb_40 = [514 426 362 321 294 273 257 244];
res_jjb_40 = res_jjb_40-contact_resistance;

%80 jjs in series
res_jjb_80 = [988 808 690 612 552 511 475 455];
res_jjb_80 = res_jjb_80-contact_resistance;

%160 jjs in series
res_jjb_160 = [1959 1583 1337 1183 1064 974 913 870];
res_jjb_160 = res_jjb_160-contact_resistance;

%% converting to Ic
Delta_Nb = 1.278e-22;
prefactor = pi*Delta_Nb/(2*p.eE);

Ic_jjb_40 = prefactor./(res_jjb_40/40);
Ic_jjb_80 = prefactor./(res_jjb_80/80);
Ic_jjb_160 = prefactor./(res_jjb_160/160);

%%

figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
plot(areas_mask_jjb,res_jjb_40,'Color',bRGY(3,:),'LineStyle','-','LineWidth',3)
hold on
plot(areas_mask_jjb,res_jjb_80,'Color',bRGY(8,:),'LineStyle','-','LineWidth',3)
plot(areas_mask_jjb,res_jjb_160,'Color',bRGY(13,:),'LineStyle','-','LineWidth',3)
legend('40 jjs','80 jjs','160 jjs')
ylabel('Room Temp. Resistance [\Omega]','FontSize',fontSize,'FontName','Times')
xlabel('Junction Area [\mu m^2]','FontSize',fontSize,'FontName','Times')
set(gca,'FontSize',fontSize,'FontName',fontName)
% k1 = gtext(figureCaptions(1:length(figureCaptions)));
% set(k1,'FontSize',fontSize_legend,'FontName','Times')

figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
plot(areas_mask_jjb,res_jjb_40/40,'Color',bRGY(3,:),'LineStyle','-','LineWidth',3)
hold on
plot(areas_mask_jjb,res_jjb_80/80,'Color',bRGY(8,:),'LineStyle','-','LineWidth',3)
plot(areas_mask_jjb,res_jjb_160/160,'Color',bRGY(13,:),'LineStyle','-','LineWidth',3)
legend('40 jjs','80 jjs','160 jjs')
ylabel('Room Temp. Resistance per JJ [\Omega]','FontSize',fontSize,'FontName','Times')
xlabel('Junction Area [\mu m^2]','FontSize',fontSize,'FontName','Times')
set(gca,'FontSize',fontSize,'FontName',fontName)
% k1 = gtext(figureCaptions(1:length(figureCaptions)));
% set(k1,'FontSize',fontSize_legend,'FontName','Times')

figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
plot(areas_mask_jjb,1e6*Ic_jjb_40,'Color',bRGY(3,:),'LineStyle','-','LineWidth',3)
hold on
plot(areas_mask_jjb,1e6*Ic_jjb_80,'Color',bRGY(8,:),'LineStyle','-','LineWidth',3)
plot(areas_mask_jjb,1e6*Ic_jjb_160,'Color',bRGY(13,:),'LineStyle','-','LineWidth',3)
legend('40 jjs','80 jjs','160 jjs')
ylabel('I_c [\mu A]','FontSize',fontSize,'FontName','Times')
xlabel('Junction Area [\mu m^2]','FontSize',fontSize,'FontName','Times')
set(gca,'FontSize',fontSize,'FontName',fontName)
% k1 = gtext(figureCaptions(1:length(figureCaptions)));
% set(k1,'FontSize',fontSize_legend,'FontName','Times')