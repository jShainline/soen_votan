%% initialize
clc; close all;
% clear all

[fontName,fontSize,fontSize_legend,bRGY,scrsz] = f_plotting;
fontSize = 16;
p = f_physicalConstants;


%% Input data

n_res__au = 16;
n_sq__au = [1 2 4 8 16];
res_au__die_1 = [15.75 17.25 19.51 27.44 47.08];

n_res__pdau = 16;
n_sq__pdau = [2 4 8];
res_pdau__die_1 = [60.340 121.639 240.67];

%% fits

p_au = polyfit(n_res__au*n_sq__au,res_au__die_1,1);
n_sq__au__dense = linspace(0,n_res__au*n_sq__au(end),100);
res_au = polyval(p_au,n_sq__au__dense);
res_per_sq__au = p_au(1);
res_cntct__au = p_au(2)/(2*n_res__au);

p_pdau = polyfit(n_res__pdau*n_sq__pdau,res_pdau__die_1,1);
n_sq__pdau__dense = linspace(0,n_res__pdau*n_sq__pdau(end),100);
res_pdau = polyval(p_pdau,n_sq__pdau__dense);
res_per_sq__pdau = p_pdau(1);
res_cntct__pdau = p_pdau(2)/(2*n_res__pdau);


%% plots

%Au
figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
plot(n_sq__au__dense,res_au,'Color',bRGY(8,:),'LineStyle','-','LineWidth',2)
hold on
plot(n_res__au*n_sq__au,res_au__die_1,'Color',bRGY(3,:),'LineStyle','-','LineWidth',2,'Marker','s','MarkerSize',4,'MarkerFaceColor',bRGY(1,:),'MarkerEdgeColor',bRGY(5,:))
legend('Fit','Data')
xlabel('Num Squares','FontSize',fontSize,'FontName','Times')
ylabel('Resistance [\Omega]','FontSize',fontSize,'FontName','Times')
set(gca,'FontSize',fontSize,'FontName',fontName)
grid on
title(sprintf('Au resistance\nResistance per square = %g Ohms; Contact resistance = %g Ohms',res_per_sq__au,res_cntct__au),'FontSize',fontSize,'FontName','Times')

%PdAu
figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
plot(n_sq__pdau__dense,res_pdau,'Color',bRGY(8,:),'LineStyle','-','LineWidth',2)
hold on
plot(n_res__pdau*n_sq__pdau,res_pdau__die_1,'Color',bRGY(3,:),'LineStyle','-','LineWidth',2,'Marker','s','MarkerSize',4,'MarkerFaceColor',bRGY(1,:),'MarkerEdgeColor',bRGY(5,:))
legend('Fit','Data')
xlabel('Num Squares','FontSize',fontSize,'FontName','Times')
ylabel('Resistance [\Omega]','FontSize',fontSize,'FontName','Times')
set(gca,'FontSize',fontSize,'FontName',fontName)
grid on
title(sprintf('PdAu resistance\nResistance per square = %g Ohms; Contact resistance = %g Ohms',res_per_sq__pdau,res_cntct__pdau),'FontSize',fontSize,'FontName','Times')
