%% initialize
clc; close all;
% clear all

[fontName,fontSize,fontSize_legend,bRGY,scrsz] = f_plotting;
fontSize = 16;
p = f_physicalConstants;


%% Input data
n1_au = 2500;
res_au__die_1 = [560.8 1105 2199 4471];
res_au__die_6 = [498 953 1893 3959];

n1_pdau = 250;
res_pdau__die_1 = [582.4 1141 2267 4602];
res_pdau__die_6 = [516 980 1941 4047];

%% fits
num_squares_vec__au = n1_au*[1 2 4 8];
p_au = polyfit(num_squares_vec__au,res_au__die_1,1);
num_squares_au__dense = linspace(0,num_squares_vec__au(end),100);
res_au = polyval(p_au,num_squares_au__dense);
res_per_sq__au = p_au(1);
res_cntct__au = p_au(2);

p_au = polyfit(num_squares_vec__au,res_au__die_6,1);
num_squares_au__dense = linspace(0,num_squares_vec__au(end),100);
res_per_sq__au_edge = p_au(1);
res_cntct__au_edge = p_au(2);

num_squares_vec__pdau = n1_pdau*[1 2 4 8];
p_pdau = polyfit(num_squares_vec__pdau,res_pdau__die_1,1);
num_squares_pdau__dense = linspace(0,num_squares_vec__pdau(end),100);
res_pdau = polyval(p_pdau,num_squares_pdau__dense);
res_per_sq__pdau = p_pdau(1);
res_cntct__pdau = p_pdau(2);

p_pdau = polyfit(num_squares_vec__pdau,res_pdau__die_6,1);
res_per_sq__pdau_edge = p_pdau(1);
res_cntct__pdau_edge = p_pdau(2);

%% plots

%Au
figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
plot(num_squares_au__dense,res_au,'Color',bRGY(8,:),'LineStyle','-','LineWidth',2)
hold on
plot(num_squares_vec__au,res_au__die_1,'Color',bRGY(3,:),'LineStyle','-','LineWidth',2,'Marker','s','MarkerSize',4,'MarkerFaceColor',bRGY(1,:),'MarkerEdgeColor',bRGY(5,:))
plot(num_squares_vec__au,res_au__die_6,'Color',bRGY(13,:),'LineStyle','-','LineWidth',2,'Marker','v','MarkerSize',4,'MarkerFaceColor',bRGY(11,:),'MarkerEdgeColor',bRGY(15,:))
legend('Fit','Center','Edge')
xlabel('Num Squares','FontSize',fontSize,'FontName','Times')
ylabel('Resistance [\Omega]','FontSize',fontSize,'FontName','Times')
set(gca,'FontSize',fontSize,'FontName',fontName)
grid on
title(sprintf('Au resistance, room temp\nCenter of wafer: resistance per square = %g Ohms; contact resistance = %g Ohms\nEdge of wafer: resistance per square = %g Ohms; contact resistance = %g Ohms\nCross wafer variation: %g%%',res_per_sq__au,res_cntct__au,res_per_sq__au_edge,res_cntct__au_edge,2*(res_per_sq__au-res_per_sq__au_edge)/(res_per_sq__au+res_per_sq__au_edge)*100),'FontSize',fontSize,'FontName','Times')

%PdAu
figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
plot(num_squares_pdau__dense,res_pdau,'Color',bRGY(8,:),'LineStyle','-','LineWidth',2)
hold on
plot(num_squares_vec__pdau,res_pdau__die_1,'Color',bRGY(3,:),'LineStyle','-','LineWidth',2,'Marker','s','MarkerSize',4,'MarkerFaceColor',bRGY(1,:),'MarkerEdgeColor',bRGY(5,:))
plot(num_squares_vec__pdau,res_pdau__die_6,'Color',bRGY(13,:),'LineStyle','-','LineWidth',2,'Marker','v','MarkerSize',4,'MarkerFaceColor',bRGY(11,:),'MarkerEdgeColor',bRGY(15,:))
legend('Fit','Center','Edge')
xlabel('Num Squares','FontSize',fontSize,'FontName','Times')
ylabel('Resistance [\Omega]','FontSize',fontSize,'FontName','Times')
set(gca,'FontSize',fontSize,'FontName',fontName)
grid on
title(sprintf('PdAu resistance\nCenter of wafer: resistance per square = %g Ohms; contact resistance = %g Ohms\nEdge of wafer: resistance per square = %g Ohms; contact resistance = %g Ohms\nCross wafer variation: %g%%',res_per_sq__pdau,res_cntct__pdau,res_per_sq__pdau_edge,res_cntct__pdau_edge,2*(res_per_sq__pdau-res_per_sq__pdau_edge)/(res_per_sq__pdau+res_per_sq__pdau_edge)*100),'FontSize',fontSize,'FontName','Times')
