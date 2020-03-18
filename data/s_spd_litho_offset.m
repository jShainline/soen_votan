%% initialize
clc; 
close all;
clear all

[fontName,fontSize,fontSize_legend,bRGY,scrsz] = f_plotting;
fontSize = 16;
p = f_physicalConstants;

%% Input data
spd_mask_width = 500:100:1300;
spd_actual_width = [278 444.3 565.1 658.6 765.7 872.8 966.2 1078 1180];

%% fit
p_fit = polyfit(spd_mask_width,spd_actual_width,1);
width_dense = linspace(0,spd_mask_width(end),100);
spd_actual_width_dense = polyval(p_fit,width_dense);
slope = p_fit(1);
offset = p_fit(2);

%% plots

%Au
figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
plot(width_dense,spd_actual_width_dense,'Color',bRGY(8,:),'LineStyle','-','LineWidth',2)
hold on
plot(spd_mask_width,spd_actual_width,'Color',bRGY(3,:),'LineStyle','-','LineWidth',2,'Marker','s','MarkerSize',4,'MarkerFaceColor',bRGY(1,:),'MarkerEdgeColor',bRGY(5,:))
legend('Fit','Measured')
xlabel('Width on mask [nm]','FontSize',fontSize,'FontName','Times')
ylabel('Width on chip [nm]','FontSize',fontSize,'FontName','Times')
set(gca,'FontSize',fontSize,'FontName',fontName)
xlim([0 1400])
ylim([0 1300])
grid on
title(sprintf('SPD litho offset\nactual\\_width = %g * mask\\_width %g\nmask\\_width = (desired\\_width+%g)/%g',slope,offset,-offset,slope),'FontSize',fontSize,'FontName','Times')
