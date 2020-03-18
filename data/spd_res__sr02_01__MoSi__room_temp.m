%% initialize
clc; close all;
clear all

[fontName,fontSize,fontSize_legend,bRGY,scrsz] = f_plotting;
fontSize = 16;
p = f_physicalConstants;


%% Input data

w_wire_vec = 0.5:0.1:1.3;%width of wire on mask
num_squares_vec__meanders = linspace(3500,10000,9);
num_squares_vec__stitches = [0.5 1 2 4 8 16];

res__varying_w_wire = [3.455 2.545 2.347 2.254 2.169 2.116 2.097 2.073 2.046];
res__varying_num_squares = [1.079 1.307 1.545 1.797 2.051 2.316 2.527 2.541 3.146];
res__consistency_test = [2.179 2.186 2.159 2.146 2.172 2.124 2.157 2.165 2.185];
res__stitches = [49.1e-3 4.791 9.52 18.85 37.46 74.4];


%% fits

p_w_wire = polyfit(w_wire_vec,res__varying_w_wire,1);
w_wire_vec__dense = linspace(w_wire_vec(1),w_wire_vec(end),100);
res_vs_width__fit = polyval(p_w_wire,w_wire_vec__dense);
res_per_sq__pdau__from_width = p_w_wire(1)/7000;
% res_cntct__pdau = p_pdau(2)/(2*n_res__pdau);


%% plots

%resistance vs width
figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
plot(w_wire_vec__dense,res_vs_width__fit,'Color',bRGY(8,:),'LineStyle','-','LineWidth',2)
hold on
plot(w_wire_vec,res__varying_w_wire,'Color',bRGY(3,:),'LineStyle','-','LineWidth',2,'Marker','s','MarkerSize',4,'MarkerFaceColor',bRGY(1,:),'MarkerEdgeColor',bRGY(5,:))
legend('Linear fit','Data')
xlabel('Wire width [\mu m]','FontSize',fontSize,'FontName','Times')
ylabel('Resistance [M\Omega]','FontSize',fontSize,'FontName','Times')
set(gca,'FontSize',fontSize,'FontName',fontName)
grid on
title(sprintf('MoSi resistance vs wire width\nWafer:sr02\\_01, num squares = 7000',res_per_sq__pdau__from_width),'FontSize',fontSize,'FontName','Times')

