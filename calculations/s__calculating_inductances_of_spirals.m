%% initialize
clc
clear all
% close all

[fontName,fontSize,fontSize_legend,bRGY,scrsz] = f_plotting;
p = f_physicalConstants;

%%
geometry = 'square';

switch geometry
    case 'square'
        k1 = 2.34;
        k2 = 2.75;
    case 'octagon'
        k1 = 2.25;
        k2 = 3.55;
end

d_out = 20e-6;
d_in = 3e-6;

d_avg = (d_out+d_in)/2;
rho = (d_out-d_in)/(d_out+d_in);

w = 1e-6;
s = 1e-6;

per = w+s;
n = floor((d_out-d_in)/2/per);%10;

L_mw = k1*p.mu0*n^2*d_avg/(1+k2*rho);

fprintf([geometry ':\n'])
fprintf('w = %g um\n',w*1e6)
fprintf('s = %g um\n',s*1e6)
fprintf('d_out = %g um\n',d_out*1e6)
fprintf('d_in = %g um\n',d_in*1e6)
fprintf('n = %g \n',n)
fprintf('L_mw = %g pH\n',L_mw*1e12)

%%
geometry = 'circle';

switch geometry
    case 'square'
        c1 = 1.27;
        c2 = 2.07;
        c3 = 0.18;
        c4 = 0.13;
    case 'circle'
        c1 = 1;
        c2 = 2.46;
        c3 = 0;
        c4 = 0.2;
end

L_gmd = 0.5*n^2*p.mu0*d_avg*c1*(log(c2/rho)+c3*rho+c4*rho^2);

fprintf([geometry ':\n'])
fprintf('w = %g um\n',w*1e6)
fprintf('s = %g um\n',s*1e6)
fprintf('d_out = %g um\n',d_out*1e6)
fprintf('d_in = %g um\n',d_in*1e6)
fprintf('n = %g \n',n)
fprintf('L_gmd = %g pH\n',L_gmd*1e12)
