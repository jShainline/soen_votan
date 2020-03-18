
close all
%open a figure with the data
uiopen('jjasinge_c5_2.fig',1)
%extract data from figure
[x,y]=ax_data(gca,1);

%%

v=x;%mV
i=y;%mA
%here I remove a small slope in the supercurrent branch due to a small
%resistance in the wires and the fact the we have a fake 4-wire measurement
ind=find(abs(i)<7e-3);
P=polyfit(i(ind),v(ind),1);
v=v-polyval(P,i);
plot(v,i);


%%
%for the initial guess of the fit
ind=find(v>0.01);
Ic=1e3*i(ind(1));
aux=(diff(v)./diff(i));
R=mean(aux(end-10:end));

T=10;
plot(v*1e3,i*1e3,IVthermal(i'*1e3,Ic,R,T),i*1e3,'r')
%% thermal fit
% Set up fittype and options.
ft = fittype( 'IVthermal(x,Ic,R,T)', 'independent', 'x','coefficient', {'Ic','R','T'});
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [Ic,R,T];
opts.Lower=[ 5  0  0];
opts.Upper=[ 20 R*2 50];

%fit wants things close to 1 so i scale the data to uA and uV
[xfit,yfit]=prepareCurveData(i*1e3,v*1e3);
% Fit model to data.
tic
[fitresult, gof] = fit( xfit, yfit, ft, opts);
toc
%%
figure(2)
plot(yfit,xfit,'.',fitresult(xfit),xfit,'linewidth',2)
axis([-Inf,Inf,-Inf,Inf])
%%
figure(3)
plot(v*1e3,i*1e3,IVthermal(i'*1e3,fitresult.Ic*1.1,fitresult.R,15),i*1e3,'r')

