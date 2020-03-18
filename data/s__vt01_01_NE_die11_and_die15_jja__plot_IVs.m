%%
clc;
close all;
% clear all
[~,bRGY] = f_colorSchemes('redWhiteBlue');
FontSize = 20;
FontName = 'Times';
scrsz = get(0,'ScreenSize');

%% switches
load_data = 'no';
save_plots = 'yes';
close_plots = 'yes';
run_thermal_fits = 'no';

%% load all
if strcmp(load_data,'yes')
    
    directory_list = {'20191125__vt01_01_NE_die11_jja','20191127__vt01_01_NE_die15_jja'};
    
    current_dir = pwd;
    data_map = containers.Map();
    for ii = 1:length(directory_list)
        cd(directory_list{ii})
        contents = dir;
        for jj = 3:length(contents)
            if strcmp(contents(jj).name(end-3:end),'.fig')
                open(contents(jj).name)
                h = get(gca,'children');
                x_data = get(h(length(h)),'xdata');
                y_data = get(h(length(h)),'ydata');
                close
                data_map(contents(jj).name(1:end-4)) = [x_data;y_data]';
            end
        end
        cd(current_dir)
    end
    
end

%% die 11, east side
data_list = {'vt01_01_NE_die11_jja_east_jj1_01','vt01_01_NE_die11_jja_east_jj2_01','vt01_01_NE_die11_jja_east_jj3_01','vt01_01_NE_die11_jja_east_jj4_01'};

%all on one
color_list = [3,8,13,18];
figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
for ii = 1:length(data_list)
    aa = data_map(data_list{ii});
    plot(aa(:,1)*1e3,aa(:,2)*1e3,'Color',bRGY(color_list(ii),:),'LineWidth',2);
    hold on
end
% ylim([0 80])
xlabel('Voltage [\mu V]','FontSize',FontSize,'FontName','Times')
ylabel('Current [\mu A]','FontSize',FontSize,'FontName','Times')
set(gca,'FontSize',FontSize,'FontName',FontName);
title(sprintf('I-V curves of JJs with different diameters\nvt01, NE quadrant, die 11, jja\neast side devices'),'FontSize',FontSize,'FontName','Times')
legend('d = 1.60\mum','d = 1.84\mum','d = 2.26\mum','d = 3.19\mum')
if strcmp(save_plots,'yes')
    saveas(gcf,'die_11_jja__east__jjIVs','png')
end
if strcmp(close_plots,'yes')
    close
end

%differential resistance
figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
for ii = 1:length(data_list)
    aa = data_map(data_list{ii});
    bb = diff(aa(:,1))./diff(aa(:,2));
    plot(aa(1:end-1,1),bb(:),'Color',bRGY(color_list(ii),:),'LineWidth',2);
    hold on
end
ylim([5 20])
xlabel('Voltage [\mu V]','FontSize',FontSize,'FontName','Times')
ylabel('Differential Resistance [\Delta V/\Delta I, \Omega]','FontSize',FontSize,'FontName','Times')
set(gca,'FontSize',FontSize,'FontName',FontName);
title(sprintf('Differential resistance of JJs with different diameters\nvt01, NE quadrant, die 11, jja\neast side devices'),'FontSize',FontSize,'FontName','Times')
legend('d = 1.60\mum','d = 1.84\mum','d = 2.26\mum','d = 3.19\mum')
if strcmp(save_plots,'yes')
    saveas(gcf,'die_11_jja__east__diffRes','png')
end
if strcmp(close_plots,'yes')
    close
end

%thermal fits
if strcmp(run_thermal_fits,'yes')
    Ic_guesses = [14e-6 22e-6 30e-6 60e-6];
    T_guesses = [10 10 10 10];
    Ic_fit = zeros(length(Ic_guesses),1);
    R_fit = zeros(length(Ic_guesses),1);
    T_fit = zeros(length(Ic_guesses),1);
    for kk = 1:length(data_list)
        
        aa = data_map(data_list{kk});
        vV = aa(:,1);
        iI = aa(:,2);
        
        %     ind = find( abs(iI) < 7e-3 );
        %     P = polyfit(iI(ind),vV(ind),1);
        %     vV = vV-polyval(P,iI);
        %     ind = find( vV > 0.01 );
        
        aux = diff(vV)./diff(iI);
        R_guess = mean(aux(end-10:end));
        
        ft = fittype( 'IVthermal(x,Ic,R,T)', 'independent', 'x','coefficient', {'Ic','R','T'});
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        opts.Display = 'Off';
        opts.StartPoint = [Ic_guesses(kk)*1e6,R_guess,T_guesses(kk)];
        opts.Lower = [0.8*Ic_guesses(kk)*1e6  0.1*R_guess  0.01*T_guesses(kk)];
        opts.Upper = [1.2*Ic_guesses(kk)*1e6 10*R_guess 10*T_guesses(kk)];
        
        %     ind = find( iI < 1.5*Ic_guesses(kk)*1e3 );
        %fit wants things close to 1, so scale the data to uA and uV
        [xfit,yfit] = prepareCurveData(iI(:)*1e3,vV(:)*1e3);
        [fitresult, gof] = fit(xfit,yfit,ft,opts);
        Ic_fit(kk) = fitresult.Ic;
        R_fit(kk) = fitresult.R;
        T_fit(kk) = fitresult.T;
        %     plot(yfit,xfit,'.',fitresult(xfit),xfit,'linewidth',2)
        figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
        plot(vV*1e3,iI*1e3,'Color',bRGY(3,:),'LineWidth',2)
        hold on
        plot(IVthermal(iI*1e3,fitresult.Ic,fitresult.R,fitresult.T)',iI*1e3,'Color',bRGY(8,:),'LineWidth',2)
        xlabel('Voltage [\mu V]','FontSize',FontSize,'FontName','Times')
        ylabel('Current [\mu A]','FontSize',FontSize,'FontName','Times')
        set(gca,'FontSize',FontSize,'FontName',FontName);
        title(sprintf('Thermal fit of JJ I-V\nvt01, NE quadrant, die 11, jja\neast side devices, jj%01d\nIc = %g uA, R = %g Ohm, T = %g K',kk,fitresult.Ic,fitresult.R,fitresult.T),'FontSize',FontSize,'FontName','Times')
        legend('Data','Fit')
        if strcmp(save_plots,'yes')
            saveas(gcf,sprintf('die_11_jja__east__thermal_fits__jj%01d',kk),'png')
        end
        if strcmp(close_plots,'yes')
            close
        end
    end
end

%Jc
diameter_vec = [1.6 1.84 2.26 3.19]*1e-6;
radius_vec = diameter_vec/2;
area_vec = (pi/4)*diameter_vec.^2;
Ic_vec = [13.7 20.5 29.2 59]*1e-6;%Ic_fit*1e-6;
Jc_vec = Ic_vec./area_vec;
figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
plot(area_vec*1e12,Jc_vec*1e-7,'Color',bRGY(color_list(ii),:),'LineWidth',2);
xlabel('Junction area [\mum^2]','FontSize',FontSize,'FontName','Times')
ylabel('Critical current density [kA/cm^2]','FontSize',FontSize,'FontName','Times')
title(sprintf('vt01, NE quadrant, die 11, jja\neast side devices'),'FontSize',FontSize,'FontName','Times')
set(gca,'FontSize',FontSize,'FontName',FontName);
if strcmp(save_plots,'yes')
    saveas(gcf,'die_11_jja__east__jc','png')
end
if strcmp(close_plots,'yes')
    close
end

p = polyfit(radius_vec,sqrt(Ic_vec),1);
% radius_vec_dense = linspace(radius_vec(1),radius_vec(end),100);
radius_vec_dense = linspace(0,radius_vec(end),100);
sqrt_Ic_vec_dense = polyval(p,radius_vec_dense);
slope = p(1);
Jc = 1e-7*slope^2/pi;
x_intercept = -p(2)/slope;
figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
plot(radius_vec_dense*1e6,sqrt_Ic_vec_dense*1e3,'Color',bRGY(3,:),'LineWidth',2);
hold on
plot(radius_vec*1e6,sqrt(Ic_vec)*1e3,'Color',bRGY(8,:),'LineWidth',2,'Marker','s','MarkerFaceColor',bRGY(6,:),'MarkerEdgeColor',bRGY(10,:),'MarkerSize',8);
% line([x_intercept x_intercept]*1e6,[min(sqrt_Ic_y_intercept_dense) max(sqrt_Ic_y_intercept_dense)]*1e3,'LineWidth',1,'LineStyle','-.','Color',bRGY(21,:))
xlabel('Junction radius [\mum]','FontSize',FontSize,'FontName','Times')
ylabel('$\sqrt{I_c}$ [$\mu$A$^{1/2}$]','FontSize',FontSize,'FontName','Times','interpreter','latex')
title(sprintf('Fit to obtain Jc\nvt01, NE quadrant, die 11, jja\neast side devices\nJc = slope$^2$/pi = %g kA/cm$^2$ (slope = %g uA$^{1/2}$/um)\nx-intercept = %g nm',Jc,1e-3*slope,x_intercept_alt*1e9),'FontSize',FontSize,'FontName','Times','interpreter','latex')
set(gca,'FontSize',FontSize,'FontName',FontName);
legend('fit','data')
ylim([0 1.1*max(sqrt_Ic_y_intercept_dense)*1e3])
xlim([0 1.1*radius_vec(end)*1e6])
if strcmp(save_plots,'yes')
    saveas(gcf,'die_11_jja__east__jc_fit','png')
end
if strcmp(close_plots,'yes')
    close
end

%% die 11, south side
data_list = {'vt01_01_NE_die11_jja_south_jj1_01','vt01_01_NE_die11_jja_south_jj2_01','vt01_01_NE_die11_jja_south_jj3_01','vt01_01_NE_die11_jja_south_jj4_01'};

%all on one
color_list = [3,8,13,18];
figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
for ii = 1:length(data_list)
    aa = data_map(data_list{ii});
    plot(aa(:,1)*1e3,aa(:,2)*1e3,'Color',bRGY(color_list(ii),:),'LineWidth',2);
    hold on
end
ylim([0 80])
ylabel('Current [\mu A]','FontSize',FontSize,'FontName','Times')
xlabel('Voltage [\mu V]','FontSize',FontSize,'FontName','Times')
set(gca,'FontSize',FontSize,'FontName',FontName);
title(sprintf('I-V curves of JJs with different diameters\nvt01, NE quadrant, die 11, jja\nsouth side devices'),'FontSize',FontSize,'FontName','Times')
legend('d = 1.60\mum','d = 1.84\mum','d = 2.26\mum','d = 3.19\mum')
if strcmp(save_plots,'yes')
    saveas(gcf,'die_11_jja__south__jjIVs','png')
end
if strcmp(close_plots,'yes')
    close
end

%differential resistance
figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
for ii = 1:length(data_list)
    aa = data_map(data_list{ii});
    bb = diff(aa(:,1))./diff(aa(:,2));
    plot(aa(1:end-1,1),bb(:),'Color',bRGY(color_list(ii),:),'LineWidth',2);
    hold on
end
ylim([5 12])
xlim([0 1])
ylabel('Differential Resistance [\Delta V/\Delta I, \Omega]','FontSize',FontSize,'FontName','Times')
xlabel('Voltage [\mu V]','FontSize',FontSize,'FontName','Times')
set(gca,'FontSize',FontSize,'FontName',FontName);
title(sprintf('Differential resistance of JJs with different diameters\nvt01, NE quadrant, die 11, jja\nsouth side devices'),'FontSize',FontSize,'FontName','Times')
legend('d = 1.60\mum','d = 1.84\mum','d = 2.26\mum','d = 3.19\mum')
if strcmp(save_plots,'yes')
    saveas(gcf,'die_11_jja__east__diffRes','png')
end
if strcmp(close_plots,'yes')
    close
end

%thermal fits
if strcmp(run_thermal_fits,'yes')
    Ic_guesses = [11.7e-6 17.7e-6 29.2e-6 59.2e-6];
    T_guesses = [10 10 10 10];
    Ic_fit = zeros(length(Ic_guesses),1);
    R_fit = zeros(length(Ic_guesses),1);
    T_fit = zeros(length(Ic_guesses),1);
    for kk = 1:length(data_list)
        
        aa = data_map(data_list{kk});
        vV = aa(:,1);
        iI = aa(:,2);
        
        %     ind = find( abs(iI) < 7e-3 );
        %     P = polyfit(iI(ind),vV(ind),1);
        %     vV = vV-polyval(P,iI);
        %     ind = find( vV > 0.01 );
        
        aux = diff(vV)./diff(iI);
        R_guess = mean(aux(end-10:end));
        
        ft = fittype( 'IVthermal(x,Ic,R,T)', 'independent', 'x','coefficient', {'Ic','R','T'});
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        opts.Display = 'Off';
        opts.StartPoint = [Ic_guesses(kk)*1e6,R_guess,T_guesses(kk)];
        opts.Lower = [0.8*Ic_guesses(kk)*1e6  0.1*R_guess  0.01*T_guesses(kk)];
        opts.Upper = [1.2*Ic_guesses(kk)*1e6 10*R_guess 10*T_guesses(kk)];
        
        %     ind = find( iI < 1.5*Ic_guesses(kk)*1e3 );
        %fit wants things close to 1, so scale the data to uA and uV
        [xfit,yfit] = prepareCurveData(iI(:)*1e3,vV(:)*1e3);
        [fitresult, gof] = fit(xfit,yfit,ft,opts);
        Ic_fit(kk) = fitresult.Ic;
        R_fit(kk) = fitresult.R;
        T_fit(kk) = fitresult.T;
        %     plot(yfit,xfit,'.',fitresult(xfit),xfit,'linewidth',2)
        figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
        plot(vV*1e3,iI*1e3,'Color',bRGY(3,:),'LineWidth',2)
        hold on
        plot(IVthermal(iI*1e3,fitresult.Ic,fitresult.R,fitresult.T)',iI*1e3,'Color',bRGY(8,:),'LineWidth',2)
        xlabel('Voltage [\mu V]','FontSize',FontSize,'FontName','Times')
        ylabel('Current [\mu A]','FontSize',FontSize,'FontName','Times')
        set(gca,'FontSize',FontSize,'FontName',FontName);
        title(sprintf('Thermal fit of JJ I-V\nvt01, NE quadrant, die 11, jja\neast side devices, jj%01d\nIc = %g uA, R = %g Ohm, T = %g K',kk,fitresult.Ic,fitresult.R,fitresult.T),'FontSize',FontSize,'FontName','Times')
        legend('Data','Fit')
        if strcmp(save_plots,'yes')
            saveas(gcf,sprintf('die_11_jja__south__thermal_fits__jj%01d',kk),'png')
        end
        if strcmp(close_plots,'yes')
            close
        end
    end
end

%Jc
diameter_vec = [1.6 1.84 2.26 3.19]*1e-6;
area_vec = (pi/4)*diameter_vec.^2;
Ic_vec = [11.7 17.7 29.2 59.2]*1e-6;
Jc_vec = Ic_vec./area_vec;
figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
plot(area_vec*1e12,Jc_vec*1e-7,'Color',bRGY(color_list(ii),:),'LineWidth',2);
ylabel('Critical current density [kA/cm^2]','FontSize',FontSize,'FontName','Times')
xlabel('Junction area [\mum^2]','FontSize',FontSize,'FontName','Times')
title(sprintf('vt01, NE quadrant, die 11, jja\nsouth side devices'),'FontSize',FontSize,'FontName','Times')
set(gca,'FontSize',FontSize,'FontName',FontName);
if strcmp(save_plots,'yes')
    saveas(gcf,'die_11_jja__south__jc','png')
end
if strcmp(close_plots,'yes')
    close
end

p = polyfit(radius_vec,sqrt(Ic_vec),1);
% radius_vec_dense = linspace(radius_vec(1),radius_vec(end),100);
radius_vec_dense = linspace(0,radius_vec(end),100);
sqrt_Ic_vec_dense = polyval(p,radius_vec_dense);
slope = p(1);
Jc = 1e-7*slope^2/pi;
x_intercept = -p(2)/slope;
figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
plot(radius_vec_dense*1e6,sqrt_Ic_vec_dense*1e3,'Color',bRGY(3,:),'LineWidth',2);
hold on
plot(radius_vec*1e6,sqrt(Ic_vec)*1e3,'Color',bRGY(8,:),'LineWidth',2,'Marker','s','MarkerFaceColor',bRGY(6,:),'MarkerEdgeColor',bRGY(10,:),'MarkerSize',8);
% line([x_intercept x_intercept]*1e6,[min(sqrt_Ic_y_intercept_dense) max(sqrt_Ic_y_intercept_dense)]*1e3,'LineWidth',1,'LineStyle','-.','Color',bRGY(21,:))
xlabel('Junction radius [\mum]','FontSize',FontSize,'FontName','Times')
ylabel('$\sqrt{I_c}$ [$\mu$A$^{1/2}$]','FontSize',FontSize,'FontName','Times','interpreter','latex')
title(sprintf('Fit to obtain Jc\nvt01, NE quadrant, die 11, jja\nsouth side devices\nJc = slope$^2$/pi = %g kA/cm$^2$ (slope = %g uA$^{1/2}$/um)\nx-intercept = %g nm',Jc,1e-3*slope,x_intercept*1e9),'FontSize',FontSize,'FontName','Times','interpreter','latex')
set(gca,'FontSize',FontSize,'FontName',FontName);
legend('fit','data')
ylim([0 1.1*max(sqrt_Ic_y_intercept_dense)*1e3])
xlim([0 1.1*radius_vec(end)*1e6])
if strcmp(save_plots,'yes')
    saveas(gcf,'die_11_jja__south__jc_fit','png')
end
if strcmp(close_plots,'yes')
    close
end

%% die 11, south and east side together
data_list_1 = {'vt01_01_NE_die11_jja_east_jj1_01','vt01_01_NE_die11_jja_east_jj2_01','vt01_01_NE_die11_jja_east_jj3_01','vt01_01_NE_die11_jja_east_jj4_01'};
data_list_2 = {'vt01_01_NE_die11_jja_south_jj1_01','vt01_01_NE_die11_jja_south_jj2_01','vt01_01_NE_die11_jja_south_jj3_01','vt01_01_NE_die11_jja_south_jj4_01'};

color_list_1 = [1,2,3,4];
color_list_2 = [6,7,8,9];

%all on one
figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
for ii = 1:length(data_list_1)
    aa = data_map(data_list_1{ii});
    plot(aa(:,1)*1e3,aa(:,2)*1e3,'Color',bRGY(color_list_1(ii),:),'LineWidth',2);
    hold on
end
for ii = 1:length(data_list_2)
    bb = data_map(data_list_2{ii});
    plot(bb(:,1)*1e3,bb(:,2)*1e3,'Color',bRGY(color_list_2(ii),:),'LineWidth',2);
    hold on
end
ylabel('Current [\mu A]','FontSize',FontSize,'FontName','Times')
xlabel('Voltage [\mu V]','FontSize',FontSize,'FontName','Times')
set(gca,'FontSize',FontSize,'FontName',FontName);
title(sprintf('I-V curves of JJs with different shunt resistors\nvt01, NE quadrant, die 11, jja\ncomparing east and south side devices'),'FontSize',FontSize,'FontName','Times')
legend('jj diameter = 1.60um, resistor = 2.05 squares','jj diameter = 1.84um, resistor = 2.05 squares','jj diameter = 2.26um, resistor = 2.05 squares','jj diameter = 3.19um, resistor = 2.05 squares','jj diameter = 1.60um, resistor = 1.15 squares','jj diameter = 1.84um, resistor = 1.15 squares','jj diameter = 2.26um, resistor = 1.15 squares','jj diameter = 3.19um, resistor = 1.15 squares')
if strcmp(save_plots,'yes')
    saveas(gcf,'die_11_jja__south_and_east__jjIVs_all','png')
end
if strcmp(close_plots,'yes')
    close
end

%each jj diameter separately
title_string_list = {'jj diameter = 1.60um','jj diameter = 1.84um','jj diameter = 2.26um','jj diameter = 3.19um'};
for ii = 1:length(data_list_1)
    figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
    aa = data_map(data_list_1{ii});
    bb = data_map(data_list_2{ii});
    plot(aa(:,1)*1e3,aa(:,2)*1e3,'Color',bRGY(3,:),'LineWidth',2);
    hold on
    plot(bb(:,1)*1e3,bb(:,2)*1e3,'Color',bRGY(8,:),'LineWidth',2);
    %     ylim([0 80])
    ylabel('Current [\mu A]','FontSize',FontSize,'FontName','Times')
    xlabel('Voltage [\mu V]','FontSize',FontSize,'FontName','Times')
    set(gca,'FontSize',FontSize,'FontName',FontName);
    title(sprintf(['I-V curves of JJs with different shunt resistors\nvt01, NE quadrant, die 11, jja\ncomparing east and south side devices, ',title_string_list{ii}]),'FontSize',FontSize,'FontName','Times')
    legend('resistor = 2.05 squares','resistor = 1.15 squares')
    if strcmp(save_plots,'yes')
        saveas(gcf,sprintf('die_11_jja__south_and_east__jjIVs_jj%01d',ii),'png')
    end
    if strcmp(close_plots,'yes')
        close
    end
end

%differential resistance
color_list_1 = [1,2,3,4];
color_list_2 = [6,7,8,9];
figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
for ii = 1:length(data_list_1)
    aa = data_map(data_list_1{ii});
    bb = diff(aa(:,1))./diff(aa(:,2));
    plot(aa(1:end-1,1),bb(:),'Color',bRGY(color_list_1(ii),:),'LineWidth',2);
    hold on
    aa = data_map(data_list_2{ii});
    bb = diff(aa(:,1))./diff(aa(:,2));
    plot(aa(1:end-1,1),bb(:),'Color',bRGY(color_list_2(ii),:),'LineWidth',2);
    hold on
end
ylim([0 20])
xlim([0 1.5])
ylabel('Differential Resistance [\Delta V/\Delta I, \Omega]','FontSize',FontSize,'FontName','Times')
xlabel('Voltage [\mu V]','FontSize',FontSize,'FontName','Times')
set(gca,'FontSize',FontSize,'FontName',FontName);
title(sprintf('Differential resistance of JJs with different diameters\nvt01, NE quadrant, die 11, jja\nsouth side devices'),'FontSize',FontSize,'FontName','Times')
legend('d = 1.60\mum, east','d = 1.60\mum, south','d = 1.84\mum, east','d = 1.84\mum, south','d = 2.26\mum, east','d = 2.26\mum, south','d = 3.19\mum, east','d = 3.19\mum, south')
if strcmp(save_plots,'yes')
    saveas(gcf,'die_11_jja__east_and_south__diffRes','png')
end
if strcmp(close_plots,'yes')
    close
end

%% die 11, north side

%series 80
data_list = {'vt01_01_NE_die11_jja_north_series80_jj1_01','vt01_01_NE_die11_jja_north_series80_jj2_01','vt01_01_NE_die11_jja_north_series80_jj3_01','vt01_01_NE_die11_jja_north_series80_jj4_01'};
color_list = [3,8,13,18];
figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
for ii = 1:length(data_list)
    aa = data_map(data_list{ii});
    plot(aa(:,1),aa(:,2)*1e3,'Color',bRGY(color_list(ii),:),'LineWidth',2);
    hold on
end
ylabel('Current [\mu A]','FontSize',FontSize,'FontName','Times')
xlabel('Voltage [mV]','FontSize',FontSize,'FontName','Times')
set(gca,'FontSize',FontSize,'FontName',FontName);
title(sprintf('I-V curves of series arrays of JJs with different diameters\nvt01, NE quadrant, die 11, jja\nnorth side devices\n80 jjs in series'),'FontSize',FontSize,'FontName','Times')
legend('d = 1.60\mum','d = 1.84\mum','d = 2.26\mum','d = 3.19\mum')
if strcmp(save_plots,'yes')
    saveas(gcf,'die_11_jja__north__series80__jjIVs','png')
end
if strcmp(close_plots,'yes')
    close
end

%comparing two that should be identical
data_list = {'vt01_01_NE_die11_jja_north_series80_jj4_01','vt01_01_NE_die11_jja_north_series80_jj4_02'};
color_list = [3,8,13,18];
figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
for ii = 1:length(data_list)
    aa = data_map(data_list{ii});
    plot(aa(:,1),aa(:,2)*1e3,'Color',bRGY(color_list(ii),:),'LineWidth',2);
    hold on
end
ylabel('Current [\mu A]','FontSize',FontSize,'FontName','Times')
xlabel('Voltage [mV]','FontSize',FontSize,'FontName','Times')
set(gca,'FontSize',FontSize,'FontName',FontName);
title(sprintf('I-V curves of 80 JJs with the same diameter\nvt01, NE quadrant, die 11, jja\nnorth side devices\n80 jjs in series'),'FontSize',FontSize,'FontName','Times')
legend('I-V sweep 1','I-V sweep 2')
if strcmp(save_plots,'yes')
    saveas(gcf,'die_11_jja__north__series80__jjIVs__two_sweeps','png')
end
if strcmp(close_plots,'yes')
    close
end

%series 160
data_list = {'vt01_01_NE_die11_jja_north_series160_jj1_01','vt01_01_NE_die11_jja_north_series160_jj2_01','vt01_01_NE_die11_jja_north_series160_jj3_01','vt01_01_NE_die11_jja_north_series160_jj4_01'};
color_list = [3,8,13,18];
figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
for ii = 1:length(data_list)
    aa = data_map(data_list{ii});
    plot(aa(:,1),aa(:,2)*1e3,'Color',bRGY(color_list(ii),:),'LineWidth',2);
    hold on
end
ylabel('Current [\mu A]','FontSize',FontSize,'FontName','Times')
xlabel('Voltage [mV]','FontSize',FontSize,'FontName','Times')
set(gca,'FontSize',FontSize,'FontName',FontName);
title(sprintf('I-V curves of series arrays of JJs with different diameters\nvt01, NE quadrant, die 11, jja\nnorth side devices\n160 jjs in series'),'FontSize',FontSize,'FontName','Times')
legend('d = 1.60\mum','d = 1.84\mum','d = 2.26\mum','d = 3.19\mum')
if strcmp(save_plots,'yes')
    saveas(gcf,'die_11_jja__north__series160__jjIVs','png')
end
if strcmp(close_plots,'yes')
    close
end

%% die 11, west side

%series 80
data_list = {'vt01_01_NE_die11_jja_west_series80_jj1_01','vt01_01_NE_die11_jja_west_series80_jj2_01','vt01_01_NE_die11_jja_west_series80_jj3_01','vt01_01_NE_die11_jja_west_series80_jj4_01'};
color_list = [3,8,13,18];
figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
for ii = 1:length(data_list)
    aa = data_map(data_list{ii});
    plot(aa(:,1),aa(:,2)*1e3,'Color',bRGY(color_list(ii),:),'LineWidth',2);
    hold on
end
ylabel('Current [\mu A]','FontSize',FontSize,'FontName','Times')
xlabel('Voltage [mV]','FontSize',FontSize,'FontName','Times')
set(gca,'FontSize',FontSize,'FontName',FontName);
title(sprintf('I-V curves of series arrays of JJs with different diameters\nvt01, NE quadrant, die 11, jja\nwest side devices\n80 jjs in series'),'FontSize',FontSize,'FontName','Times')
legend('d = 1.60\mum','d = 1.84\mum','d = 2.26\mum','d = 3.19\mum')
if strcmp(save_plots,'yes')
    saveas(gcf,'die_11_jja__west__series80__jjIVs','png')
end
if strcmp(close_plots,'yes')
    close
end

%series 160
data_list = {'vt01_01_NE_die11_jja_west_series160_jj1_01','vt01_01_NE_die11_jja_west_series160_jj2_01','vt01_01_NE_die11_jja_west_series160_jj3_01','vt01_01_NE_die11_jja_west_series160_jj4_01'};
color_list = [3,8,13,18];
figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
for ii = 1:length(data_list)
    aa = data_map(data_list{ii});
    plot(aa(:,1),aa(:,2)*1e3,'Color',bRGY(color_list(ii),:),'LineWidth',2);
    hold on
end
ylabel('Current [\mu A]','FontSize',FontSize,'FontName','Times')
xlabel('Voltage [mV]','FontSize',FontSize,'FontName','Times')
set(gca,'FontSize',FontSize,'FontName',FontName);
title(sprintf('I-V curves of series arrays of JJs with different diameters\nvt01, NE quadrant, die 11, jja\nwest side devices\n160 jjs in series'),'FontSize',FontSize,'FontName','Times')
legend('d = 1.60\mum','d = 1.84\mum','d = 2.26\mum','d = 3.19\mum')
if strcmp(save_plots,'yes')
    saveas(gcf,'die_11_jja__west__series160__jjIVs','png')
end
if strcmp(close_plots,'yes')
    close
end

%comparing 80 to 160
data_list_1 = {'vt01_01_NE_die11_jja_west_series80_jj1_01','vt01_01_NE_die11_jja_west_series80_jj2_01','vt01_01_NE_die11_jja_west_series80_jj3_01','vt01_01_NE_die11_jja_west_series80_jj4_01'};
data_list_2 = {'vt01_01_NE_die11_jja_west_series160_jj1_01','vt01_01_NE_die11_jja_west_series160_jj2_01','vt01_01_NE_die11_jja_west_series160_jj3_01','vt01_01_NE_die11_jja_west_series160_jj4_01'};

title_string_list = {'jj diameter = 1.60um','jj diameter = 1.84um','jj diameter = 2.26um','jj diameter = 3.19um'};
for ii = 1:length(data_list_1)
    figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
    aa = data_map(data_list_1{ii});
    bb = data_map(data_list_2{ii});
    plot(aa(:,1)*1e3,aa(:,2)*1e3,'Color',bRGY(3,:),'LineWidth',2);
    hold on
    plot(bb(:,1)*1e3,bb(:,2)*1e3,'Color',bRGY(8,:),'LineWidth',2);
    %     ylim([0 80])
    ylabel('Current [\mu A]','FontSize',FontSize,'FontName','Times')
    xlabel('Voltage [\mu V]','FontSize',FontSize,'FontName','Times')
    set(gca,'FontSize',FontSize,'FontName',FontName);
    title(sprintf(['I-V curves of series arrays of JJs with different shunt resistors\nvt01, NE quadrant, die 11, jja\nwest side devices, comparing series arrays of 80 to 160\n',title_string_list{ii}]),'FontSize',FontSize,'FontName','Times')
    legend('80 jjs in series','160 jjs in series')
    if strcmp(save_plots,'yes')
        saveas(gcf,sprintf('die_11_jja__west__series__jjIVs_jj%01d',ii),'png')
    end
    if strcmp(close_plots,'yes')
        close
    end
end

%% die 15, east side
data_list = {'vt01_01_NE_die15_jja_east_jj1_01','vt01_01_NE_die15_jja_east_jj2_01','vt01_01_NE_die15_jja_east_jj3_01','vt01_01_NE_die15_jja_east_jj4_01'};

%all on one
color_list = [3,8,13,18];
figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
for ii = 1:length(data_list)
    aa = data_map(data_list{ii});
    plot(aa(:,1)*1e3,aa(:,2)*1e3,'Color',bRGY(color_list(ii),:),'LineWidth',2);
    hold on
end
ylim([0 80])
ylabel('Current [\mu A]','FontSize',FontSize,'FontName','Times')
xlabel('Voltage [\mu V]','FontSize',FontSize,'FontName','Times')
set(gca,'FontSize',FontSize,'FontName',FontName);
title(sprintf('I-V curves of JJs with different diameters\nvt01, NE quadrant, die 15, jja\neast side devices'),'FontSize',FontSize,'FontName','Times')
legend('d = 1.60\mum','d = 1.84\mum','d = 2.26\mum','d = 3.19\mum')
if strcmp(save_plots,'yes')
    saveas(gcf,'die_15_jja__east__jjIVs','png')
end
if strcmp(close_plots,'yes')
    close
end

%differential resistance
figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
for ii = 1:length(data_list)
    aa = data_map(data_list{ii});
    bb = diff(aa(:,1))./diff(aa(:,2));
    plot(aa(1:end-1,1),bb(:),'Color',bRGY(color_list(ii),:),'LineWidth',2);
    hold on
end
ylim([-60 100])
xlim([0 1])
ylabel('Differential Resistance [\Delta V/\Delta I, \Omega]','FontSize',FontSize,'FontName','Times')
xlabel('Voltage [\mu V]','FontSize',FontSize,'FontName','Times')
set(gca,'FontSize',FontSize,'FontName',FontName);
title(sprintf('Differential resistance of JJs with different diameters\nvt01, NE quadrant, die 15, jja\neast side devices'),'FontSize',FontSize,'FontName','Times')
legend('d = 1.60\mum','d = 1.84\mum','d = 2.26\mum','d = 3.19\mum')
if strcmp(save_plots,'yes')
    saveas(gcf,'die_15_jja__east__diffRes','png')
end
if strcmp(close_plots,'yes')
    close
end

%thermal fits
if strcmp(run_thermal_fits,'yes')
    Ic_guesses = [15.6 24 30.5 67.5]*1e-6;
    T_guesses = [10 10 10 10];
    Ic_fit = zeros(length(Ic_guesses),1);
    R_fit = zeros(length(Ic_guesses),1);
    T_fit = zeros(length(Ic_guesses),1);
    for kk = 1:length(data_list)
        
        aa = data_map(data_list{kk});
        vV = aa(:,1);
        iI = aa(:,2);
        
        %     ind = find( abs(iI) < 7e-3 );
        %     P = polyfit(iI(ind),vV(ind),1);
        %     vV = vV-polyval(P,iI);
        %     ind = find( vV > 0.01 );
        
        aux = diff(vV)./diff(iI);
        R_guess = mean(aux(end-10:end));
        
        ft = fittype( 'IVthermal(x,Ic,R,T)', 'independent', 'x','coefficient', {'Ic','R','T'});
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        opts.Display = 'Off';
        opts.StartPoint = [Ic_guesses(kk)*1e6,R_guess,T_guesses(kk)];
        opts.Lower = [0.8*Ic_guesses(kk)*1e6  0.1*R_guess  0.01*T_guesses(kk)];
        opts.Upper = [1.2*Ic_guesses(kk)*1e6 10*R_guess 10*T_guesses(kk)];
        
        %     ind = find( iI < 1.5*Ic_guesses(kk)*1e3 );
        %fit wants things close to 1, so scale the data to uA and uV
        [xfit,yfit] = prepareCurveData(iI(:)*1e3,vV(:)*1e3);
        [fitresult, gof] = fit(xfit,yfit,ft,opts);
        Ic_fit(kk) = fitresult.Ic;
        R_fit(kk) = fitresult.R;
        T_fit(kk) = fitresult.T;
        %     plot(yfit,xfit,'.',fitresult(xfit),xfit,'linewidth',2)
        figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
        plot(vV*1e3,iI*1e3,'Color',bRGY(3,:),'LineWidth',2)
        hold on
        plot(IVthermal(iI*1e3,fitresult.Ic,fitresult.R,fitresult.T)',iI*1e3,'Color',bRGY(8,:),'LineWidth',2)
        xlabel('Voltage [\mu V]','FontSize',FontSize,'FontName','Times')
        ylabel('Current [\mu A]','FontSize',FontSize,'FontName','Times')
        set(gca,'FontSize',FontSize,'FontName',FontName);
        title(sprintf('Thermal fit of JJ I-V\nvt01, NE quadrant, die 11, jja\neast side devices, jj%01d\nIc = %g uA, R = %g Ohm, T = %g K',kk,fitresult.Ic,fitresult.R,fitresult.T),'FontSize',FontSize,'FontName','Times')
        legend('Data','Fit')
        if strcmp(save_plots,'yes')
            saveas(gcf,sprintf('die_15_jja__east__thermal_fits__jj%01d',kk),'png')
        end
        if strcmp(close_plots,'yes')
            close
        end
    end
end

%Jc
diameter_vec = [1.6 1.84 2.26 3.19]*1e-6;
area_vec = (pi/4)*diameter_vec.^2;
Ic_vec = [15.5 24 30.5 67.5]*1e-6;
Jc_vec = Ic_vec./area_vec;
figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
plot(area_vec*1e12,Jc_vec*1e-7,'Color',bRGY(color_list(ii),:),'LineWidth',2);
ylabel('Critical current density [kA/cm^2]','FontSize',FontSize,'FontName','Times')
xlabel('Junction area [\mum^2]','FontSize',FontSize,'FontName','Times')
title(sprintf('vt01, NE quadrant, die 15, jja\neast side devices'),'FontSize',FontSize,'FontName','Times')
set(gca,'FontSize',FontSize,'FontName',FontName);
if strcmp(save_plots,'yes')
    saveas(gcf,'die_15_jja__east__jc','png')
end
if strcmp(close_plots,'yes')
    close
end

p = polyfit(radius_vec,sqrt(Ic_vec),1);
% radius_vec_dense = linspace(radius_vec(1),radius_vec(end),100);
radius_vec_dense = linspace(0,radius_vec(end),100);
sqrt_Ic_vec_dense = polyval(p,radius_vec_dense);
slope = p(1);
Jc = 1e-7*slope^2/pi;
x_intercept = -p(2)/slope;
figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
plot(radius_vec_dense*1e6,sqrt_Ic_vec_dense*1e3,'Color',bRGY(3,:),'LineWidth',2);
hold on
plot(radius_vec*1e6,sqrt(Ic_vec)*1e3,'Color',bRGY(8,:),'LineWidth',2,'Marker','s','MarkerFaceColor',bRGY(6,:),'MarkerEdgeColor',bRGY(10,:),'MarkerSize',8);
% line([x_intercept x_intercept]*1e6,[min(sqrt_Ic_y_intercept_dense) max(sqrt_Ic_y_intercept_dense)]*1e3,'LineWidth',1,'LineStyle','-.','Color',bRGY(21,:))
xlabel('Junction radius [\mum]','FontSize',FontSize,'FontName','Times')
ylabel('$\sqrt{I_c}$ [$\mu$A$^{1/2}$]','FontSize',FontSize,'FontName','Times','interpreter','latex')
title(sprintf('Fit to obtain Jc\nvt01, NE quadrant, die 15, jja\neast side devices\nJc = slope$^2$/pi = %g kA/cm$^2$ (slope = %g uA$^{1/2}$/um)\nx-intercept = %g nm',Jc,1e-3*slope,x_intercept*1e9),'FontSize',FontSize,'FontName','Times','interpreter','latex')
set(gca,'FontSize',FontSize,'FontName',FontName);
legend('fit','data')
ylim([0 1.1*max(sqrt_Ic_y_intercept_dense)*1e3])
xlim([0 1.1*radius_vec(end)*1e6])
if strcmp(save_plots,'yes')
    saveas(gcf,'die_15_jja__east__jc_fit','png')
end
if strcmp(close_plots,'yes')
    close
end

%% die 15, east side, hysteresis
data_list = {'vt01_01_NE_die15_jja_east_jj1_02','vt01_01_NE_die15_jja_east_jj2_02','vt01_01_NE_die15_jja_east_jj3_02','vt01_01_NE_die15_jja_east_jj4_02'};
color_list = [3,8,13,18];
title_str = {'d = 1.60um','d = 1.84um','d = 2.26um','d = 3.19um'};
for ii = 1:length(data_list)
    figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
    aa = data_map(data_list{ii});
    plot(aa(:,1)*1e3,aa(:,2)*1e3,'Color',bRGY(color_list(ii),:),'LineWidth',2);
    %     ylim([0 80])
    ylabel('Current [\mu A]','FontSize',FontSize,'FontName','Times')
    xlabel('Voltage [\mu V]','FontSize',FontSize,'FontName','Times')
    set(gca,'FontSize',FontSize,'FontName',FontName);
    title(sprintf(['I-V curve of JJ with ',title_str{ii},'\nforward and reverse sweep\nvt01, NE quadrant, die 15, jja\neast side devices']),'FontSize',FontSize,'FontName','Times')
    if strcmp(save_plots,'yes')
        saveas(gcf,sprintf('die_15_jja__east__jjIVs_hysteresis__jj%01d',ii),'png')
    end
    if strcmp(close_plots,'yes')
        close
    end
end

%% die 15, south side

data_list = {'vt01_01_NE_die15_jja_south_jj1_01','vt01_01_NE_die15_jja_south_jj2_01','vt01_01_NE_die15_jja_south_jj3_01','vt01_01_NE_die15_jja_south_jj4_01'};

%all on one
color_list = [3,8,13,18];
figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
for ii = 1:length(data_list)
    aa = data_map(data_list{ii});
    plot(aa(:,1)*1e3,aa(:,2)*1e3,'Color',bRGY(color_list(ii),:),'LineWidth',2);
    hold on
end
ylim([0 80])
ylabel('Current [\mu A]','FontSize',FontSize,'FontName','Times')
xlabel('Voltage [\mu V]','FontSize',FontSize,'FontName','Times')
set(gca,'FontSize',FontSize,'FontName',FontName);
title(sprintf('I-V curves of JJs with different diameters\nvt01, NE quadrant, die 15, jja\nsouth side devices'),'FontSize',FontSize,'FontName','Times')
legend('d = 1.60\mum','d = 1.84\mum','d = 2.26\mum','d = 3.19\mum')
if strcmp(save_plots,'yes')
    saveas(gcf,'die_15_jja__south__jjIVs','png')
end
if strcmp(close_plots,'yes')
    close
end

%differential resistance
figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
for ii = 1:length(data_list)
    aa = data_map(data_list{ii});
    bb = diff(aa(:,1))./diff(aa(:,2));
    plot(aa(1:end-1,1),bb(:),'Color',bRGY(color_list(ii),:),'LineWidth',2);
    hold on
end
ylim([-60 100])
xlim([0 1])
ylabel('Differential Resistance [\Delta V/\Delta I, \Omega]','FontSize',FontSize,'FontName','Times')
xlabel('Voltage [\mu V]','FontSize',FontSize,'FontName','Times')
set(gca,'FontSize',FontSize,'FontName',FontName);
title(sprintf('Differential resistance of JJs with different diameters\nvt01, NE quadrant, die 15, jja\nsouth side devices'),'FontSize',FontSize,'FontName','Times')
legend('d = 1.60\mum','d = 1.84\mum','d = 2.26\mum','d = 3.19\mum')
if strcmp(save_plots,'yes')
    saveas(gcf,'die_15_jja__south__diffRes','png')
end
if strcmp(close_plots,'yes')
    close
end

%thermal fits
if strcmp(run_thermal_fits,'yes')
    Ic_guesses = [15.6 24 30.5 67.5]*1e-6;
    T_guesses = [10 10 10 10];
    Ic_fit = zeros(length(Ic_guesses),1);
    R_fit = zeros(length(Ic_guesses),1);
    T_fit = zeros(length(Ic_guesses),1);
    for kk = 1:length(data_list)
        
        aa = data_map(data_list{kk});
        vV = aa(:,1);
        iI = aa(:,2);
        
        %     ind = find( abs(iI) < 7e-3 );
        %     P = polyfit(iI(ind),vV(ind),1);
        %     vV = vV-polyval(P,iI);
        %     ind = find( vV > 0.01 );
        
        aux = diff(vV)./diff(iI);
        R_guess = mean(aux(end-10:end));
        
        ft = fittype( 'IVthermal(x,Ic,R,T)', 'independent', 'x','coefficient', {'Ic','R','T'});
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        opts.Display = 'Off';
        opts.StartPoint = [Ic_guesses(kk)*1e6,R_guess,T_guesses(kk)];
        opts.Lower = [0.8*Ic_guesses(kk)*1e6  0.1*R_guess  0.01*T_guesses(kk)];
        opts.Upper = [1.2*Ic_guesses(kk)*1e6 10*R_guess 10*T_guesses(kk)];
        
        %     ind = find( iI < 1.5*Ic_guesses(kk)*1e3 );
        %fit wants things close to 1, so scale the data to uA and uV
        [xfit,yfit] = prepareCurveData(iI(:)*1e3,vV(:)*1e3);
        [fitresult, gof] = fit(xfit,yfit,ft,opts);
        Ic_fit(kk) = fitresult.Ic;
        R_fit(kk) = fitresult.R;
        T_fit(kk) = fitresult.T;
        %     plot(yfit,xfit,'.',fitresult(xfit),xfit,'linewidth',2)
        figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
        plot(vV*1e3,iI*1e3,'Color',bRGY(3,:),'LineWidth',2)
        hold on
        plot(IVthermal(iI*1e3,fitresult.Ic,fitresult.R,fitresult.T)',iI*1e3,'Color',bRGY(8,:),'LineWidth',2)
        xlabel('Voltage [\mu V]','FontSize',FontSize,'FontName','Times')
        ylabel('Current [\mu A]','FontSize',FontSize,'FontName','Times')
        set(gca,'FontSize',FontSize,'FontName',FontName);
        title(sprintf('Thermal fit of JJ I-V\nvt01, NE quadrant, die 11, jja\neast side devices, jj%01d\nIc = %g uA, R = %g Ohm, T = %g K',kk,fitresult.Ic,fitresult.R,fitresult.T),'FontSize',FontSize,'FontName','Times')
        legend('Data','Fit')
        if strcmp(save_plots,'yes')
            saveas(gcf,sprintf('die_15_jja__south__thermal_fits__jj%01d',kk),'png')
        end
        if strcmp(close_plots,'yes')
            close
        end
    end
end

%Jc
diameter_vec = [1.6 1.84 2.26 3.19]*1e-6;
area_vec = (pi/4)*diameter_vec.^2;
Ic_vec = [13.7 19.5 33 61.5]*1e-6;

Jc_vec = Ic_vec./area_vec;
h = figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
plot(area_vec*1e12,Jc_vec*1e-7,'Color',bRGY(color_list(ii),:),'LineWidth',2);
ylabel('Critical current density [kA/cm^2]','FontSize',FontSize,'FontName','Times')
xlabel('Junction area [\mum^2]','FontSize',FontSize,'FontName','Times')
title(sprintf('vt01, NE quadrant, die 15, jja\nsouth side devices'),'FontSize',FontSize,'FontName','Times')
set(gca,'FontSize',FontSize,'FontName',FontName);
saveas(gcf,'die_15_jja__south__jc','png')
if strcmp(close_plots,'yes')
    close
end

p = polyfit(radius_vec,sqrt(Ic_vec),1);
% radius_vec_dense = linspace(radius_vec(1),radius_vec(end),100);
radius_vec_dense = linspace(0,radius_vec(end),100);
sqrt_Ic_vec_dense = polyval(p,radius_vec_dense);
slope = p(1);
Jc = 1e-7*slope^2/pi;
x_intercept = -p(2)/slope;
figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
plot(radius_vec_dense*1e6,sqrt_Ic_vec_dense*1e3,'Color',bRGY(3,:),'LineWidth',2);
hold on
plot(radius_vec*1e6,sqrt(Ic_vec)*1e3,'Color',bRGY(8,:),'LineWidth',2,'Marker','s','MarkerFaceColor',bRGY(6,:),'MarkerEdgeColor',bRGY(10,:),'MarkerSize',8);
% line([x_intercept x_intercept]*1e6,[min(sqrt_Ic_y_intercept_dense) max(sqrt_Ic_y_intercept_dense)]*1e3,'LineWidth',1,'LineStyle','-.','Color',bRGY(21,:))
xlabel('Junction radius [\mum]','FontSize',FontSize,'FontName','Times')
ylabel('$\sqrt{I_c}$ [$\mu$A$^{1/2}$]','FontSize',FontSize,'FontName','Times','interpreter','latex')
title(sprintf('Fit to obtain Jc\nvt01, NE quadrant, die 15, jja\nsouth side devices\nJc = slope$^2$/pi = %g kA/cm$^2$ (slope = %g uA$^{1/2}$/um)\nx-intercept = %g nm',Jc,1e-3*slope,x_intercept*1e9),'FontSize',FontSize,'FontName','Times','interpreter','latex')
set(gca,'FontSize',FontSize,'FontName',FontName);
legend('fit','data')
ylim([0 1.1*max(sqrt_Ic_y_intercept_dense)*1e3])
xlim([0 1.1*radius_vec(end)*1e6])
if strcmp(save_plots,'yes')
    saveas(gcf,'die_15_jja__south__jc_fit','png')
end
if strcmp(close_plots,'yes')
    close
end

%% die 15, south side, hysteresis
data_list = {'vt01_01_NE_die15_jja_south_jj1_02','vt01_01_NE_die15_jja_south_jj2_02','vt01_01_NE_die15_jja_south_jj3_02','vt01_01_NE_die15_jja_south_jj4_02'};
color_list = [3,8,13,18];
title_str = {'d = 1.60um','d = 1.84um','d = 2.26um','d = 3.19um'};
for ii = 1:length(data_list)
    figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
    aa = data_map(data_list{ii});
    plot(aa(:,1)*1e3,aa(:,2)*1e3,'Color',bRGY(color_list(ii),:),'LineWidth',2);
    %     ylim([0 80])
    ylabel('Current [\mu A]','FontSize',FontSize,'FontName','Times')
    xlabel('Voltage [\mu V]','FontSize',FontSize,'FontName','Times')
    set(gca,'FontSize',FontSize,'FontName',FontName);
    title(sprintf(['I-V curve of JJ with ',title_str{ii},'\nforward and reverse sweep\nvt01, NE quadrant, die 15, jja\nsouth side devices']),'FontSize',FontSize,'FontName','Times')
    if strcmp(save_plots,'yes')
        saveas(gcf,sprintf('die_15_jja__south__jjIVs_hysteresis__jj%01d',ii),'png')
    end
    if strcmp(close_plots,'yes')
        close
    end
end

%% die 15, south and east side together
data_list_1 = {'vt01_01_NE_die15_jja_east_jj1_01','vt01_01_NE_die15_jja_east_jj2_01','vt01_01_NE_die15_jja_east_jj3_01','vt01_01_NE_die15_jja_east_jj4_01'};
data_list_2 = {'vt01_01_NE_die15_jja_south_jj1_01','vt01_01_NE_die15_jja_south_jj2_01','vt01_01_NE_die15_jja_south_jj3_01','vt01_01_NE_die15_jja_south_jj4_01'};

color_list_1 = [1,2,3,4];
color_list_2 = [6,7,8,9];

figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
for ii = 1:length(data_list_1)
    aa = data_map(data_list_1{ii});
    plot(aa(:,1)*1e3,aa(:,2)*1e3,'Color',bRGY(color_list_1(ii),:),'LineWidth',2);
    hold on
end
for ii = 1:length(data_list_2)
    bb = data_map(data_list_2{ii});
    plot(bb(:,1)*1e3,bb(:,2)*1e3,'Color',bRGY(color_list_2(ii),:),'LineWidth',2);
    hold on
end
ylabel('Current [\mu A]','FontSize',FontSize,'FontName','Times')
xlabel('Voltage [\mu V]','FontSize',FontSize,'FontName','Times')
set(gca,'FontSize',FontSize,'FontName',FontName);
title(sprintf(['I-V curves of JJs with different shunt resistors\nvt01, NE quadrant, die 15, jja\ncomparing east and south side devices, ',title_string_list{ii}]),'FontSize',FontSize,'FontName','Times')
legend('jj diameter = 1.60um, resistor = 2.05 squares','jj diameter = 1.84um, resistor = 2.05 squares','jj diameter = 2.26um, resistor = 2.05 squares','jj diameter = 3.19um, resistor = 2.05 squares','jj diameter = 1.60um, resistor = 1.15 squares','jj diameter = 1.84um, resistor = 1.15 squares','jj diameter = 2.26um, resistor = 1.15 squares','jj diameter = 3.19um, resistor = 1.15 squares')
if strcmp(save_plots,'yes')
    saveas(gcf,'die_15_jja__south_and_east__jjIVs_all','png')
end
if strcmp(close_plots,'yes')
    close
end

title_string_list = {'jj diameter = 1.60um','jj diameter = 1.84um','jj diameter = 2.26um','jj diameter = 3.19um'};
for ii = 1:length(data_list_1)
    figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
    aa = data_map(data_list_1{ii});
    bb = data_map(data_list_2{ii});
    plot(aa(:,1)*1e3,aa(:,2)*1e3,'Color',bRGY(3,:),'LineWidth',2);
    hold on
    plot(bb(:,1)*1e3,bb(:,2)*1e3,'Color',bRGY(8,:),'LineWidth',2);
    %     ylim([0 80])
    ylabel('Current [\mu A]','FontSize',FontSize,'FontName','Times')
    xlabel('Voltage [\mu V]','FontSize',FontSize,'FontName','Times')
    set(gca,'FontSize',FontSize,'FontName',FontName);
    title(sprintf(['I-V curves of JJs with different shunt resistors\nvt01, NE quadrant, die 15, jja\ncomparing east and south side devices, ',title_string_list{ii}]),'FontSize',FontSize,'FontName','Times')
    legend('resistor = 2.05 squares','resistor = 1.15 squares')
    if strcmp(save_plots,'yes')
        saveas(gcf,sprintf('die_15_jja__south_and_east__jjIVs_jj%01d',ii),'png')
    end
    if strcmp(close_plots,'yes')
        close
    end
end

%% comparing die 15 and die 11, east side
data_list_1 = {'vt01_01_NE_die11_jja_east_jj1_01','vt01_01_NE_die11_jja_east_jj2_01','vt01_01_NE_die11_jja_east_jj3_01','vt01_01_NE_die11_jja_east_jj4_01'};
data_list_2 = {'vt01_01_NE_die15_jja_east_jj1_01','vt01_01_NE_die15_jja_east_jj2_01','vt01_01_NE_die15_jja_east_jj3_01','vt01_01_NE_die15_jja_east_jj4_01'};

color_list_1 = [1,2,3,4];
color_list_2 = [6,7,8,9];
figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
for ii = 1:length(data_list_1)
    aa = data_map(data_list_1{ii});
    plot(aa(:,1)*1e3,aa(:,2)*1e3,'Color',bRGY(color_list_1(ii),:),'LineWidth',2);
    hold on
end
for ii = 1:length(data_list_1)
    bb = data_map(data_list_2{ii});
    plot(bb(:,1)*1e3,bb(:,2)*1e3,'Color',bRGY(color_list_2(ii),:),'LineWidth',2);
    hold on
end
ylabel('Current [\mu A]','FontSize',FontSize,'FontName','Times')
xlabel('Voltage [\mu V]','FontSize',FontSize,'FontName','Times')
set(gca,'FontSize',FontSize,'FontName',FontName);
title('I-V curves of JJs from different chips\nvt01, NE quadrant, jja, east side','FontSize',FontSize,'FontName','Times')
legend('jj diameter = 1.60um, die 11','jj diameter = 1.84um, die 11','jj diameter = 2.26um, die 11','jj diameter = 3.19um, die 11','jj diameter = 1.60um, die 15','jj diameter = 1.84um, die 15','jj diameter = 2.26um, die 15','jj diameter = 3.19um, die 15');
if strcmp(save_plots,'yes')
    saveas(gcf,'die_11_and_15_jja__east__jjIVs','png')
end
if strcmp(close_plots,'yes')
    close
end

title_string_list = {'jj diameter = 1.60um','jj diameter = 1.84um','jj diameter = 2.26um','jj diameter = 3.19um'};
for ii = 1:length(data_list_1)
    figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
    aa = data_map(data_list_1{ii});
    bb = data_map(data_list_2{ii});
    plot(aa(:,1)*1e3,aa(:,2)*1e3,'Color',bRGY(3,:),'LineWidth',2);
    hold on
    plot(bb(:,1)*1e3,bb(:,2)*1e3,'Color',bRGY(8,:),'LineWidth',2);
    %     ylim([0 80])
    ylabel('Current [\mu A]','FontSize',FontSize,'FontName','Times')
    xlabel('Voltage [\mu V]','FontSize',FontSize,'FontName','Times')
    set(gca,'FontSize',FontSize,'FontName',FontName);
    title(sprintf(['I-V curves of JJs from different chips\nvt01, NE quadrant, jja, east side\n',title_string_list{ii}]),'FontSize',FontSize,'FontName','Times')
    legend('die 11','die 15')
    if strcmp(save_plots,'yes')
        saveas(gcf,sprintf('die_11_and_15_jja__east__jjIVs_jj%01d',ii),'png')
    end
    if strcmp(close_plots,'yes')
        close
    end
end

%% comparing die 15 and die 11, south side
data_list_1 = {'vt01_01_NE_die11_jja_south_jj1_01','vt01_01_NE_die11_jja_south_jj2_01','vt01_01_NE_die11_jja_south_jj3_01','vt01_01_NE_die11_jja_south_jj4_01'};
data_list_2 = {'vt01_01_NE_die15_jja_south_jj1_01','vt01_01_NE_die15_jja_south_jj2_01','vt01_01_NE_die15_jja_south_jj3_01','vt01_01_NE_die15_jja_south_jj4_01'};

color_list_1 = [1,2,3,4];
color_list_2 = [6,7,8,9];
figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
for ii = 1:length(data_list_1)
    aa = data_map(data_list_1{ii});
    plot(aa(:,1)*1e3,aa(:,2)*1e3,'Color',bRGY(color_list_1(ii),:),'LineWidth',2);
    hold on
end
for ii = 1:length(data_list_1)
    bb = data_map(data_list_2{ii});
    plot(bb(:,1)*1e3,bb(:,2)*1e3,'Color',bRGY(color_list_2(ii),:),'LineWidth',2);
    hold on
end
ylabel('Current [\mu A]','FontSize',FontSize,'FontName','Times')
xlabel('Voltage [\mu V]','FontSize',FontSize,'FontName','Times')
set(gca,'FontSize',FontSize,'FontName',FontName);
title('I-V curves of JJs from different chips\nvt01, NE quadrant, jja, south side','FontSize',FontSize,'FontName','Times')
legend('jj diameter = 1.60um, die 11','jj diameter = 1.84um, die 11','jj diameter = 2.26um, die 11','jj diameter = 3.19um, die 11','jj diameter = 1.60um, die 15','jj diameter = 1.84um, die 15','jj diameter = 2.26um, die 15','jj diameter = 3.19um, die 15');
if strcmp(save_plots,'yes')
    saveas(gcf,'die_11_and_15_jja__south__jjIVs','png')
end
if strcmp(close_plots,'yes')
    close
end

title_string_list = {'jj diameter = 1.60um','jj diameter = 1.84um','jj diameter = 2.26um','jj diameter = 3.19um'};
for ii = 1:length(data_list_1)
    figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
    aa = data_map(data_list_1{ii});
    bb = data_map(data_list_2{ii});
    plot(aa(:,1)*1e3,aa(:,2)*1e3,'Color',bRGY(3,:),'LineWidth',2);
    hold on
    plot(bb(:,1)*1e3,bb(:,2)*1e3,'Color',bRGY(8,:),'LineWidth',2);
    %     ylim([0 80])
    ylabel('Current [\mu A]','FontSize',FontSize,'FontName','Times')
    xlabel('Voltage [\mu V]','FontSize',FontSize,'FontName','Times')
    set(gca,'FontSize',FontSize,'FontName',FontName);
    title(sprintf(['I-V curves of JJs from different chips\nvt01, NE quadrant, jja, south side\n',title_string_list{ii}]),'FontSize',FontSize,'FontName','Times')
    legend('die 11','die 15')
    if strcmp(save_plots,'yes')
        saveas(gcf,sprintf('die_11_and_15_jja__south__jjIVs_jj%01d',ii),'png')
    end
    if strcmp(close_plots,'yes')
        close
    end
end

%% all jj 4wires from both die
data_list = {'vt01_01_NE_die11_jja_south_jj1_01','vt01_01_NE_die11_jja_south_jj2_01','vt01_01_NE_die11_jja_south_jj3_01','vt01_01_NE_die11_jja_south_jj4_01',...             
             'vt01_01_NE_die15_jja_south_jj1_01','vt01_01_NE_die15_jja_south_jj2_01','vt01_01_NE_die15_jja_south_jj3_01','vt01_01_NE_die15_jja_south_jj4_01',...
             'vt01_01_NE_die11_jja_east_jj1_01','vt01_01_NE_die11_jja_east_jj2_01','vt01_01_NE_die11_jja_east_jj3_01','vt01_01_NE_die11_jja_east_jj4_01',...
             'vt01_01_NE_die15_jja_east_jj1_01','vt01_01_NE_die15_jja_east_jj2_01','vt01_01_NE_die15_jja_east_jj3_01','vt01_01_NE_die15_jja_east_jj4_01'};

color_list = [1,2,3,4,6,7,8,9,11,12,13,14,16,17,18,19];
figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
for ii = 1:length(data_list)
    aa = data_map(data_list{ii});
    plot(aa(:,1)*1e3,aa(:,2)*1e3,'Color',bRGY(color_list(ii),:),'LineWidth',2);
    hold on
end
ylabel('Current [\mu A]','FontSize',FontSize,'FontName','Times')
xlabel('Voltage [\mu V]','FontSize',FontSize,'FontName','Times')
set(gca,'FontSize',FontSize,'FontName',FontName);
title(sprintf('I-V curves of JJs \nvt01, NE quadrant, jja\nall jj 4wire'),'FontSize',FontSize,'FontName','Times')
legend('die 11, south, jj1','die 11, south, jj2','die 11, south, jj3','die 11, south, jj4',...
       'die 15, south, jj1','die 15, south, jj2','die 15, south, jj3','die 15, south, jj4',...
       'die 11, east, jj1','die 11, east, jj2','die 11, east, jj3','die 11, east, jj4',...
       'die 15, east, jj1','die 15, east, jj2','die 15, east, jj3','die 15, east, jj4');
if strcmp(save_plots,'yes')
    saveas(gcf,'die_11_and_15_jja__all_jj4wireIVs','png')
end
if strcmp(close_plots,'yes')
    close
end
