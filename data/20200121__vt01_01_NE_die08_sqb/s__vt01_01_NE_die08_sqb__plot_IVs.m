%%
clc;
close all;
% clear all
[~,bRGY] = f_colorSchemes('redWhiteBlue');
FontSize = 20;
FontName = 'Times';
scrsz = get(0,'ScreenSize');

%% switches
load_data = 'yes';
save_plots = 'yes';
close_plots = 'no';
run_thermal_fits = 'no';

%% load all
if strcmp(load_data,'yes')
    
%     directory_list = {'20200121__vt01_01_NE_die08_sqb'};
    
    current_dir = pwd;
    data_map = containers.Map();
%     for ii = 1:length(directory_list)
%         cd(directory_list{ii})
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
%     end
    
end

%% die 8, east side, squid 3
% data_list = {'20200121__vt01_01_ne_die08__device3__sweep_I_flux_04','20200121__vt01_01_ne_die08__device3__sweep_I_flux_05'};
data_list = {'20200121__vt01_01_ne_die08__device3__sweep_I_flux_04__zoom'};

%voltage versus flux
color_list = [3,8,13,18];
figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
legend_str = 'legend(';
for ii = 100:160
    legend_str = [legend_str sprintf('''I_bias = %g uA'',',ii)];
end
legend_str = [legend_str(1:end-1) ');'];
for ii = 1:length(data_list)
    aa = data_map(data_list{ii});
    for jj = 1:length(aa(:,1))
        aa(jj,2) = aa(jj,2)-jj*0.0001;
    end
    plot(aa(:,1),aa(:,2),'Color',bRGY(3,:),'LineWidth',2);
    hold on
end
% ylim([0 80])
xlabel('Applied flux [mA]','FontSize',FontSize,'FontName','Times')
ylabel('SQUID voltage [mV]','FontSize',FontSize,'FontName','Times')
set(gca,'FontSize',FontSize,'FontName',FontName);
title(sprintf('Voltage versus applied flux for several values of current bias\nvt01, NE quadrant, die 8, sqb\neast side, device 3'),'FontSize',FontSize,'FontName','Times')
% legend('d = 1.60\mum','d = 1.84\mum','d = 2.26\mum','d = 3.19\mum')
if strcmp(save_plots,'yes')
    saveas(gcf,'die08_sqb__east_device3__voltage_v_flux','png')
end
if strcmp(close_plots,'yes')
    close
end

% %differential resistance
% figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
% for ii = 1:length(data_list)
%     aa = data_map(data_list{ii});
%     bb = diff(aa(:,1))./diff(aa(:,2));
%     plot(aa(1:end-1,1),bb(:),'Color',bRGY(color_list(ii),:),'LineWidth',2);
%     hold on
% end
% ylim([5 20])
% xlabel('Voltage [\mu V]','FontSize',FontSize,'FontName','Times')
% ylabel('Differential Resistance [\Delta V/\Delta I, \Omega]','FontSize',FontSize,'FontName','Times')
% set(gca,'FontSize',FontSize,'FontName',FontName);
% title(sprintf('Differential resistance of JJs with different diameters\nvt01, NE quadrant, die 11, jja\neast side devices'),'FontSize',FontSize,'FontName','Times')
% legend('d = 1.60\mum','d = 1.84\mum','d = 2.26\mum','d = 3.19\mum')
% if strcmp(save_plots,'yes')
%     saveas(gcf,'die_11_jja__east__diffRes','png')
% end
% if strcmp(close_plots,'yes')
%     close
% end
% 
% %thermal fits
% if strcmp(run_thermal_fits,'yes')
%     Ic_guesses = [14e-6 22e-6 30e-6 60e-6];
%     T_guesses = [10 10 10 10];
%     Ic_fit = zeros(length(Ic_guesses),1);
%     R_fit = zeros(length(Ic_guesses),1);
%     T_fit = zeros(length(Ic_guesses),1);
%     for kk = 1:length(data_list)
%         
%         aa = data_map(data_list{kk});
%         vV = aa(:,1);
%         iI = aa(:,2);
%         
%         %     ind = find( abs(iI) < 7e-3 );
%         %     P = polyfit(iI(ind),vV(ind),1);
%         %     vV = vV-polyval(P,iI);
%         %     ind = find( vV > 0.01 );
%         
%         aux = diff(vV)./diff(iI);
%         R_guess = mean(aux(end-10:end));
%         
%         ft = fittype( 'IVthermal(x,Ic,R,T)', 'independent', 'x','coefficient', {'Ic','R','T'});
%         opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
%         opts.Display = 'Off';
%         opts.StartPoint = [Ic_guesses(kk)*1e6,R_guess,T_guesses(kk)];
%         opts.Lower = [0.8*Ic_guesses(kk)*1e6  0.1*R_guess  0.01*T_guesses(kk)];
%         opts.Upper = [1.2*Ic_guesses(kk)*1e6 10*R_guess 10*T_guesses(kk)];
%         
%         %     ind = find( iI < 1.5*Ic_guesses(kk)*1e3 );
%         %fit wants things close to 1, so scale the data to uA and uV
%         [xfit,yfit] = prepareCurveData(iI(:)*1e3,vV(:)*1e3);
%         [fitresult, gof] = fit(xfit,yfit,ft,opts);
%         Ic_fit(kk) = fitresult.Ic;
%         R_fit(kk) = fitresult.R;
%         T_fit(kk) = fitresult.T;
%         %     plot(yfit,xfit,'.',fitresult(xfit),xfit,'linewidth',2)
%         figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
%         plot(vV*1e3,iI*1e3,'Color',bRGY(3,:),'LineWidth',2)
%         hold on
%         plot(IVthermal(iI*1e3,fitresult.Ic,fitresult.R,fitresult.T)',iI*1e3,'Color',bRGY(8,:),'LineWidth',2)
%         xlabel('Voltage [\mu V]','FontSize',FontSize,'FontName','Times')
%         ylabel('Current [\mu A]','FontSize',FontSize,'FontName','Times')
%         set(gca,'FontSize',FontSize,'FontName',FontName);
%         title(sprintf('Thermal fit of JJ I-V\nvt01, NE quadrant, die 11, jja\neast side devices, jj%01d\nIc = %g uA, R = %g Ohm, T = %g K',kk,fitresult.Ic,fitresult.R,fitresult.T),'FontSize',FontSize,'FontName','Times')
%         legend('Data','Fit')
%         if strcmp(save_plots,'yes')
%             saveas(gcf,sprintf('die_11_jja__east__thermal_fits__jj%01d',kk),'png')
%         end
%         if strcmp(close_plots,'yes')
%             close
%         end
%     end
% end
% 
% %Jc
% diameter_vec = [1.6 1.84 2.26 3.19]*1e-6;
% radius_vec = diameter_vec/2;
% area_vec = (pi/4)*diameter_vec.^2;
% Ic_vec = [13.7 20.5 29.2 59]*1e-6;%Ic_fit*1e-6;
% Jc_vec = Ic_vec./area_vec;
% figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
% plot(area_vec*1e12,Jc_vec*1e-7,'Color',bRGY(color_list(ii),:),'LineWidth',2);
% xlabel('Junction area [\mum^2]','FontSize',FontSize,'FontName','Times')
% ylabel('Critical current density [kA/cm^2]','FontSize',FontSize,'FontName','Times')
% title(sprintf('vt01, NE quadrant, die 11, jja\neast side devices'),'FontSize',FontSize,'FontName','Times')
% set(gca,'FontSize',FontSize,'FontName',FontName);
% if strcmp(save_plots,'yes')
%     saveas(gcf,'die_11_jja__east__jc','png')
% end
% if strcmp(close_plots,'yes')
%     close
% end
% 
% p = polyfit(radius_vec,sqrt(Ic_vec),1);
% % radius_vec_dense = linspace(radius_vec(1),radius_vec(end),100);
% radius_vec_dense = linspace(0,radius_vec(end),100);
% sqrt_Ic_vec_dense = polyval(p,radius_vec_dense);
% slope = p(1);
% Jc = 1e-7*slope^2/pi;
% x_intercept = -p(2)/slope;
% figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
% plot(radius_vec_dense*1e6,sqrt_Ic_vec_dense*1e3,'Color',bRGY(3,:),'LineWidth',2);
% hold on
% plot(radius_vec*1e6,sqrt(Ic_vec)*1e3,'Color',bRGY(8,:),'LineWidth',2,'Marker','s','MarkerFaceColor',bRGY(6,:),'MarkerEdgeColor',bRGY(10,:),'MarkerSize',8);
% % line([x_intercept x_intercept]*1e6,[min(sqrt_Ic_y_intercept_dense) max(sqrt_Ic_y_intercept_dense)]*1e3,'LineWidth',1,'LineStyle','-.','Color',bRGY(21,:))
% xlabel('Junction radius [\mum]','FontSize',FontSize,'FontName','Times')
% ylabel('$\sqrt{I_c}$ [$\mu$A$^{1/2}$]','FontSize',FontSize,'FontName','Times','interpreter','latex')
% title(sprintf('Fit to obtain Jc\nvt01, NE quadrant, die 11, jja\neast side devices\nJc = slope$^2$/pi = %g kA/cm$^2$ (slope = %g uA$^{1/2}$/um)\nx-intercept = %g nm',Jc,1e-3*slope,x_intercept_alt*1e9),'FontSize',FontSize,'FontName','Times','interpreter','latex')
% set(gca,'FontSize',FontSize,'FontName',FontName);
% legend('fit','data')
% ylim([0 1.1*max(sqrt_Ic_y_intercept_dense)*1e3])
% xlim([0 1.1*radius_vec(end)*1e6])
% if strcmp(save_plots,'yes')
%     saveas(gcf,'die_11_jja__east__jc_fit','png')
% end
% if strcmp(close_plots,'yes')
%     close
% end

