%% initialize
clc
clear all
% close all

[fontName,fontSize,fontSize_legend,bRGY,scrsz] = f_plotting;
p = f_physicalConstants;

%%

Jc_vec = [0.5e3 1e3 2e3 4e3]*1e4;%critical current density in A/m^2
d_vec = linspace(0.1e-6,3.5e-6,100);%junction width in meters
a_vec = (pi/4)*d_vec.^2;%junction area in meters for square junction
target_Ic = 40e-6;

Ic_vec = zeros(length(Jc_vec),length(d_vec));
target_vec = zeros(length(Jc_vec),1);
for ii = 1:length(Jc_vec)
    Ic_vec(ii,:) = Jc_vec(ii)*a_vec(:);
    [~,ind] = min( abs(Ic_vec(ii,:)-target_Ic) );
    target_vec(ii) = ind;
end



%%
figureCaptions = {sprintf('For Ic = %2.1f uA and J_{c} = %2.1f kA/cm^2, d = %1.2f um',target_Ic*1e6,Jc_vec(1)*1e-7,d_vec(target_vec(1))*1e6),...
    sprintf('For Ic = %2.1f uA and J_{c} = %2.1f kA/cm^2, d = %1.2f um',target_Ic*1e6,Jc_vec(2)*1e-7,d_vec(target_vec(2))*1e6),...
    sprintf('For Ic = %2.1f uA and J_{c} = %2.1f kA/cm^2, d = %1.2f um',target_Ic*1e6,Jc_vec(3)*1e-7,d_vec(target_vec(3))*1e6),...
    sprintf('For Ic = %2.1f uA and J_{c} = %2.1f kA/cm^2, d = %1.2f um',target_Ic*1e6,Jc_vec(4)*1e-7,d_vec(target_vec(4))*1e6),...
    };

legend_str = 'legend(';
color_vec = [3 8 13 18];
figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
for ii = 1:length(Jc_vec)
    plot(d_vec*1e6,Ic_vec(ii,:)*1e6,'Color',bRGY(color_vec(ii),:),'LineStyle','-','LineWidth',3)
    hold on
    legend_str = [legend_str sprintf('''Jc = %1.1f kA/cm^2'',',Jc_vec(ii)*1e-7)];
end
line([d_vec(1) d_vec(end)]*1e6,[target_Ic target_Ic]*1e6,'LineStyle','-.','LineWidth',1,'Color',bRGY(end,:))
for ii = 1:length(Jc_vec)
    line([d_vec(target_vec(ii)) d_vec(target_vec(ii))]*1e6,[min(Ic_vec(ii,:)) max(Ic_vec(ii,:))]*1e6,'LineStyle','-.','LineWidth',1,'Color',bRGY(color_vec(ii),:))
end
legend_str = [legend_str(1:end-1) ');'];
eval(legend_str)
xlabel('junction diameter [um]','FontSize',fontSize,'FontName',fontName)
ylabel('I_c [uA]','FontSize',fontSize,'FontName',fontName)
set(gca,'FontSize',fontSize,'FontName',fontName)
xlim([d_vec(1) d_vec(end)]*1e6)
ylim([0 1.1*target_Ic]*1e6)
k1 = gtext(figureCaptions(1:length(figureCaptions)));
set(k1,'FontSize',fontSize_legend,'FontName',fontName)