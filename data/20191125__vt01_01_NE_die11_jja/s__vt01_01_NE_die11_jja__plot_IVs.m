%%
clc;
close all;
% clear all
[~,bRGY] = f_colorSchemes('redWhiteBlue');
FontSize = 20;
FontName = 'Times';
scrsz = get(0,'ScreenSize');

%% east side
fig_name_list = {'vt01_01_NE_die11_jja_east_jj1_01','vt01_01_NE_die11_jja_east_jj2_01','vt01_01_NE_die11_jja_east_jj3_01','vt01_01_NE_die11_jja_east_jj4_01'};
x_data = cell(length(fig_name_list),1);
y_data = cell(length(fig_name_list),1);
for ii = 1:length(fig_name_list)
    open([fig_name_list{ii} '.fig'])
    [x_data{ii},y_data{ii}] = ax_data(gca,1);
% [x_data,y_data] = ax_data(gca,1);
    close
%     figure;
%     plot(x_data,y_data)
%     pause(5)
end

color_list = [3,8,13,18];
h = figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
for ii = 1:length(fig_name_list)
    plot(x_data{ii}*1e3,y_data{ii}*1e3,'Color',bRGY(color_list(ii),:),'LineWidth',2);
    hold on
end
ylim([0 80])
ylabel('Current [\mu A]','FontSize',FontSize,'FontName','Times')
xlabel('Voltage [\mu V]','FontSize',FontSize,'FontName','Times')
set(gca,'FontSize',FontSize,'FontName',FontName);
title(sprintf('I-V curves of JJs with different diameters\nvt01, NE quadrant, die 11, jja\neast side devices'),'FontSize',FontSize,'FontName','Times')
legend('d = 1.60\mum','d = 1.84\mum','d = 2.26\mum','d = 3.19\mum')

diameter_vec = [1.6 1.84 2.26 3.19]*1e-6;
area_vec = (pi/4)*diameter_vec.^2;
Ic_vec = [14 20 29 59]*1e-6;
Jc_vec = Ic_vec./area_vec;
h = figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
plot(area_vec*1e12,Jc_vec*1e-7,'Color',bRGY(color_list(ii),:),'LineWidth',2);
ylabel('Critical current density [kA/cm^2]','FontSize',FontSize,'FontName','Times')
xlabel('Junction area [\mum^2]','FontSize',FontSize,'FontName','Times')
title(sprintf('vt01, NE quadrant, die 11, jja\neast side devices'),'FontSize',FontSize,'FontName','Times')
set(gca,'FontSize',FontSize,'FontName',FontName);

%% south side
fig_name_list = {'vt01_01_NE_die11_jja_south_jj1_01','vt01_01_NE_die11_jja_south_jj2_01','vt01_01_NE_die11_jja_south_jj3_01','vt01_01_NE_die11_jja_south_jj4_01'};
x_data = cell(length(fig_name_list),1);
y_data = cell(length(fig_name_list),1);
for ii = 1:length(fig_name_list)
    open([fig_name_list{ii} '.fig'])
    [x_data{ii},y_data{ii}] = ax_data(gca,1);
% [x_data,y_data] = ax_data(gca,1);
    close
%     figure;
%     plot(x_data,y_data)
%     pause(5)
end

color_list = [3,8,13,18];
h = figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
for ii = 1:length(fig_name_list)
    plot(x_data{ii}*1e3,y_data{ii}*1e3,'Color',bRGY(color_list(ii),:),'LineWidth',2);
    hold on
end
ylim([0 80])
ylabel('Current [\mu A]','FontSize',FontSize,'FontName','Times')
xlabel('Voltage [\mu V]','FontSize',FontSize,'FontName','Times')
set(gca,'FontSize',FontSize,'FontName',FontName);
title(sprintf('I-V curves of JJs with different diameters\nvt01, NE quadrant, die 11, jja\nsouth side devices'),'FontSize',FontSize,'FontName','Times')
legend('d = 1.60\mum','d = 1.84\mum','d = 2.26\mum','d = 3.19\mum')

diameter_vec = [1.6 1.84 2.26 3.19]*1e-6;
area_vec = (pi/4)*diameter_vec.^2;
Ic_vec = [11 17 29 59]*1e-6;
Jc_vec = Ic_vec./area_vec;
h = figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
plot(area_vec*1e12,Jc_vec*1e-7,'Color',bRGY(color_list(ii),:),'LineWidth',2);
ylabel('Critical current density [kA/cm^2]','FontSize',FontSize,'FontName','Times')
xlabel('Junction area [\mum^2]','FontSize',FontSize,'FontName','Times')
title(sprintf('vt01, NE quadrant, die 11, jja\nsouth side devices'),'FontSize',FontSize,'FontName','Times')
set(gca,'FontSize',FontSize,'FontName',FontName);

%% south and east side together
fig_name_list_1 = {'vt01_01_NE_die11_jja_east_jj1_01','vt01_01_NE_die11_jja_east_jj2_01','vt01_01_NE_die11_jja_east_jj3_01','vt01_01_NE_die11_jja_east_jj4_01'};
fig_name_list_2 = {'vt01_01_NE_die11_jja_south_jj1_01','vt01_01_NE_die11_jja_south_jj2_01','vt01_01_NE_die11_jja_south_jj3_01','vt01_01_NE_die11_jja_south_jj4_01'};

x_data_1 = cell(length(fig_name_list_1),1);
y_data_1 = cell(length(fig_name_list_1),1);
x_data_2 = cell(length(fig_name_list_2),1);
y_data_2 = cell(length(fig_name_list_2),1);
for ii = 1:length(fig_name_list)
    open([fig_name_list_1{ii} '.fig'])
    [x_data_1{ii},y_data_1{ii}] = ax_data(gca,1);
    close
    open([fig_name_list_2{ii} '.fig'])
    [x_data_2{ii},y_data_2{ii}] = ax_data(gca,1);
    close
end

title_string_list = {'jj diameter = 1.60um','jj diameter = 1.84um','jj diameter = 2.26um','jj diameter = 3.19um'};
for ii = 1:length(fig_name_list_1)
    h = figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
    plot(x_data_1{ii}*1e3,y_data_1{ii}*1e3,'Color',bRGY(3,:),'LineWidth',2);
    hold on
    plot(x_data_2{ii}*1e3,y_data_2{ii}*1e3,'Color',bRGY(8,:),'LineWidth',2);
%     ylim([0 80])
    ylabel('Current [\mu A]','FontSize',FontSize,'FontName','Times')
    xlabel('Voltage [\mu V]','FontSize',FontSize,'FontName','Times')
    set(gca,'FontSize',FontSize,'FontName',FontName);
    title(sprintf(['I-V curves of JJs with different shunt resistors\nvt01, NE quadrant, die 11, jja\ncomparing east and south side devices, ',title_string_list{ii}]),'FontSize',FontSize,'FontName','Times')
    legend('resistor = 2.05 squares','resistor = 1.15 squares')
end
