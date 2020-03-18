%% initialize
clc; close all;
% clear all

[fontName,fontSize,fontSize_legend,bRGY,scrsz] = f_plotting;

p = f_physicalConstants;

%% input parameters
%from WRSpice
p.fileNames = {'_dat__squid__Ic_20uA__Ib_36uA__L1_500pH__L2_25pH__beta_c_0p3',...
               '_dat__squid__Ic_20uA__Ib_38uA__L1_500pH__L2_25pH__beta_c_0p3',...
               '_dat__squid__Ic_20uA__Ib_40uA__L1_500pH__L2_25pH__beta_c_0p3',...
               '_dat__squid__Ic_20uA__Ib_42uA__L1_500pH__L2_25pH__beta_c_0p3',...
               '_dat__squid__Ic_20uA__Ib_44uA__L1_500pH__L2_25pH__beta_c_0p3',...
               };

%for jjs
I_c = 20e-6;
beta_c = 0.3;
I_b_vec = [36 38 40 42 44]*1e-6;

%for mutual inductors
k = 0.5;
L1 = 500e-12;
L2 = 25e-12;
L_parasitic = 2*1e-12;%factor of two for WR symmetry
mutual_inductance = 2*k*sqrt(L1*L2);%factor of two because WRSpice simulation has two MIs for symmetry

%% switches
s.loadData = 'no';

%% read and format data
if strcmp(s.loadData,'yes')
    
    dataMat = cell(length(p.fileNames),1);
    for kk = 1:length(p.fileNames)
        
        fprintf('\nReading file %d of %d\n',kk,length(p.fileNames))
        
        fileID = fopen(p.fileNames{kk},'r');
        C = textscan(fileID,'%s');
        A = C{1};
        numRows = length(A);
        tN = find(strcmp(A,'Variables:'));
        tN2 = tN(2);
        numVars = str2double(A{tN(1)+1});
        
        tN = find(strcmp(A,'Points:'));
        numPts = str2double(A{tN(1)+1});
        
        varList = cell(numVars,1);
        for ii = 1:numVars
            varList{ii} = A{tN2+2+3*(ii-1)};
        end
        
        tN = find(strcmp(A,'Values:'));
        dataMat{kk} = zeros(numPts,numVars);
        tN = tN+2;
        for ii = 1:numPts
            for jj = 1:numVars
                dataMat{kk}(ii,jj) = str2double(A{tN+jj-1});
            end
            tN = tN+1+numVars;
        end
        fclose(fileID);
        
    end
    
end

%% further processing

num_files = length(dataMat);
time_vec = cell(num_files,1);
current_vec = cell(num_files,1);
voltage_vec = cell(num_files,1);
avg_time_vec = cell(num_files,1);
avg_current_vec = cell(num_files,1);
avg_voltage_vec = cell(num_files,1);
for ii = 1:num_files
    
    time_vec{ii} = dataMat{ii}(:,1);
    current_vec{ii} = dataMat{ii}(:,2)/2;%dividing by 2 because current is split to two MIs for symmetry
    voltage_vec{ii} = dataMat{ii}(:,8);
    
%     norm_current_vec = current_vec/max(current_vec);
%     norm_voltage_vec = voltage_vec/max(voltage_vec);
    
    %perform averaging
    time_bin_duration = 1e-9;
    num_t_step = length(time_vec{ii});
    % temp_index = 1;
    start_time = 1e-9;
    [~,ind_1] = min( abs( time_vec{ii} - start_time ) );
    ind_main = 1;
    while ind_1 < num_t_step
        fprintf('ind_main = %g\n',ind_main)
        fprintf('ind_1 = %g of num_t_step = %g\n\n',ind_1,num_t_step)
        [~,ind_2] = min( abs( time_vec{ii} - ( time_vec{ii}(ind_1)+time_bin_duration ) ) );
        temp_vec_1 = current_vec{ii}(ind_1:ind_2);
        temp_vec_2 = ( time_vec{ii}(ind_1:ind_2) - time_vec{ii}((ind_1-1):(ind_2-1)) );
        avg_current_vec{ii}(ind_main) = (temp_vec_1'*temp_vec_2)/(time_vec{ii}(ind_2)-time_vec{ii}(ind_1));
        temp_vec_1 = voltage_vec{ii}(ind_1:ind_2);
        avg_voltage_vec{ii}(ind_main) = (temp_vec_1'*temp_vec_2)/(time_vec{ii}(ind_2)-time_vec{ii}(ind_1));
        avg_time_vec{ii}(ind_main) = time_vec{ii}(ind_1) + ( time_vec{ii}(ind_2) - time_vec{ii}(ind_1) )/2;
        
        ind_main = ind_main+1;
        ind_1 = ind_2;
    end
    
    % plots
    
    % figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
    % plot(time_vec*1e9,norm_voltage_vec,'Color',bRGY(8,:),'LineStyle','-','LineWidth',2)
    % hold on
    % plot(time_vec*1e9,norm_current_vec,'Color',bRGY(3,:),'LineStyle','-','LineWidth',2)
    % xlabel('Time [ns]','FontSize',fontSize,'FontName','Times')
    % ylabel('Signal [Arb.]','FontSize',fontSize,'FontName','Times')
    % set(gca,'FontSize',fontSize,'FontName',fontName)
    % % ylim([0 11])
    % % xlim([0 160])
    % legend('SQUID Voltage','Applied Current')
    % grid on
    %
    % figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
    % plot(avg_time_vec*1e9,avg_current_vec*1e6,'Color',bRGY(8,:),'LineStyle','-','LineWidth',2)
    % xlabel('Time [ns]','FontSize',fontSize,'FontName','Times')
    % ylabel('Applied Current [\mu A]','FontSize',fontSize,'FontName','Times')
    % set(gca,'FontSize',fontSize,'FontName',fontName)
    % grid on
    %
    % figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
    % plot(avg_time_vec*1e9,avg_voltage_vec*1e6,'Color',bRGY(8,:),'LineStyle','-','LineWidth',2)
    % xlabel('Time [ns]','FontSize',fontSize,'FontName','Times')
    % ylabel('SQUID voltage [\mu V]','FontSize',fontSize,'FontName','Times')
    % set(gca,'FontSize',fontSize,'FontName',fontName)
    % grid on
    
%     figureCaptions = {sprintf('JJ I_c = %g uA',I_c*1e6),...
%         sprintf('JJ beta_c = %g',beta_c),...
%         sprintf('bias to SQUID I_b = %g uA',I_b_vec(ii)*1e6),...
%         sprintf('I_b/2I_c = %g',I_b_vec(ii)/(2*I_c)),...
%         sprintf('L_1 = %g pH',L1*1e12),...
%         sprintf('L_2 = %g pH',L2*1e12),...
%         sprintf('k = %g',k),...
%         sprintf('M = 2k(L1 L2)^{1/2} = %g pH',mutual_inductance*1e12),...
%         sprintf('L_{parasitic} = %g pH',L_parasitic*1e12),...
%         sprintf('SQUID beta_L = %g',2*I_c*(2*L2+L_parasitic)/p.Phi0),...
%         sprintf('factor of 2 in M comes from two MIs for symmetry'),...
%         };
%     
%     figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
%     plot(avg_current_vec{ii}*1e6,avg_voltage_vec{ii}*1e6,'Color',bRGY(3,:),'LineStyle','-','LineWidth',2)
%     xlabel('Applied Current [\mu A]','FontSize',fontSize,'FontName','Times')
%     ylabel('Time-averaged SQUID Voltage [\mu V]','FontSize',fontSize,'FontName','Times')
%     set(gca,'FontSize',fontSize,'FontName',fontName)
%     grid on
%     k1 = gtext(figureCaptions(1:length(figureCaptions)));
%     set(k1,'FontSize',fontSize_legend,'FontName','Times')
%     
%     figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
%     plot(mutual_inductance*avg_current_vec{ii}/p.Phi0,avg_voltage_vec{ii}*1e6,'Color',bRGY(3,:),'LineStyle','-','LineWidth',2)
%     xlabel('Applied Flux [\Phi/\Phi_0]','FontSize',fontSize,'FontName','Times')
%     ylabel('Time-averaged SQUID Voltage [\mu V]','FontSize',fontSize,'FontName','Times')
%     set(gca,'FontSize',fontSize,'FontName',fontName)
%     grid on
%     k1 = gtext(figureCaptions(1:length(figureCaptions)));
%     set(k1,'FontSize',fontSize_legend,'FontName','Times')
    
end

%% figure comparing Ib values
figureCaptions = {sprintf('JJ I_c = %g uA',I_c*1e6),...
    sprintf('JJ beta_c = %g',beta_c),...
    sprintf('L_1 = %g pH',L1*1e12),...
    sprintf('L_2 = %g pH',L2*1e12),...
    sprintf('k = %g',k),...
    sprintf('M = 2k(L1 L2)^{1/2} = %g pH',mutual_inductance*1e12),...
    sprintf('L_{parasitic} = %g pH',L_parasitic*1e12),...
    sprintf('SQUID beta_L = %g',2*I_c*(2*L2+L_parasitic)/p.Phi0),...
    sprintf('factor of 2 in M comes from two MIs for symmetry'),...
    };
    
figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
plot(mutual_inductance*avg_current_vec{1}/p.Phi0,avg_voltage_vec{1}*1e6,'Color',bRGY(15,:),'LineStyle','-','LineWidth',2)
hold on
plot(mutual_inductance*avg_current_vec{2}/p.Phi0,avg_voltage_vec{2}*1e6,'Color',bRGY(14,:),'LineStyle','-','LineWidth',2)
plot(mutual_inductance*avg_current_vec{3}/p.Phi0,avg_voltage_vec{3}*1e6,'Color',bRGY(13,:),'LineStyle','-','LineWidth',2)
plot(mutual_inductance*avg_current_vec{4}/p.Phi0,avg_voltage_vec{4}*1e6,'Color',bRGY(12,:),'LineStyle','-','LineWidth',2)
plot(mutual_inductance*avg_current_vec{5}/p.Phi0,avg_voltage_vec{5}*1e6,'Color',bRGY(11,:),'LineStyle','-','LineWidth',2)
xlabel('Applied Flux [\Phi/\Phi_0]','FontSize',fontSize,'FontName','Times')
ylabel('Time-averaged SQUID Voltage [\mu V]','FontSize',fontSize,'FontName','Times')
legend('I_b = 36 \mu A','I_b = 38 \mu A','I_b = 40 \mu A','I_b = 42 \mu A','I_b = 44 \mu A')
title('Comparing SQUID response for several values of I_b')
set(gca,'FontSize',fontSize,'FontName',fontName)
grid on
k1 = gtext(figureCaptions(1:length(figureCaptions)));
set(k1,'FontSize',fontSize_legend,'FontName','Times')
xlim([0 1])
