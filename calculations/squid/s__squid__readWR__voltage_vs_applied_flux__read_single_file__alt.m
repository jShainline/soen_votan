%% initialize
clc; close all;
% clear all

[fontName,fontSize,fontSize_legend,bRGY,scrsz] = f_plotting;

p = f_physicalConstants;

%% input parameters
%from WRSpice
p.fileNames = {'_dat__squid__Ic_20uA__Ib_40uA__L1_500pH__L2_25pH__beta_c_0p3__alt_layout'};

%for jjs
I_c = 20e-6;
beta_c = 0.3;
I_b = 40e-6;

%for mutual inductors
k = 0.5;
L1 = 500e-12;
L2 = 25e-12;
mutual_inductance = 2*k*sqrt(L1*L2);%factor of two because WRSpice simulation has two MIs for symmetry

time_bin_duration = 2e-9;
start_time = 1e-9;

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

%% further formatting

time_vec = dataMat{1}(:,1);
current_vec = dataMat{1}(:,2)/2;%dividing by 2 because current is split to two MIs for symmetry
voltage_vec = dataMat{1}(:,3);

norm_current_vec = current_vec/max(current_vec);
norm_voltage_vec = voltage_vec/max(voltage_vec);

%perform averaging
num_t_step = length(time_vec);
% temp_index = 1;
[~,ind_1] = min( abs( time_vec - start_time ) );
ind_main = 1;
while ind_1 < num_t_step
    fprintf('ind_main = %g\n',ind_main)
    fprintf('ind_1 = %g of num_t_step = %g\n\n',ind_1,num_t_step)
    [~,ind_2] = min( abs( time_vec - ( time_vec(ind_1)+time_bin_duration ) ) );
    temp_vec_1 = current_vec(ind_1:ind_2);
    temp_vec_2 = ( time_vec(ind_1:ind_2) - time_vec((ind_1-1):(ind_2-1)) );
    avg_current_vec(ind_main) = (temp_vec_1'*temp_vec_2)/(time_vec(ind_2)-time_vec(ind_1));
    temp_vec_1 = voltage_vec(ind_1:ind_2);
    avg_voltage_vec(ind_main) = (temp_vec_1'*temp_vec_2)/(time_vec(ind_2)-time_vec(ind_1));    
    avg_time_vec(ind_main) = time_vec(ind_1) + ( time_vec(ind_2) - time_vec(ind_1) )/2;
    
    ind_main = ind_main+1;
    ind_1 = ind_2;
end

%% plot
figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
plot(time_vec*1e9,norm_voltage_vec,'Color',bRGY(8,:),'LineStyle','-','LineWidth',2)
hold on
plot(time_vec*1e9,norm_current_vec,'Color',bRGY(3,:),'LineStyle','-','LineWidth',2)
xlabel('Time [ns]','FontSize',fontSize,'FontName','Times')
ylabel('Signal [Arb.]','FontSize',fontSize,'FontName','Times')
set(gca,'FontSize',fontSize,'FontName',fontName)
% ylim([0 11])
% xlim([0 160])
legend('SQUID Voltage','Applied Current')
grid on

figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
plot(avg_time_vec*1e9,avg_current_vec*1e6,'Color',bRGY(8,:),'LineStyle','-','LineWidth',2)
xlabel('Time [ns]','FontSize',fontSize,'FontName','Times')
ylabel('Applied Current [\mu A]','FontSize',fontSize,'FontName','Times')
set(gca,'FontSize',fontSize,'FontName',fontName)
grid on

figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
plot(avg_time_vec*1e9,avg_voltage_vec*1e6,'Color',bRGY(8,:),'LineStyle','-','LineWidth',2)
xlabel('Time [ns]','FontSize',fontSize,'FontName','Times')
ylabel('SQUID voltage [\mu V]','FontSize',fontSize,'FontName','Times')
set(gca,'FontSize',fontSize,'FontName',fontName)
grid on

%final plots
figureCaptions = {sprintf('JJ I_c = %g uA',I_c*1e6),...
                  sprintf('JJ beta_c = %g',beta_c),...
                  sprintf('bias to SQUID I_b = %g uA',I_b*1e6),...
                  sprintf('I_b/2I_c = %g',I_b/(2*I_c)),...
                  sprintf('L_1 = %g pH',L1*1e12),...
                  sprintf('L_2 = %g pH',L2*1e12),...
                  sprintf('k = %g',k),...
                  sprintf('M = 2k(L1 L2)^{1/2} = %g pH',mutual_inductance*1e12),...
                  sprintf('factor of 2 in M comes from two MIs for symmetry'),...
                  };

figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
plot(avg_current_vec*1e6,avg_voltage_vec*1e6,'Color',bRGY(3,:),'LineStyle','-','LineWidth',2)
xlabel('Applied Current [\mu A]','FontSize',fontSize,'FontName','Times')
ylabel('SQUID Voltage [\mu V]','FontSize',fontSize,'FontName','Times')
set(gca,'FontSize',fontSize,'FontName',fontName)
grid on
k1 = gtext(figureCaptions(1:length(figureCaptions)));
set(k1,'FontSize',fontSize_legend,'FontName','Times')

figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
plot(mutual_inductance*avg_current_vec/p.Phi0,avg_voltage_vec*1e6,'Color',bRGY(3,:),'LineStyle','-','LineWidth',2)
xlabel('Applied Flux [\Phi/\Phi_0]','FontSize',fontSize,'FontName','Times')
ylabel('SQUID Voltage [\mu V]','FontSize',fontSize,'FontName','Times')
set(gca,'FontSize',fontSize,'FontName',fontName)
grid on
k1 = gtext(figureCaptions(1:length(figureCaptions)));
set(k1,'FontSize',fontSize_legend,'FontName','Times')
