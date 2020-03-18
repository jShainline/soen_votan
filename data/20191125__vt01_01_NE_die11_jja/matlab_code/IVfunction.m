function [I_values,V_values]=IVfunction(yokoobj,nanovoltobj)
%close all
I_values = -2.5e-3:2e-6:2.5e-3;
%ramp downc
%I_values=[I_values(1:end-1),fliplr(I_values)];
V_values = [];

%% settings
nplc = 1;
VDwellTime = 0;
Irange = max(I_values);
%Vrange = 10*Irange;
Vrange = 0.001;
compliance_voltage = 10;
%path = 'C:\Users\mac3\Documents\DATA\IMS\2018_04_18_cooldown\Device_1\2018_04_20\';



%% Initialize Current source
fprintf(yokoobj,':SOUR:FUNC CURR')
fprintf(yokoobj,[':SOUR:RANG ',num2str(Irange)])  % sets current range
fprintf(yokoobj,[':SOUR:PROT:VOLT ',num2str(compliance_voltage)])  % sets complaince voltage
fprintf(yokoobj,'OUTP ON')
%reset current to 0:
fprintf(yokoobj, [':SOUR:LEV:FIX ',num2str(0)])

% Initialize Voltmeter
fprintf(nanovoltobj,'*RST');
fprintf(nanovoltobj,'status:preset');
fprintf(nanovoltobj,'*cls');
fprintf(nanovoltobj,['conf:volt:dc ',num2str(Vrange)]);
fprintf(nanovoltobj,'*sre 32') %not sure what this does either
fprintf(nanovoltobj,'input:filter off')
%The following line sets the rate at which the voltmeter averages the data over (ie sampling rate)
%0.02, 0.2, 1, 2, 10, 20, 100, 200 where 0.02 is the fastest sampling rate and lowest averaging
fprintf(nanovoltobj,['SENS:VOLT:DC:NPLC ',num2str(nplc)]) %number of power line cycles to average over
fprintf(nanovoltobj,'trigger:source bus')
fprintf(nanovoltobj,'trigger:delay 0')
fprintf(nanovoltobj,'initiate')
fprintf(nanovoltobj,'*TRG')
offset = str2num(query(nanovoltobj,'fetch?'));

%%
sprintf('Starting sweep up')
for n = 1:length(I_values)
    
    
    fprintf(yokoobj, [':SOUR:LEV:FIX ',num2str(I_values(n))])
    fprintf(nanovoltobj,'initiate')
    fprintf(nanovoltobj,'*TRG')
    if n==1
        pause(0.2);
    end
    V_values(n) = str2num(query(nanovoltobj,'fetch?'));
end


fprintf(yokoobj, [':SOUR:LEV:FIX ',num2str(0e-3)])

%%
%figure(1);
close all
%uiopen('C:\Users\mac3\Documents\DATA\IMS\2018_04_18_cooldown\chip35\Device_3\2018_05_16\IV_4_5k_withOptiLab_on.fig',1)
hold all
plot(V_values-mean(V_values(1:10)),I_values)
%  cd(path)
%  R=diff(V_values)./diff(I_values);  figure(2)
% plot(I_values(1:end-1),R)

