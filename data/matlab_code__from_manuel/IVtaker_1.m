%% quick IV taker, independent
I_values = 0:0.2e-6:100e-6;
%ramp downc
I_values=[I_values(1:end-1),fliplr(I_values)];
V_values = [];

%% settings
VDwellTime = 0;
Irange = max(I_values);
Vrange = 300e-3;%50*Irange;
compliance_voltage = 10;
nplc = 1;
path = 'E:\DATA\Jeff_chips\jja\2019_11_27\';

%%
%Initialize Current source
fprintf(yokoobj,':SOUR:FUNC CURR')
fprintf(yokoobj,[':SOUR:RANG ',num2str(Irange)])  % sets current range
fprintf(yokoobj,[':SOUR:PROT:VOLT ',num2str(compliance_voltage)])  % sets complaince voltage
fprintf(yokoobj,'OUTP ON')
%reset current to 0:
fprintf(yokoobj, [':SOUR:LEV:FIX ',num2str(0)])
% turn off signal
%set8195_runenable(obj8195,false)
%pause(1)
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
%turn on signal
%set8195_runenable(obj8195,true)
%pause(1)

%%

for n = 1:length(I_values)
    
 %   sprintf('Starting sweep up')
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
close all
phi=2.067833831e-15;
vexp=phi*1500*16e9;

% uiopen('E:\DATA\cooldown_march_2019\wafer180327B_Chip24\2019-05-17_ajs\IV_4_7K_NoPower.fig',1)
% uiopen('E:\DATA\cooldown_march_2019\wafer180327B_Chip24\2019-05-17_ajs\IV_3_1K_NoPower.fig',1)
hold on
%offset=mean(V_values(and(I_values>-2.7e-3,I_values<-1.0e-3) ) );
offset=mean(V_values(1:5));
%plot(V_values-0*mean(V_values(1:2))+7.5e-6,I_values,-vexp*ones(size(I_values)),I_values)
%plot(V_values-offset,I_values,vexp*ones(size(I_values)),I_values)
%plot(V_values-offset,I_values,vexp*ones(size(I_values)),I_values)

R=diff(V_values)./diff(I_values);  
%plot(V_values-I_values*mean(R((I_values<5e-3))),I_values,'m')
plot((V_values-offset)*1e3,I_values*1e3,'.-b')
xlabel('mV')
ylabel('mA')
%axis([0,150e-3,0,Inf]) 
cd(path)
% %%
 
 
 R1=diff(V_values)./diff(I_values);  
%% 
 figure(2)
 
 hold on
 plot(I_values(1:end-1),R1,'b')
 
