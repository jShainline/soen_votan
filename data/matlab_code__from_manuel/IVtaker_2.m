%% quick IV taker, independent
% this assumes there is a small resitance due to measuring outside with no
% inductive taps
%
I_values = 0e-3:200e-6:20e-3;
%ramp downc
%I_values=[I_values(1:end-1),fliplr(I_values)];
V_values = [];

%% settings
VDwellTime = 0;
Irange = 0.015;
Vrange = 0.06;
compliance_voltage = 10;
nplc = 1;
path = 'C:\Users\mac3\Documents\DATA\IMS\dunktest\2018_06_20\';


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
phi=2.067833831e-15;
figure(1);
hold all
vexp=phi*16e9*480;
offsetind=and(I_values>-7.5e-3,I_values<7.5e-3) ;
plot(V_values-offset,I_values,'b',vexp*ones(size(I_values)),I_values)
% 
R=diff(V_values)./diff(I_values); 
figure(2)

Roffset=mean(R(offsetind));
plot(I_values(1:end-1),R-Roffset)

%
Vnew=V_values-I_values*Roffset;
offset=mean(Vnew(offsetind));
Vnew=Vnew-offset;
figure(1)
plot(Vnew,I_values,'k')

cd(path)
