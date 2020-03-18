%takemargins
%simple program to take margins. I will ramp the yoko current and get the
%traces. 
%in this case i also take an iv curve since then i can map the cleanliness
%of the margins to the flat area in the iv curve. 
%% quick IV taker, independent
filename='Margins_200MHz_BP_NONWA.mat'
path = 'C:\Users\mac3\Documents\DATA\IMS\2018_04_18_cooldown\chip35\Device_4\';
I_values = [-1.0e-3:100e-6:0.0e-3,3.5e-3:100e-6:7.5e-3];
%I_values_rand=I_values(randperm(length(I_values)));
%I_values=[-1.5e-3:100e-6:-1.0e-3,I_values_rand];
%I_values = [6.5e-3:50e-6:7.25e-3];
amplitude_8195=1000e-3;
%ramp downc
%I_values=[I_values(1:end-1),fliplr(I_values)];
V_values = [];

%% settings
VDwellTime = 1;
Irange = 0.010;
Vrange = 0.05;
compliance_voltage = 10;
nplc = 10;


%%
%Initialize Current source
fprintf(yokoobj,':SOUR:FUNC CURR')
fprintf(yokoobj,[':SOUR:RANG ',num2str(Irange)])  % sets current range
fprintf(yokoobj,[':SOUR:PROT:VOLT ',num2str(compliance_voltage)])  % sets complaince voltage
fprintf(yokoobj,'OUTP ON')
%reset current to 0:
fprintf(yokoobj, [':SOUR:LEV:FIX ',num2str(0)])
% turn off signal
set8195_runenable(obj8195,false)
pause(1)
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
set8195_runenable(obj8195,true)
pause(1)

%%
[outData] = SaveScopeData( scopeobj, 1,false);

margindata=zeros(length(I_values),length(outData));
t=outData(1,:);

for n = 1:length(I_values)
    
 %   sprintf('Starting sweep up')
    fprintf(yokoobj, [':SOUR:LEV:FIX ',num2str(I_values(n))])
    fprintf(nanovoltobj,'initiate')
    fprintf(nanovoltobj,'*TRG')
    if n==1
        pause(0.2);
    end
    
    %fetch voltage from multimeter
    V_values(n) = str2num(query(nanovoltobj,'fetch?'));
    %fetch the scope
    [outData] = SaveScopeData( scopeobj, 1,true);
    margindata(n,:)=outData(2,:);
    title(I_values(n));
end

fprintf(yokoobj, [':SOUR:LEV:FIX ',num2str(0e-3)])

%%
phi=2.067833831e-15;
%vexp=phi*1500*16e9;
vexp=  0.008287241007531;
figure(1);
hold all
%plot(V_values -mean(V_values(1:2)),I_values,vexp*ones(size(I_values)),I_values)
plot(V_values ,I_values,'o',vexp*ones(size(I_values)),I_values)
% 
cd(path)
save(filename)
