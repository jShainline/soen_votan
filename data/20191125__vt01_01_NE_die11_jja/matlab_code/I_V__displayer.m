%% quick IV taker, independent
I_val =3.6e-3;
%ramp downc
%I_values=[I_values(1:end-1),fliplr(I_values)];
V_values = [];

%% settings
VDwellTime = 0;
Irange = 0.015;
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
fprintf(yokoobj, [':SOUR:LEV:FIX ',num2str(I_val)])
for k = 1:1000
    
    %   sprintf('Starting sweep up')
    
    fprintf(nanovoltobj,'initiate')
    fprintf(nanovoltobj,'*TRG')
    V_values(k) = str2num(query(nanovoltobj,'fetch?'));
    
    
    %%
    phi=2.067833831e-15;
    %vexp=phi*1500*16e9;
    vexp= 0.008572621392121;
    figure(2);
    plot(V_values)
    title(V_values(end))
    pause(0.5)
end