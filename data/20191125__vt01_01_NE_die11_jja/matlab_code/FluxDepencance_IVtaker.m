%Flux depencance iv curves
phi=2.067833831e-15;
path = 'E:\DATA\Challenger_2\11_21_2019\';
filename='IV_vs_flux_moat_2.mat';
cd(path)
I={};
V={};
flux=linspace(28e-3,32e-3,41);
%flux=[linspace(-7e-3,0e-3,21),linspace(21e-3,25e-3,21)];
setIrange_keithley(50e-3,keithobj)
keithoutput('ON',keithobj)
setI_keithley(0,keithobj)

[Io,Vo]=IVfunction(yokoobj,nanovoltobj);
for k=1:length(flux)
    if k==1
        RampI_keithley(0,flux(k),keithobj)
    else
        setI_keithley(flux(k),keithobj)
    end
    pause(5)
    tic
    [I_values,V_values]=IVfunction(yokoobj,nanovoltobj);
    t=toc
    fprintf('Time remaining %f  seconds\n',(t*(length(flux)+1-k)))
    
    plot(Vo-mean(Vo(1:10)),Io)
    I{k}=I_values;
    V{k}=V_values;
    
    save(filename,'I','V','Io','Vo','flux')
    
end

RampI_keithley(flux(k),0,keithobj)
keithoutput('OFF',keithobj)