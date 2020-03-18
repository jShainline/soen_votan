function RampI_keithley(Io,If,keithobj)
Ilist=linspace(Io,If,21);
for k=1:length(Ilist)
    
    setI_keithley(Ilist(k),keithobj)
    pause(0.25)
end


