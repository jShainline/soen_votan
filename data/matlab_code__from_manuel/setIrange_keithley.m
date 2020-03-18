function setIrange_keithley(Irange,keithobj)

if abs(Irange)>100e-3
    disp('range too large, set to default')
    Irange=100e-3;
end
fprintf(keithobj, sprintf('CURR:RANG %f ',Irange));
% Disconnect from instrument object, obj1.


