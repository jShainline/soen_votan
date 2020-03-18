function setI_keithley(Ibias,keithobj)


fprintf(keithobj, sprintf('CURR %f ',Ibias));
% Disconnect from instrument object, obj1.


