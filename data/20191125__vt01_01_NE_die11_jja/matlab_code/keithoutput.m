function keithoutput(stout,keithobj)

if strcmp(stout,'ON')
    fprintf(keithobj, sprintf('OUTP %s',stout));
elseif strcmp(stout,'OFF')
    fprintf(keithobj, sprintf('OUTP %s',stout));
else
    error('wrong output state')
end



