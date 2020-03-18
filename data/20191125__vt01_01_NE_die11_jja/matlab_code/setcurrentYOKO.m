function setcurrentYOKO(yokoobj,Ivalue)


if abs(Ivalue)>20e-3
    
    error('value too high')
end

fprintf(yokoobj, [':SOUR:LEV:FIX ',num2str(Ivalue)])
