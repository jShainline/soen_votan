%create all the gpib and ethernet objects I need
instrreset
%% load yoko

yoko_GS200 = 1;
yokoobj = instrfind('Type', 'gpib', 'BoardIndex', 1, 'PrimaryAddress', yoko_GS200, 'Tag', '');
if isempty(yokoobj)
    yokoobj = gpib('NI', 1, yoko_GS200);
else
    fclose(yokoobj);
    yokoobj = yokoobj(1);
end
fopen(yokoobj);

%% load yoko flux

yoko_GS200FLUX = 5;
yokofluxobj = instrfind('Type', 'gpib', 'BoardIndex', 1, 'PrimaryAddress', yoko_GS200FLUX, 'Tag', '');
if isempty(yokofluxobj)
    yokofluxobj = gpib('NI', 1, yoko_GS200FLUX);
else
    fclose(yokofluxobj);
    yokofluxobj = yokofluxobj(1);
end
fopen(yokofluxobj);


%% load nanovoltmeter

address_V = 24;
nanovoltobj = instrfind('Type', 'gpib', 'BoardIndex', 1, 'PrimaryAddress', address_V, 'Tag', '');
if isempty(nanovoltobj)
    nanovoltobj = gpib('NI', 1, address_V);
else
    fclose(nanovoltobj);
    nanovoltobj = nanovoltobj(1);
end
fopen(nanovoltobj);
