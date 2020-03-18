%create all the gpib and ethernet objects I need
instrreset
%% load yoko

yoko_GS200 = 1;
yokoobj = instrfind('Type', 'gpib', 'BoardIndex', 0, 'PrimaryAddress', yoko_GS200, 'Tag', '');
if isempty(yokoobj)
    yokoobj = gpib('NI', 0, yoko_GS200);
else
    fclose(yokoobj);
    yokoobj = yokoobj(1);
end
fopen(yokoobj);

%% load nanovoltmeter

address_V = 24;
nanovoltobj = instrfind('Type', 'gpib', 'BoardIndex', 0, 'PrimaryAddress', address_V, 'Tag', '');
if isempty(nanovoltobj)
    nanovoltobj = gpib('NI', 0, address_V);
else
    fclose(nanovoltobj);
    nanovoltobj = nanovoltobj(1);
end
fopen(nanovoltobj);

%% keithley
% Find a GPIB object.
keithobj = instrfind('Type', 'gpib', 'BoardIndex', 0, 'PrimaryAddress', 23, 'Tag', '');

% Create the GPIB object if it does not exist
% otherwise use the object that was found.
if isempty(keithobj)
    keithobj = gpib('NI', 0, 23);
else
    fclose(keithobj);
    keithobj = keithobj(1);
end

% Connect to instrument object, obj1.
fopen(keithobj);

