function setpower(pow)


% Find a GPIB object.
obj1 = instrfind('Type', 'gpib', 'BoardIndex', 0, 'PrimaryAddress', 19, 'Tag', '');

% Create the GPIB object if it does not exist
% otherwise use the object that was found.
if isempty(obj1)
    obj1 = gpib('ni', 0, 19);
else
    fclose(obj1);
    obj1 = obj1(1);
end

% Connect to instrument object, obj1.
fopen(obj1);

%% Disconnect and Clean Up
if pow> 25 %dbm
    error('power too large')
end
fprintf(obj1, sprintf(':POW %f dBm',pow));
% Disconnect from instrument object, obj1.

 fclose(obj1);


