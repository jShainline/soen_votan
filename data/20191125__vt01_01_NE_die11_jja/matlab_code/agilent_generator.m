function agilent_generator()


% Find a GPIB object.
obj1 = instrfind('Type', 'gpib', 'BoardIndex', 7, 'PrimaryAddress', 1, 'Tag', '');

% Create the GPIB object if it does not exist
% otherwise use the object that was found.
if isempty(obj1)
    obj1 = gpib('KEYSIGHT', 7, 1);
else
    fclose(obj1);
    obj1 = obj1(1);
end

% Connect to instrument object, obj1.
fopen(obj1);

%% Disconnect and Clean Up





% Disconnect from instrument object, obj1.

 fclose(obj1);


