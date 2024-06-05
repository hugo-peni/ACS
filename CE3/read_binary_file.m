function [y, r, u, d, anyRef, sample] = read_binary_file(path)
%READ_BINARY_FILE Reads binary (.bin) file generated during the experiment
% 
%   READ_BINARY_FILE(path) reads the file from path and returns the experiment
%       data.
% 
%   ----------------------------------------------------------------------------
%   Output:
%       y       % Measurements
%       r       % Reference signals
%       u       % Actuator commands
%       d       % Disturbance
%       anyRef  % Reference applied?    (true whenever reference is applied)
%       sample  % Sample number         (Always increases by 1)
% 
%   ----------------------------------------------------------------------------
%   Copyright 2022 Vaibhav Gupta, DDMAC, EPFL (MIT License)
% 

    arguments
        path
    end

    fileID = fopen(path, 'r');
    
    % Read data as 'double' in little-endian ordering
    data = fread(fileID, 'double', 'l');

    % Reshape data into (:, n) array
    n = 8;  % Number of columns in data
    data = reshape(data', [n, numel(data)/n])';

    % Extract
    y = data(:, 1:2);
    r = data(:, 3:4);
    u = data(:, 5);
    d = data(:, 6);
    anyRef = logical(data(:,7));
    sample = data(:,8);

    fclose(fileID);

end