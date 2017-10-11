function config = ParseConfigOptions(filename)
% ParseConfigOptions is executed by ScanArchives to open the config file
% and update the application settings.
% 
% Author: Mark Geurts, mark.w.geurts@gmail.com
% Copyright (C) 2017 University of Wisconsin Board of Regents
%
% This program is free software: you can redistribute it and/or modify it 
% under the terms of the GNU General Public License as published by the  
% Free Software Foundation, either version 3 of the License, or (at your 
% option) any later version.
%
% This program is distributed in the hope that it will be useful, but 
% WITHOUT ANY WARRANTY; without even the implied warranty of 
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General 
% Public License for more details.
% 
% You should have received a copy of the GNU General Public License along 
% with this program. If not, see http://www.gnu.org/licenses/.

% Log event and start timer
t = tic;
Event(['Opening file handle to config file ', filename]);

% Initialize return variable
config = struct;

% Open file handle to config.txt file
fid = fopen(filename, 'r');

% Verify that file handle is valid
if fid < 3
    
    % If not, throw an error
    Event(['The ', filename, ' file could not be opened. Verify that this ', ...
        'file exists in the working directory. See documentation for ', ...
        'more information.'], 'ERROR');
end

% Scan config file contents
c = textscan(fid, '%s', 'Delimiter', '=');

% Close file handle
fclose(fid);

% Loop through textscan array, separating key/value pairs into array
for i = 1:2:length(c{1})
    config.(strtrim(c{1}{i})) = strtrim(c{1}{i+1});
end

% Clear temporary variables
clear c i fid;

% Log event and completion
Event(sprintf('Configuration options loaded successfully in %0.3f seconds', ...
    toc(t)));

% Remove similarity or registration inputs if they are invalid
if isfield(config, 'REGISTRATION_METHOD') && ...
        ~strcmp(config.REGISTRATION_METHOD, 'MATLAB') && ...
        ~strcmp(config.REGISTRATION_METHOD, 'PLASTIMATCH')
    config = rmfield(config, 'REGISTRATION_METHOD');
end

if isfield(config, 'SIMILARITY_METRIC') && ...
        ~strcmp(config.SIMILARITY_METRIC, 'SSI') && ...
        ~strcmp(config.SIMILARITY_METRIC, 'MSE')
    config = rmfield(config, 'SIMILARITY_METRIC');
end

% Reformat numerical and logical inputs
if isfield(config, 'ANON_RESULTS')
    config.ANON_RESULTS = logical(str2double(config.ANON_RESULTS));
else
    config.ANON_RESULTS = false;
end

if isfield(config, 'REGISTER_BONE')
    config.REGISTER_BONE = logical(str2double(config.REGISTER_BONE));
else
    config.REGISTER_BONE = false;
end

if isfield(config, 'ALLOW_ROTATIONS')
    config.ALLOW_ROTATIONS = logical(str2double(config.ALLOW_ROTATIONS));
else
    config.ALLOW_ROTATIONS = true;
end

if isfield(config, 'REGISTRATION_LEVELS')
    config.REGISTRATION_LEVELS = str2double(config.REGISTRATION_LEVELS);
else
    config.REGISTRATION_LEVELS = 3;
end

if isfield(config, 'REGISTRATION_ITER')
    config.REGISTRATION_ITER = str2double(config.REGISTRATION_ITER);
else
    config.REGISTRATION_ITER = 30;
end
