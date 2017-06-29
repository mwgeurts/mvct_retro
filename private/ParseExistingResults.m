function results = ParseExistingResults(filename)
% ParseExistingResults attempts to load an existing set of results into a
% return variable. If not found, it will create an empty CSV file to begin
% writing results to.
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

% Open file handle to current results .csv set
fid = fopen(filename, 'r');

% If a valid file handle was returned
if fid > 0
    
    % Log loading of existing results
    Event('Found results file');
    
    % Scan results .csv file for the following format of columns (see
    % documentation above for the results file format)
    results = textscan(fid, '%s %s %s %s %s %s %s %s %s %s %s %s %s %s', ...
        'Delimiter', {','}, 'commentStyle', '#');
    
    % Close the file handle
    fclose(fid);
    
    % Log completion
    Event(sprintf('%i results loaded from %s', size(results{1}, 1) - 1, ...
        filename));

% Otherwise, create new results file, saving column headers
else
    
    % Log generation of new file
    Event(['Generating new results file ', filename]);
    
    % Open write file handle to current results set
    fid = fopen(filename, 'w');
    
    % Print version information
    fprintf(fid, '# TomoTherapy MVCT Retrospective Analysis Tool\n');
    fprintf(fid, '# Author: Mark Geurts <mark.w.geurts@gmail.com>\n');
    fprintf(fid, ['# See ScanArchives.m and README.md for ', ...
        'more information on the format of this results file\n']);
    
    % Print column headers
    fprintf(fid, 'Archive,');
    fprintf(fid, 'SHA1,');
    fprintf(fid, 'MVCT_UID,');
    fprintf(fid, 'Plan_Name,');
    fprintf(fid, 'Scan_Length,');
    fprintf(fid, 'User_Pitch,');
    fprintf(fid, 'User_Yaw,');
    fprintf(fid, 'User_Roll,');
    fprintf(fid, 'User_X,');
    fprintf(fid, 'User_Y,');
    fprintf(fid, 'User_Z,');
    fprintf(fid, 'Version,');
    fprintf(fid, 'Reg_Method,');
    fprintf(fid, 'Reg_Pitch,');
    fprintf(fid, 'Reg_Yaw,');
    fprintf(fid, 'Reg_Roll,');
    fprintf(fid, 'Reg_X,');
    fprintf(fid, 'Reg_Y,');
    fprintf(fid, 'Reg_Z,');
    fprintf(fid, 'Similarity_Metric,');
    fprintf(fid, 'User_Similarity,');
    fprintf(fid, 'Reg_Similarity\n');

    % Close the file handle
    fclose(fid);
    
    % Return empty array
    results = [];
end

% Clear file hande
clear fid;