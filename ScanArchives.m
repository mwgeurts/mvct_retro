function ScanArchives()
% ScanArchives scans for MVCT scans from patients archives in a  directory 
% (given by the confgi.txt variable ARCHIVE_PATH) and stores information
% about each MVCT to the CSV file specified by RESULTS_CSV. If a 
% configuration option is set for REGISTRATION_METHOD then for each MVCT 
% this function will re-register the images and store the new registration 
% adjustments as well.
%
% If an entry already exists for the patient archive (determined by SHA1
% signature), registration method, similarity metric, and MVCT (determined 
% by UID), the workflow will be skipped. In this manner, ScanArchives can 
% be run multiple times to analyze a large directory of archives.
%
% The RESULTS_CSV file contains the following columns. The registration and
% similarity columns will only be filled out if configured in config.txt.
%
%   {1}:  Full path to patient archive _patient.xml.  However, if 
%         the config option ANON_RESULTS is set to 1, will be empty.
%   {2}:  SHA1 signature of _patient.xml file
%   {3}:  MVCT Timestamp (MATLAB datenum)
%   {4}:  MVCT UID
%   {5}:  Plan Name
%   {6}:  Scan Length (cm)
%   {7}:  Time from scan to treatment (minutes)
%   {8}:  Number of MVCT scans performed on same day
%   {9}:  User Registered Pitch (degrees)
%   {10}:  User Registered Yaw (degrees)
%   {11}:  User Registered Roll (degrees)
%   {12}: User Registered X Translation (cm)
%   {13}: User Registered Y Translation (cm)
%   {14}: User Registered Z Translation (cm)
%   {15}: Tool Version
%   {16}: Registration Method, if provided
%   {17}: Re-Registered Pitch (degrees)
%   {18}: Re-Registered Yaw (degrees)
%   {19}: Re-Registered Roll (degrees)
%   {20}: Re-Registered X Translation (cm)
%   {21}: Re-Registered Y Translation (cm)
%   {22}: Re-Registered Z Translation (cm)
%   {23}: Similarity Metric, if provided
%   {24}: User Registration Similarity
%   {25}: Re-Registration Similarity
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

%% Set runtime variables
% Turn off MATLAB warnings
warning('off','all');

% Set version handle
version = '1.0.3';

% Determine path of current application
[path, ~, ~] = fileparts(mfilename('fullpath'));

% Set current directory to location of this application
cd(path);

% Clear temporary variable
clear path;

%% Initialize Log
% Set version information.  See LoadVersionInfo for more details.
versionInfo = LoadVersionInfo;

% Store program and MATLAB/etc version information as a string cell array
string = {'TomoTherapy MVCT Retrospective Analysis Tool'
    sprintf('Version: %s (%s)', version, versionInfo{6});
    sprintf('Author: Mark Geurts <mark.w.geurts@gmail.com>');
    sprintf('MATLAB Version: %s', versionInfo{2});
    sprintf('MATLAB License Number: %s', versionInfo{3});
    sprintf('Operating System: %s', versionInfo{1});
    sprintf('CUDA: %s', versionInfo{4});
    sprintf('Java Version: %s', versionInfo{5})
};

% Add dashed line separators      
separator = repmat('-', 1,  size(char(string), 2));
string = sprintf('%s\n', separator, string{:}, separator);

% Log information
Event(string, 'INIT');

% Clear temporary variables
clear string separator versionInfo;

% Add archive extraction submodule to search path
addpath('./tomo_extract');

% Check if MATLAB can find FindMVCTScans
if exist('FindMVCTScans', 'file') ~= 2
    
    % If not, throw an error
    Event(['The tomo_extract submodule does not exist in the ', ...
        'search path. Use git clone --recursive or git submodule init ', ...
        'followed by git submodule update to fetch all submodules'], ...
        'ERROR');
end

% Add fusion_tools submodule to search path
addpath('./fusion_tools');

% Check if MATLAB can find RigidRegister
if exist('RigidRegister', 'file') ~= 2
    
    % If not, throw an error
    Event(['The fusion_tools submodule does not exist in the ', ...
        'search path. Use git clone --recursive or git submodule init ', ...
        'followed by git submodule update to fetch all submodules'], ...
        'ERROR');
end

% Load configuration settings
config = ParseConfigOptions('config.txt');

% Load Results .csv
results = ParseExistingResults(config.RESULTS_CSV, config);

%% Start scanning for archives
% Note beginning execution
Event(['ScanArchives beginning search of ', config.ARCHIVE_PATH, ...
    ' for patient archives']);

% Retrieve folder contents of input directory
folderList = dir(config.ARCHIVE_PATH);

% Shuffle random number generator seed
rng shuffle;

% Randomize order of folder list
folderList = folderList(randperm(size(folderList, 1)), :);

% Initialize folder counter
i = 0;

% Initialize plan counter
count = 0;

% Start AutoSystematicError timer
t = tic;

% Start recursive loop through each folder, subfolder
while i < size(folderList, 1)
    
    % Increment current folder being analyzed
    i = i + 1;
    
    % If the folder content is . or .., skip to next folder in list
    if strcmp(folderList(i).name, '.') || strcmp(folderList(i).name, '..')
        continue
        
    % Otherwise, if the folder content is a subfolder    
    elseif folderList(i).isdir == 1
        
        % Retrieve the subfolder contents
        subFolderList = dir(fullfile(config.ARCHIVE_PATH, ...
            folderList(i).name));
        
        % Randomize order of subfolder list
        subFolderList = subFolderList(randperm(size(subFolderList, 1)), :);
        
        % Look through the subfolder contents
        for j = 1:size(subFolderList, 1)
            
            % If the subfolder content is . or .., skip to next subfolder 
            if strcmp(subFolderList(j).name, '.') || ...
                    strcmp(subFolderList(j).name, '..')
                continue
            else
                
                % Otherwise, replace the subfolder name with its full
                % reference
                subFolderList(j).name = fullfile(folderList(i).name, ...
                    subFolderList(j).name);
            end
        end
        
        % Append the subfolder contents to the main folder list
        folderList = vertcat(folderList, subFolderList); %#ok<AGROW>
        
        % Clear temporary variable
        clear subFolderList;
        
    % Otherwise, if the folder content is a patient archive
    elseif size(strfind(folderList(i).name, '_patient.xml'), 1) > 0
        
        % Generate a SHA1 signature for the archive patient XML file using
        % the shasum system command on Unix/Mac, or sha1sum on Windows
        % (provided as part of this repository)
        if ispc
            [~, cmdout] = system(['sha1sum "', ...
                fullfile(config.ARCHIVE_PATH, folderList(i).name), '"']);
        else
            [~, cmdout] = system(['shasum "', ...
                fullfile(config.ARCHIVE_PATH, folderList(i).name), '"']);
        end
        
        % Save just the 40-character signature
        sha = cmdout(1:40);
        
        % Log patient XML and SHA1 signature
        Event(['Found patient archive ', folderList(i).name, ...
            ' with SHA1 signature ', sha]);
        
        % Clear temporary variable
        clear cmdout;

        % Generate separate path and name variables for XML
        [path, name, ext] = ...
            fileparts(fullfile(config.ARCHIVE_PATH, folderList(i).name));
        name = strcat(name, ext);
        
        % Clear temporary variable
        clear ext;
        
        % Attempt to find MVCT scans
        try 
        
            % Search for all MVCT scans in the archive
            scans = FindMVCTScans(path, name);
            
            % Search for all treatments in the archive
            txs = FindTreatments(path, name);
        
        % If an error is thrown, catch
        catch exception

            % Report exception to error log
            Event(getReport(exception, 'extended', 'hyperlinks', ...
                'off'), 'CATCH');

            % Continue to next file
            continue;
        end
                
        % Loop through each plan
        Event('Looping through each MVCT scan');
        for j = 1:length(scans)
            
            % If a registration method is set
            if isfield(config, 'REGISTRATION_METHOD') || ...
                    isfield(config, 'SIMILARITY_METRIC')
                
                % Attempt to load plan
                try
                    
                    % Load the planning CT
                    reference = LoadImage(path, name, scans{j}.planUID);
                    
                % If an error is thrown, catch
                catch exception
                    
                    % Report exception to error log
                    Event(getReport(exception, 'extended', 'hyperlinks', ...
                        'off'), 'CATCH');

                    % Continue to next image set
                    continue;
                end
            end
            
            % Loop through each scan
            for k = 1:length(scans{j}.scanUIDs)
            
                % Initialize flag to indicate whether the current plan
                % already contains contents in RESULTS_CSV
                found = false;

                % If the results .csv exists and was loaded above
                if exist('results', 'var') && ~isempty(results)

                    % Loop through each result
                    for l = 2:size(results{1}, 1)

                        % If the XML MVCT UID & versions match
                        if strcmp(results{3}{l}, scans{j}.scanUIDs{k}) && ...
                                strcmp(results{12}{l}, version)

                            % Set the flag to true, since a match was found
                            found = true;

                            % Break the loop to stop searching
                            break;
                        end
                    end

                    % Clear temporary variable
                    clear l;
                end

                % If results exist for this daily image, continue
                if found
                    
                    % Log matching UID
                    Event(['UID ', scans{j}.scanUIDs{k}, ...
                        ' skipped as results were found in ', ...
                        config.RESULTS_CSV]);
                    
                    continue;
                end

                % Attempt to run registration
                try 

                    % If a registration method is set
                    if isfield(config, 'REGISTRATION_METHOD') || ...
                            isfield(config, 'SIMILARITY_METRIC')
                        
                        % Load the daily image
                        daily = LoadDailyImage(path, 'ARCHIVE', name, ...
                            scans{j}.scanUIDs{k});
                    end

                    % If a registration method is set
                    if isfield(config, 'REGISTRATION_METHOD')

                        % Register the daily image
                        daily.rigid = RigidRegister(reference, daily, ...
                            'method', config.REGISTRATION_METHOD, ...
                            'levels', config.REGISTRATION_LEVELS, ...
                            'iterations', config.REGISTRATION_ITER, ...
                            'metric', config.REGISTRATION_METRIC, ...
                            'bone', config.REGISTER_BONE);
                    end
                    
                    % If a similarity metric is set
                    if isfield(config, 'SIMILARITY_METRIC')
                        
                        % Convert reference image to equivalent daily-IVDT 
                        % image
                        r = interp1(daily.ivdt(:,2), daily.ivdt(:,1), ...
                            interp1(reference.ivdt(:,1), ...
                            reference.ivdt(:,2), reference.data, ...
                            'linear', 'extrap'), 'linear', 'extrap');

                        % Merge the daily image to the reference image
                        % coordinates using the user
                        m = MergeImages(reference, daily, ...
                            daily.registration);

                        % Find only the slices that correspond to the
                        % MVCT
                        s = find(sum(sum(m.mask, 1), 2));
                        
                        % Log action
                        Event('Computing similarity metrics');
                        
                        % Calculate image similarity metric
                        switch config.SIMILARITY_METRIC
                            
                            % Structural Similarity Index
                            case 'SSI'

                                % Compute the SSI on only the masked slices
                                daily.user_similarity = ...
                                    ssim(r(:,:,min(s):max(s)), ...
                                    m.data(:,:,min(s):max(s)));
                                
                                % Log result
                                Event(sprintf(['User registration structural', ...
                                    ' simiarlity index = %f'], ...
                                    daily.user_similarity));
                            
                            % Mean square error
                            case 'MSE'

                                % Compute the MSE on only the masked slices
                                daily.user_similarity = ...
                                    immse(r(:,:,min(s):max(s)), ...
                                    m.data(:,:,min(s):max(s)));
                                
                                % Log result
                                Event(sprintf(['User registration mean', ...
                                    ' square error = %f'], ...
                                    daily.user_similarity));
                        end
                        
                        % If re-registration data does exists
                        if isfield(daily, 'rigid')
                        
                            % Merge the daily image to the reference image
                            % coordinates using the user
                            m = MergeImages(reference, daily, ...
                                daily.rigid);

                            % Find only the slices that correspond to the
                            % MVCT
                            s = find(sum(sum(m.mask, 1), 2));

                            % Log action
                            Event('Computing similarity metrics');
                            
                            % Calculate image similarity metric
                            switch config.SIMILARITY_METRIC

                                % Structural Similarity Index
                                case 'SSI'

                                    % Compute the SSI on only the masked slices
                                    daily.re_similarity = ...
                                        ssim(r(:,:,min(s):max(s)), ...
                                        m.data(:,:,min(s):max(s)));
                                    
                                    % Log result
                                    Event(sprintf(['Re-registration structural', ...
                                        ' simiarlity index = %f'], ...
                                        daily.re_similarity));

                                % Mean square error
                                case 'MSE'

                                    % Compute the SSI on only the masked slices
                                    daily.re_similarity = ...
                                        immse(r(:,:,min(s):max(s)), ...
                                        m.data(:,:,min(s):max(s)));
                                    
                                    % Log result
                                    Event(sprintf(['Re-registration mean', ...
                                        ' square error = %f'], ...
                                        daily.re_similarity));
                            end
                        end
                        
                        % Clear temporary variables
                        clear r m s;
                    end

                    %% Append Results
                    % Log start
                    Event(['Writing results to ', config.RESULTS_CSV]);

                    % Open append file handle to results .csv
                    fid = fopen(config.RESULTS_CSV, 'a');

                    % If anon is TRUE, do not store the XML name and 
                    % location in column 1
                    if config.ANON_RESULTS

                        % Instead, replace with 'ANON'
                        fprintf(fid,'ANON,');
                    else

                        % Otherwise, write relative path location 
                        fprintf(fid, '%s,', ...
                            strrep(folderList(i).name, ',', ''));
                    end

                    % Write XML SHA1 signature in column 2
                    fprintf(fid, '%s,', sha);
                    
                    % Write MVCT imtestamp in column 3
                    fprintf(fid, '%f,', datenum([scans{j}.date{k}, ...
                        'T', scans{j}.time{k}], 'yyyymmddTHHMMSS'));

                    % Write MVCT UID in column 4
                    fprintf(fid, '%s,', scans{j}.scanUIDs{k});

                    % Write plan name in column 5
                    fprintf(fid, '%s,', ...
                        strrep(scans{j}.planName, ',', ' '));

                    % Write scan length in column 6
                    fprintf(fid, '%f,', ...
                        abs(diff(scans{j}.scanLengths(k,:))));

                    % Write time to treatment in column 7, in minutes
                    t = 0;
                    for l = 1:length(txs{j}.date)
                        d = datenum([txs{j}.date{l}, 'T', ...
                            txs{j}.time{l}], 'yyyymmddTHHMMSS') - ...
                            datenum([scans{j}.date{k}, 'T', ...
                            scans{j}.time{k}], 'yyyymmddTHHMMSS');
                        if (d > 0 && d < 1) && (t == 0 || d < t)
                            t = d;
                        end
                    end
                    fprintf(fid, '%f,', t*24*60);
                    clear t l d;
                    
                    % Write number of MVCT scans in column 8
                    c = 1;
                    for l = 1:length(scans{j}.date)
                        if strcmp(scans{j}.date{k}, scans{j}.date{l})
                            c = c + 1;
                        end
                    end
                    for l = 1:length(txs{j}.date)
                        if strcmp(scans{j}.date{k}, txs{j}.date{l})
                            c = c - 1;
                        end
                    end
                    fprintf(fid, '%f,', max(c, 1));
                    clear c l;
                
                    % Write user registration in columns 9-14
                    fprintf(fid, '%f,%f,%f,%f,%f,%f,', ...
                        scans{j}.registration(k,:));
                    
                    % Write version in column 15
                    fprintf(fid, '%s,', version);
                    
                    % If a new registration was performed
                    if isfield(daily, 'rigid')

                        % Write registration method in column 16
                        fprintf(fid, '%s,', [config.REGISTRATION_METHOD, ...
                            '_', config.REGISTRATION_METRIC]);

                        % Write new registration values in columns 17-22
                        fprintf(fid, '%f,%f,%f,%f,%f,%f,', ...
                            daily.rigid);
                    end
                    
                    % If a similarity metric was calculated
                    if isfield(daily, 'user_similarity')
                        
                        % Write similarity metric in column 23
                        fprintf(fid, '%s,', config.SIMILARITY_METRIC);
                        
                        % Write user metric in column 24
                        fprintf(fid, '%f,', daily.user_similarity);
                    end

                    % If a new similarity metric was calculated
                    if isfield(daily, 're_similarity')
                        
                        % Write user metric in column 25
                        fprintf(fid, '%f,', daily.re_similarity);
                    end

                    % Write new line
                    fprintf(fid, '\n');
                    
                    % Close file handle
                    fclose(fid);

                    % Clear temporary variables
                    clear fid daily;

                    % Log start
                    Event(sprintf(['Completed workflow', ...
                        ' on MVCT UID %s'], scans{j}.scanUIDs{k}));

                    % Increment the count of processed images
                    count = count + 1;

                % If an error is thrown, catch
                catch exception

                    % Report exception to error log
                    Event(getReport(exception, 'extended', 'hyperlinks', ...
                        'off'), 'CATCH');

                    % Continue to next image set
                    continue;
                end
            end
            
            % Clear temporary variables
            clear reference;
        end
        
        % Clear temporary variables
        clear path name scans sha txs;
    end 
end

% Log completion of script
Event(sprintf(['ScanArchives completed in %0.0f minutes, ', ...
    'processing %i MVCT scans'], toc(t)/60, count));

% Clear temporary variables
clear i j k t count;
