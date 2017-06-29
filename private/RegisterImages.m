function rigid = RegisterImages(reference, daily, varargin)
% RegisterImages rigidly registers a daily image to a reference image. 
% The reference image is converted to daily-IVDT equivalent Hounsfield 
% Units by first interpolating to density using the reference IVDT and then 
% subsequently interpolating back to HU using the daily IVDT.  The final 
% merged image is therefore in daily-equivalent Hounsfield Units.
%
% The following variables are required for proper execution: 
%
%   reference:  structure containing the image data, dimensions, width,
%               start coordinates, structure set UID, couch checksum and 
%               IVDT.  See LoadImage.
%   daily:      structure containing the image data, dimensions, width,
%               start coordinates and IVDT.  See LoadDailyImage.
%
% In addition, the following options may be provided as name/value pairs:
%
%   method:     string contianing the algorithm to use for registration.  
%               Can be 'PLASTIMATCH' or 'MATLAB'
%   metric:     string containing metric. Can be 'MSE' 'GM' (plastimatch 
%               only) or 'MI'. If not provided, will default to 'MSE'
%   bone:       logical indicating whether to mask registration to only 
%               bony anatomy (true) or full image (false). If not provided,
%               will default to full image.
%   levels:     the number of levels to use during optimization.
%               Plastimatch is configured to support up to 3 levels.
%   iterations: number of iterations to run at each level.
%
% The following variables are returned upon succesful completion:
%
%   rigid:      six element vector of registration parameters [pitch yaw 
%               roll x y z], where angles are in radians and distances are
%               in cm.
%
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

% Initialize default options
method = 'MATLAB';
metric = 'MSE';
bone = false;
levels = 3;
iterations = 30;

% Loop through remaining input arguments
for i = 3:2:length(varargin)

    if strcmpi(varargin{i}, 'method')
        method = varargin{i+1};
    elseif strcmpi(varargin{i}, 'metric')
        metric = varargin{i+1};
    elseif strcmpi(varargin{i}, 'bone')
        bone = varargin{i+1};
    elseif strcmpi(varargin{i}, 'levels')
        levels = varargin{i+1};
    elseif strcmpi(varargin{i}, 'iterations')
        iterations = varargin{i+1};
    end
end

% If the bone flag is set
if bone
    
    % Note use of bony anatomy to event log
    Event('Merging reference and daily images using bony anatomy');
else
    
    % Note use of full image to event log
    Event('Merging reference and daily images using full image');
end

% Start timer
t = tic;
    
% Log which registration method was chosen
Event(['Method ', method, ' selected for registration']);

% Convert reference image to equivalent daily-IVDT image
reference.data = interp1(daily.ivdt(:,2), daily.ivdt(:,1), ...
    interp1(reference.ivdt(:,1), reference.ivdt(:,2), ...
    reference.data, 'linear', 'extrap'), 'linear', 'extrap');

% Note conversion in log
Event(['Reference image converted to daily-equivalent Hounsfield ', ...
    ' Units using IVDT']);

% Execute registration based on method variable
switch method

%% Use Plastimatch 6-DOF rigid registration
% This rigid registration technique requires plastimatch be installed
% on this workstation.  Data is passed to/from plastimatch using the ITK 
% .mha file format and text command file. See 
% http://iopscience.iop.org/0031-9155/55/21/001 for additional details
% on the algorithm.
case 'PLASTIMATCH'
    
    %% Build reference MHA file
    % Generate a temprary filename for the reference image
    refFilename = [tempname, '.mha'];
    
    % Open a write file handle to the temporary reference image
    fid = fopen(refFilename, 'w', 'l');
    
    % Start writing the ITK header
    fprintf(fid, 'ObjectType=Image\n');
    fprintf(fid, 'NDims=3\n');
    
    % Specify the dimensions of the reference image
    fprintf(fid, 'DimSize=%i %i %i\n', reference.dimensions);
    
    % Specify the data format (USHORT referring to unsigned 16-bit integer)
    fprintf(fid,'ElementType=MET_USHORT\n');
    
    % Specify the byte order as little
    fprintf(fid,'ElementByteOrderMSB=False\n');
    
    % Specify the reference voxel widths (in mm)
    fprintf(fid, 'ElementSize=%i %i %i\n', reference.width*10);
    
    % Specify the reference voxel spacing to equal the widths (in mm)
    fprintf(fid, 'ElementSpacing=%i %i %i\n', reference.width*10);
    
    % Specify the coordinate frame origin (in mm)
    fprintf(fid, 'Origin=%i %i %i\n', reference.start*10);
    
    % Complete the .mha file header
    fprintf(fid, 'ElementDataFile=LOCAL\n');
    
    % Write the reference image data to the temporary file as uint16
    fwrite(fid, reference.data, 'ushort', 0, 'l');
    
    % Close the file handle
    fclose(fid);
    
    % Clear the temporary variable
    clear fid;
    
    % Log where the reference file was saved
    Event(['Reference image written to ', refFilename]);
    
    %% Build mask for reference image (excluding outside MVCT FOV)
    % Initialize null array of the same size as the reference image
    refMask = zeros(reference.dimensions);
    
    % Create meshgrid the same size as one image
    [x,y] = meshgrid(reference.start(1):reference.width(1):...
        reference.start(1) + reference.width(1) * ...
        (reference.dimensions(1) - 1), reference.start(2):...
        reference.width(2):reference.start(2) + ...
        reference.width(2) * (reference.dimensions(2) - 1));
    
    % Loop through each reference image slice
    for i = 1:reference.dimensions(3)
        
        % If the reference slice IEC-Y coordinate value is within the daily
        % image slice range
        if reference.start(3)+(i*reference.width(3)) > ...
                daily.start(3) && reference.start(3) + (i * ...
                reference.width(3)) < daily.start(3) + ...
                daily.dimensions(3) * daily.width(3)
            
            % Set the mask to 1 within the daily image FOV
            refMask(:,:,i) = sqrt(x.^2+y.^2) < dailyFOV/2 - 0.1;
        end
    end
    
    % If the bone flag is enabled
    if bone
        
        % Update the mask to only include values above 176 HU
        refMask = refMask .* ...
            ceil((reference.data - 1200) / 65535);
    else
        % Otherwise, set mask to exclude image noise (values below -824 HU)
        refMask = refMask .* ...
            ceil((reference.data - 200) / 65535);
    end
    
    % Generate a temporary file name for the reference image mask
    refMaskFilename = [tempname, '.mha'];
    
    % Open a write file handle to the temporary reference image mask
    fid = fopen(refMaskFilename, 'w', 'l');
    
    % Start writing the ITK header
    fprintf(fid,'ObjectType=Image\n');
    fprintf(fid,'NDims=3\n');
    
    % Specify the dimensions of the reference image
    fprintf(fid, 'DimSize=%i %i %i\n', reference.dimensions);
    
    % Specify the data format (USHORT referring to unsigned 16-bit integer)
    fprintf(fid,'ElementType=MET_USHORT\n');
    
    % Specify the byte order as little
    fprintf(fid,'ElementByteOrderMSB=False\n');
    
    % Specify the reference voxel widths (in mm)
    fprintf(fid, 'ElementSize=%i %i %i\n', reference.width*10);
    
    % Specify the merged voxel spacing to equal the widths (in mm)
    fprintf(fid, 'ElementSpacing=%i %i %i\n', reference.width*10);
    
    % Specify the coordinate frame origin (in mm)
    fprintf(fid, 'Origin=%i %i %i\n', reference.start*10);
    
    % Complete the .mha file header
    fprintf(fid, 'ElementDataFile=LOCAL\n');
    
    % Write the reference image mask data to the temporary file as uint16
    fwrite(fid, refMask, 'ushort', 0, 'l');
    fclose(fid);
    Event(['Reference mask image written to ', refMaskFilename]); 
    
    %% Build daily MHA file
    % Generate a temporary file name for the daily image
    dailyFilename = [tempname, '.mha'];
    
    % Open a write file handle to the temporary daily image
    fid = fopen(dailyFilename, 'w', 'l');
    
    % Start writing the ITK header
    fprintf(fid, 'ObjectType=Image\n');
    fprintf(fid, 'NDims=3\n');
    
    % Specify the dimensions of the merged image
    fprintf(fid, 'DimSize=%i %i %i\n', daily.dimensions);
    
    % Specify the data format (USHORT referring to unsigned 16-bit integer)
    fprintf(fid, 'ElementType=MET_USHORT\n');
    
    % Specify the byte order as little
    fprintf(fid, 'ElementByteOrderMSB=False\n');
    
    % Specify the daily voxel widths (in mm)
    fprintf(fid, 'ElementSize=%i %i %i\n', daily.width*10);
    
    % Specify the daily voxel spacing to equal the widths (in mm)
    fprintf(fid, 'ElementSpacing=%i %i %i\n', daily.width*10);
    
    % Specify the coordinate frame origin (in mm)
    fprintf(fid, 'Origin=%i %i %i\n', daily.start*10);
    
    % Complete the .mha file header
    fprintf(fid, 'ElementDataFile=LOCAL\n');
    
    % Write the merged image data to the temporary file as uint16
    fwrite(fid, daily.data, 'uint16', 0, 'l');
    
    % Close the file handle
    fclose(fid);
    
    % Clear the temporary variable
    clear fid;
    
    % Log where the daily file was saved
    Event(['Daily image written to ', dailyFilename]);
    
    %% Build mask for daily image (excluding outside FOV)
    % Initialize null array of the same size as the daily image
    dailyMask = zeros(daily.dimensions);
    
    % Create meshgrid the same size as one image
    [x,y] = meshgrid(daily.start(1):daily.width(1):...
        daily.start(1) + daily.width(1) * ...
        (daily.dimensions(1) - 1), daily.start(2):...
        daily.width(2):daily.start(2) + ...
        daily.width(2) * (daily.dimensions(2) - 1));
    
    % Set the first mask slice to one within the FOV
    dailyMask(:,:,1) = sqrt(x.^2+y.^2) < dailyFOV/2 - 0.1;
    
    % Loop through each slice
    for i = 2:daily.dimensions(3)
        % Copy the daily mask to each slice
        dailyMask(:,:,i) = dailyMask(:,:,1);
    end
    
    % If the bone flag is enabled
    if bone
        % Update the mask to only include values above 176 HU
        dailyMask = dailyMask .* ceil((daily.data - 1200) / 65535);
    else
        % Otherwise, set mask to exclude image noise (values below -824 HU)
        dailyMask = dailyMask .* ceil((daily.data - 200) / 65535);
    end
    
    % Generate a temporary file name for the daily image mask
    dailyMaskFilename = [tempname, '.mha'];
    
    % Open a write file handle to the temporary daily image mask
    fid = fopen(dailyMaskFilename, 'w', 'l');
    
    % Start writing the ITK header
    fprintf(fid, 'ObjectType=Image\n');
    fprintf(fid, 'NDims=3\n');
    
    % Specify the dimensions of the daily image mask
    fprintf(fid, 'DimSize=%i %i %i\n', daily.dimensions);
    
    % Specify the data format (USHORT referring to unsigned 16-bit integer)
    fprintf(fid, 'ElementType=MET_USHORT\n');
    
    % Specify the byte order as little
    fprintf(fid,'ElementByteOrderMSB=False\n');
    
    % Specify the daily voxel widths (in mm)
    fprintf(fid, 'ElementSize=%i %i %i\n', daily.width*10);
    
    % Specify the daily voxel spacing to equal the widths (in mm)
    fprintf(fid, 'ElementSpacing=%i %i %i\n', daily.width*10);
    
    % Specify the coordinate frame origin (in mm)
    fprintf(fid, 'Origin=%i %i %i\n', daily.start*10);
    
    % Complete the .mha file header
    fprintf(fid, 'ElementDataFile=LOCAL\n');
    
    % Write the merged image mask data to the temporary file as uint16
    fwrite(fid, dailyMask, 'ushort', 0, 'l');
    
    % Close the file handle
    fclose(fid);
    
    % Clear the temporary variable
    clear fid;
    
    % Log where the daily mask file was saved
    Event(['Daily mask image written to ', dailyMaskFilename]);
    
    %% Build plastimatch command file
    % Generate a temporary file name for the command file
    commandFile = [tempname, '.txt'];
    
    % Open a write file handle to the temporary command file
    fid = fopen(commandFile, 'w');
    
    % Specify the inputs to the registration
    fprintf(fid, '[GLOBAL]\n');
    fprintf(fid, 'fixed=%s\n', refFilename);
    fprintf(fid, 'moving=%s\n', dailyFilename);
    fprintf(fid, 'fixed_mask=%s\n', refMaskFilename);
    fprintf(fid, 'moving_mask=%s\n', dailyMaskFilename);
    
    % Generate a temporary filename for the resulting coefficients
    adjustments = [tempname, '.txt'];
    
    % Specify the output file
    fprintf(fid, 'xform_out=%s\n', adjustments);
    
    % If at least two levels are requested
    if levels >= 3
    
        % Specify stage 1 deformable image registration parameters.  Refer to 
        % http://plastimatch.org/registration_command_file_reference.html for
        % more information on these parameters
        fprintf(fid, '[STAGE]\n');
        fprintf(fid, 'xform=align_center\n');
    end
    
    % If at least two levels are requested
    if levels >= 2
        
        % Specify stage 2 parameters
        fprintf(fid, '[STAGE]\n');
        fprintf(fid, 'impl=plastimatch\n');
        fprintf(fid, 'xform=rigid\n');
        fprintf(fid, 'optim=versor\n');
        switch metric
            case 'MI'
                fprintf(fid, 'metric=mi\n');
            otherwise
                fprintf(fid, 'metric=mse\n');
        end
        fprintf(fid, 'max_its=%i\n', iterations);
        fprintf(fid, 'min_step=0.1\n');
        fprintf(fid, 'res=4 4 2\n');
        fprintf(fid, 'threading=cuda\n');
    end

    % Specify stage 3 parameters
    fprintf(fid, '[STAGE]\n');
    fprintf(fid, 'impl=plastimatch\n');
    fprintf(fid, 'xform=rigid\n');
    fprintf(fid, 'optim=versor\n');
    switch metric
        case 'MI'
            fprintf(fid, 'metric=mi\n');
        case 'GM'
            fprintf(fid, 'metric=gm\n');
        otherwise
            fprintf(fid, 'metric=mse\n');
    end
    fprintf(fid, 'max_its=%i\n', iterations);
    fprintf(fid, 'min_step=0.1\n');
    fprintf(fid, 'res=1 1 1\n');
    fprintf(fid, 'threading=cuda\n');
    
    %% Run plastimatch
    % Log execution of system call
    Event(['Executing plastimatch register ', commandFile]);

    % Execute plastimatch using system call, saving the output and status
    [status, cmdout] = system(['plastimatch register ', commandFile]);
    
    % If the status == 0, the command completed successfully
    if status == 0
        
        % Log output
        Event(cmdout);
    else
        % Otherwise, plastimatch didn't complete succesfully, so log the 
        % resulting command output as an error
        Event(cmdout, 'ERROR');
    end
    
    % Clear temporary variables
    clear status cmdout commandFile;
    
    %% Read in registration result
    % Open file handle to temporary file
    fid = fopen(adjustments, 'r');
    
    % Retrieve the first line of the result text
    tline = fgetl(fid);
    
    % Initialize temporary variables to flag if the results are found
    flag1 = 0;
    flag2 = 0;
    
    % Start a while loop to read in result text
    while ischar(tline)
        
        % Search for the text line containing the rigid registration,
        % storing the results and flag if the results are found
        [rigid, flag1] = sscanf(tline, ...
            '\nParameters: %f %f %f %f %f %f\n');
        
        % Search for the text line containing the rigid registration
        % origin, storing the origin and flag if the origin is found
        [origin, flag2] = sscanf(tline, '\nFixedParameters: %f %f %f\n');
        
        % Read in the next line of the results file
        tline = fgetl(fid);
    end
    
    % Close the file handle
    fclose(fid);
    
    % Clear the file handle
    clear fid;
    
    % If both flags are set, the results were successfully found
    if flag1 > 0 && flag2 > 0
        
        % Log an error indicating the the results were not parsed
        % correctly.  This usually indicates the registration failed
        Event(['Unable to parse plastimatch results from ', adjustments], ...
            'ERROR'); 
    else
        % Otherwise, log success
        Event(['Plastimatch results read from ', adjustments]);
    end
    
    % If the registration origin is not equal to the DICOM center
    if ~isequal(origin, [0 0 0])
        
        % Log an error
        Event(['Error: non-zero centers of rotation are not supported', ...
            ' at this time'], 'ERROR'); 
    end
    
    % Clear temporary variables
    clear flag1 flag2 origin;
    
    % Report registration adjustments.  Note angles are stored in radians
    Event(sprintf(['Rigid registration matrix [pitch yaw roll x y z] ', ...
        'computed as [%E %E %E %E %E %E] in %0.3f seconds'], rigid, toc(t)));
     
    % Clear temporary variables
    clear referenceFilename dailyFilename referenceMaskFilename ...
        dailyMaskFilename adjustments commandFile;
    
%% Use MATLAB 6-DOF rigid registration
% This rigid registration technique requires MATLAB's Image Processing
% Toolbox imregtform function.
case 'MATLAB'
    
    % If the bone flag is enabled
    if bone
        
        % Update the fixed image to only include values above 176 HU
        fixed = reference.data .* ...
            ceil((reference.data - 1200) / 65535);
    else
        % Otherwise, set fixed image to exclude image noise (values below 
        % -824 HU)
        fixed = reference.data .* ...
            ceil((reference.data - 200) / 65535);
    end
    
    %% Set fixed (reference) image and reference coordinates
    % Generate a reference meshgrid in the x, y, and z dimensions using the
    % start and width structure fields
    Rfixed = imref3d(reference.dimensions, ...
        [reference.start(1) reference.start(1) + ...
        reference.width(1) * (reference.dimensions(1)-1)], ...
        [reference.start(2) reference.start(2) + ...
        reference.width(2) * (reference.dimensions(2)-1)], ...
        [reference.start(3) reference.start(3) + ...
        reference.width(3) * (reference.dimensions(3)-1)]);
    
    %% Set moving (daily) image and reference coordinates
    % If the bone flag is enabled
    if bone
        
        % Update the moving image to only include values above 176 HU
        moving = daily.data .* ceil((daily.data - 1200) / 65535);
    else
        
        % Otherwise, set moving image to exclude image noise (values below 
        % -824 HU)
        moving = daily.data .* ceil((daily.data - 200) / 65535);
    end
    
    % Generate a reference meshgrid in the x, y, and z dimensions using the
    % start and width structure fields
    Rmoving = imref3d(daily.dimensions, ...
        [daily.start(1) daily.start(1) + ...
        daily.width(1) * (daily.dimensions(1)-1)], ...
        [daily.start(2) daily.start(2) + ...
        daily.width(2) * (daily.dimensions(2)-1)], ...
        [daily.start(3) daily.start(3) + ...
        daily.width(3) * (daily.dimensions(3)-1)]);
    
    %% Run rigid registration
    % Initialize Regular Step Gradient Descent MATLAB object
    optimizer = registration.optimizer.RegularStepGradientDescent();
    
    % Set number of iterations to run
    optimizer.MaximumIterations = iterations;
    
    switch metric
        case 'MI'
            
            % Initialize Mattes Mutual Information metric MATLAB object
            metric = registration.metric.MattesMutualInformation;

            % Log start of optimization
            Event('Executing imregtform rigid using Mattes mutual information');
        otherwise
            % Initialize Mattes Mutual Information metric MATLAB object
            metric = registration.metric.MeanSquares;

            % Log start of optimization
            Event('Executing imregtform rigid using mean square error');
    end
    
    % Execute imregtform using 3 resampling levels
    tform = imregtform(moving, Rmoving, fixed, Rfixed, 'rigid', optimizer, ...
        metric, 'DisplayOptimization', 1, 'PyramidLevels', levels);
    
    % Clear temporary variables
    clear moving fixed Rmoving Rfixed metric optimizer;
    
    % Verify resulting transformation matrix is valid (the values (1,1) and
    % (3,3) must not be zero for atan2 to compute correctly)
    if tform.T(1,1) ~= 0 || tform.T(3,3) ~= 0
        
        % Compute yaw
        rigid(2) = atan2(tform.T(1,2), tform.T(1,1));
        
        % Compute pitch
        rigid(1) = atan2(-tform.T(1,3), ...
            sqrt(tform.T(2,3)^2 + tform.T(3,3)^2));
        
        % Compute roll
        rigid(3) = -atan2(tform.T(2,3), tform.T(3,3));
    else
        % Otherwise, atan2 cannot compute, so throw an error
        Event('Error: incompatible registration matrix determined', ...
            'ERROR');
    end
    
    % Set x, y, and z values
    rigid(4) = tform.T(4,2);
    rigid(5) = tform.T(4,3);
    rigid(6) = tform.T(4,1);
    
    % Clear transformation array
    clear tform;
    
    % Report registration adjustments.  Note angles are stored in radians
    Event(sprintf(['Rigid registration matrix [pitch yaw roll x y z] ', ...
        'computed as [%E %E %E %E %E %E] in %0.3f seconds'], ...
        rigid, toc(t)));

% Otherwise, the method passed to RegisterImages was not supported    
otherwise
    
    % Throw error
    Event(['Unsupported method ', method, ' passed to RegisterImages'], ...
        'ERROR');
end
