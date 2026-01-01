function [cellShellList, cellNodesArray, ui32CornersIDs] = IcosahedronSampleS2Sphere(ui32DensityArray, dSphereRadii, charSamplingRoutinePath)
arguments
    ui32DensityArray
    dSphereRadii
    charSamplingRoutinePath
end
%% PROTOTYPE
% [cellShellList, cellNodesArray, ui32CornersIDs] = IcosahedronSampleS2Sphere(ui32DensityArray, dSphereRadii, charSamplingRoutinePath)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% This function is a MATLAB wrapper of the "Icosahedron" FORTRAN routine 
% for Sphere pixelization. It produces a uniformly spaced discretization of
% the sphere, the number of points being dependent on the request density.
% (See code and reference for more details).
% REFERENCE:
% 1) Max Tegmark, Max-Planck-Institut fuer Physik, Munich, April 1996
%
% NOTE: Recompile the .f source if program execution returns missing gcc
% library error. Compiler gfortran is required.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% densityArray: [N,1]   Array containing the integer values of the density
%                       for each requested discretization. If not integer,
%                       round is applied.
% sphereRadii:  [N,1]   Array containing the radius of the sphere to scale
%                       the 3D coordinates of the grid points.
% path2Icomain: [1]     String specifying the path to icomain.exe. 
%                       Default: retrieved using "which" function, assuming
%                       it is on the MATLAB path
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% shellList:      [N,1]     String array containing the names of each shell
% nodesArrayCell: [1,N]     Cell array containing the N discretized grids
%                           as [Npoints, 4] arrays: [id,x,y,z] entries
% cornersIDs:     [N,12]    2D array specifying the ids of the corners of
%                           the icosahedron (having only 5 adjacent points)
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 18-01-2024      Pietro Califano    Adapted from work on AOCS Path Planner
%                                    for IAC2023 paper (Perico et al.)
% 04-04-2024      Pietro Califano    Linux x86_64 version added.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% icomain.exe function compiled from icosahedron.f (using gfortran). By
% default it is in the same folder as the function.
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code

% Number of nodes function
handleComputeNumNodes = @(densityValue) 40*densityValue*(densityValue - 1) + 12; 

%% ID map generation
% Generation of the ID maps by call to exe with n_density vector specifying
% required node density. Cycle only computes maps with densities not found
% in nodes_maps folder.

% TODO: update function

if strcmp(computer, 'PCWIN64')
    charProgramName = which("icomain.exe");

elseif strcmp(computer, 'GLNXA64')
    charProgramName = "icomain";
end

if not(exist('path2Icomain', 'var'))
    charSamplingRoutinePath = fullfile(charProgramName);
else
    charSamplingRoutinePath = fullfile(charSamplingRoutinePath, charProgramName);
end

% Static allocation
Nshells = length(ui32DensityArray);
filename = strings(Nshells, 1);
cellShellList = strings(Nshells, 1);
Nnodes = zeros(Nshells, 1);
ui32CornersIDs = zeros(Nshells, 12);

cellNodesArray = cell(1, Nshells);

%% FOR WINDOWS
% BAT FILE WRITING and CALL
for shellID = 1:Nshells

    if strcmp(computer, 'PCWIN64')
        batchfile_code = cell(2);

        % Write number of nodes of the ith layer and name
        Nnodes(shellID) = handleComputeNumNodes(round(ui32DensityArray(shellID)));
        fprintf('\nNumber of nodes = %4.0f', Nnodes(shellID));
        NewShellName = "Shell" + num2str(shellID);
        cellShellList(shellID) = NewShellName;

        filename(shellID) = "nodesArray_density" + num2str(ui32DensityArray(shellID)) + ".dat";

        % Check if nodedensity map already exist or not
        if ~exist(filename(shellID) , "file")
            %   Write batch file to call icomain
            command_string = "ECHO   "   + num2str(ui32DensityArray(shellID) )  + "  | " + '"' + charSamplingRoutinePath + '"';
            batchfile_code{1} = "ECHO CALL ICOMAIN.EXE";
            batchfile_code{2} = command_string;
            filetowrite = fopen('call_icomain.bat', 'w');
            for rowid = 1:length(batchfile_code)
                fprintf(filetowrite,'%s \n', batchfile_code{rowid});
            end
            fclose('all');

            % Create ico_ID_map.dat file by calling exe
            system('call_icomain.bat');
        end


    elseif strcmp(computer, 'GLNXA64')

        commandString = strcat("echo", num2str(ui32DensityArray(shellID) ), "| ./icomain");
        system(commandString);

    else
        error('OS architecture not supported.')
    end

    % Rename ico_ID_map.dat as nodedensity<n_density>.dat and move to nodes_maps folder
    name = "ico_ID_map.dat";
    movefile(name, filename(shellID));

    formatSpec = '%6f%6f%9f%9f%f%[^\n\r]';

    fileID = fopen(filename(shellID),'r');

    % Read columns of data according to the format.
    % This call is based on the structure of the file used to generate this code. If an error occurs for a different file, try regenerating the code from the Import Tool.
    dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'EmptyValue', NaN,  'ReturnOnError', false);
    fclose(fileID);

    NodesTable = table(dataArray{1:end-1}, 'VariableNames', {'VarName1','VarName2','VarName3','VarName4','VarName5'});

    % Columns 3, 4, 5 corresponds to the node coordinates X, Y, Z
    % ACHTUNG: The numbering of the nodes starts from 0
    NodesArrayTmp = table2array(NodesTable);
    % Nodes_array cannot be 3D array --> structure with variable field name
    % taken from shell_lit

    % ID shift to start indexing from 1 instead of 0
    NodesArrayTmp(:, 1:2) = NodesArrayTmp(:, 1:2) + 1;

    % DEVNOTE: check what are the first and second columns and get what is
    % needed only.

    % Apply sphere scaling (multiplication by radius) and assign discretization nodes to cell
    cellNodesArray{shellID} = dSphereRadii(shellID) * NodesArrayTmp(:, 2:end);

    % Create corners ids
    nCorners = 0:11;
    ui32CornersIDs(shellID, :) = (Nnodes(shellID) - nCorners);

    % Remove temporary files
    delete(filename(shellID))
    delete("call_icomain.bat")
end

fprintf('\n');


end
