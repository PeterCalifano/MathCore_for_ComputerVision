close all
clear
clc

LoadUserPathConfig; % Does nothing for now
SetupOptions;
charProjectDir = pwd;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SCENARIO SETUP
%% Options and input specification
enumTrajName                = EnumTrajectoryNames.SSTO_1p4; %EnumTrajectoryNames.RTO_4t1_J13p0; % RTO_4t1_J11p0_60dt, SSTO_14_60dt, SSTO_12_60dt

% Timegrid
ui32InitialTimeID           = 1;
dFramesDeltaTime            = 120;    % [s]
dTotalTimeDuration          = 12*86400;  % [s]
ui32ImgAcquisitionStepFreq  = uint32(1);
dInitialRelTimestamp = 0; 
bUseKernelInitialTimestamp = true;
bIsTimeGridRelative = true;

% Output units
bConvertData_m2mBU_beforeRender = true;
charKernelLengthUnits           = "km";
charOutputLengthUnits           = "m";

% Frames
enumWorldFrame                  = "J2000";
enumTargetFrame                 = "APOPHIS_FIXED_HF"; %EnumFrameName_RCS1.APOPHIS_FIXED; % "IAU_EARTH"
% enumTargetFrame                 = EnumFrameName_RCS1.APOPHIS_FIXED; % "IAU_EARTH"

%%% Trajectory kernel loader input specifications
% Target point, body and frames
ui32RenderingTargetBodyID       = uint32(0); % Keep 0 if target is the main body
varTargetID                     = int32(-19920605); % -10003001
% Target frames 
varReferenceCentre          = "APOPHIS";
varTargetBodyID             = "APOPHIS";

% Kernel name and mk path
enumTrajectKernelName       = enumTrajName; % kernelFUT_h600

% charTrajKernelFolderPath    = "/home/peterc/devDir/projects-DART/data/future/phase-C/kernels/trajectories";
charTrajKernelFolderPath = "/home/peterc/devDir/projects-DART/data/rcs-1/phase-B/kernels/trajectories";

% Timescale for timegrids from kernel
charKernelTimescale  = "ET";
charUserDefTimescale = "ET";

% Additional bodies
% cellAdditionalTargetsID          = {"MOON"};
% cellAdditionalTargetNames        = {};
% bRequireAdditionalBodiesAttitude = true(1,1);
% cellAdditionalTargetsFrames      = {"IAU_MOON"};

cellAdditionalTargetsID          = {};
cellAdditionalTargetNames        = {};
bRequireAdditionalBodiesAttitude = false(1,1);
cellAdditionalTargetsFrames      = {};

dTimegrid = 0.0:120:8*86400;
bIsTimeGridRelative = true;

charKernelPathRoot = "/home/peterc/devDir/projects-DART/data/rcs-1/phase-B/kernels/mk";

if not(strcmpi(charKernelPathRoot, ""))
    charCurrentDir = pwd;
    cd(charKernelPathRoot)
    cspice_furnsh('kernels.mk')
    cd(charCurrentDir);
end

% Load dataset data from kernels
[objDataset, ~, ~] = LoadReferenceDataFromKernels(varTargetID, ...
                                                enumTrajectKernelName, ...
                                                dTimegrid, ...
                                                enumWorldFrame, ...
                                                varReferenceCentre, ...
                                                enumTargetFrame, ...
                                                "bLoadManoeuvres", false, ...
                                                "cellAdditionalTargetsID", cellAdditionalTargetsID, ...
                                                "cellAdditionalTargetNames", cellAdditionalTargetNames, ...
                                                "charKernelLengthUnits", charKernelLengthUnits,...
                                                "charTrajKernelFolderPath", charTrajKernelFolderPath, ...
                                                "charOutputLengthUnits", charOutputLengthUnits, ...
                                                "varTargetBodyID", varTargetBodyID, ...
                                                "bIsTimeGridRelative", bIsTimeGridRelative, ...
                                                "bAdditionalBodiesRequireAttitude", bRequireAdditionalBodiesAttitude, ...
                                                "cellAdditionalTargetFrames", cellAdditionalTargetsFrames, ...
                                                "charKernelTimescale", charKernelTimescale, ...
                                                "charUserDefTimescale", charUserDefTimescale, ...
                                                "bUseKernelInitialTimestamp", bUseKernelInitialTimestamp);

% Generate new attitude data for the target, starting from initial condition
dAbsIntergTimegrid = objDataset.dTimestamps;
dQuat0      = DCM2quat(objDataset.dDCM_TBfromW(:,:,1), false);
% i32CKHandle = cspice_cklpf( char(charCKKernelName) );

% Documentation: https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/MATLAB/mice/cspice_ckgpav.html
% i32ApophisFrame_ID = -200999422; % int32 ID of the frame 
% i32ApophisFrame_ID = 20099942;
% [dDCM, dRefAngVelSeq, dTimeOuts, bFound] = cspice_ckgpav(int32(i32CKHandle), ...
%                                                            dAbsIntergTimegrid, ...
%                                                            1E-6, ...
%                                                            'APOPHIS_FIXED_HF');

% Get 6x6 rotation matrices
dSxDCM_TFfromW = cspice_sxform('J2000', 'APOPHIS_FIXED_HF', dAbsIntergTimegrid );

% Convert into DCM and angular velocity
[ dDCM_GT_TFfromW,  dRefAngVelSeq_TFfromW ] = cspice_xf2rav( dSxDCM_TFfromW );

%% Integrate using constant angular velocity
% Compute perturbed angular velocity at initial time instant
dPertubAngVel0 = dRefAngVelSeq_TFfromW(:,1);

% Integrate attitude kinematics 
objQuatIntegr = CQuatKinematicsIntegrator();
dTimestep = 1;

[dQuatEndConstant0, dTimegridOut, dQuatSeqConstant0] = objQuatIntegr.integrate(dQuat0, ...
                                                                dPertubAngVel0, ...
                                                                dTimegrid, ...
                                                                'rkmk4', ...
                                                                dTimestep);
% Plot visualization
objDataset.plotDatasetData("dTargetAttitudeSet2", QuatSeq2DCM(dQuatSeqConstant0, false), ...
                            "bPlotTargetAttitude", true, ...
                            "charLblTargetAttitudeSet2", "Integrated", ...
                            "charTargetAttitudePlotTitle", 'Ang. vel. fixed as initial time');
return

%% Integrate using angular velocity profile (internally defined spline)
% Integrate attitude kinematics 
objQuatIntegr = CQuatKinematicsIntegrator(); %#ok<UNRCH>
dTimestep = 5;

tic
[dQuatEndConstant0, dTimegridOut, dQuatSeqConstant0] = objQuatIntegr.integrate(dQuat0, ...
                                                                dRefAngVelSeq_TFfromW, ...
                                                                dTimegrid, ...
                                                                'rkmk4', ...
                                                                dTimestep, ...
                                                                5.0, ...
                                                                dAbsIntergTimegrid, ...
                                                                "linear");
toc

% Plot visualization
objDataset.plotDatasetData("dTargetAttitudeSet2", QuatSeq2DCM(dQuatSeqConstant0, false), ...
                            "bPlotTargetAttitude", true, ...
                            "charLblTargetAttitudeSet2", "Integrated", ...
                            "charTargetAttitudePlotTitle", 'Ang. vel. spline as reference');
return
%% test_TorqueFreeMotion_SymmetricTop
% Torque-free motion for symmetric top (I1=I2 != I3)
dI1 = 10;
dI2 = 10;
dI3 = 50;

dI = diag([dI1, dI2, dI3]);
dPertubAngVel0 = dRefAngVelSeq_TFfromW(:,1);

objIntegrator = CRigidBodyDynamicsIntegrator(dI, CQuatKinematicsIntegrator());

% Time grid
dTimegrid = 0:0.1:100.0;

% Initial angular velocity: small transverse and dominant spin about symmetry axis
dOmega0 = [0.05; 0.0; 2.0];
dQuat0 = [1; 0; 0; 0];

varTorque = zeros(3,1);
[dQuatSeq, dOmegaSeq] = objIntegrator.integrate(dTimegrid, ...
    dQuat0, dOmega0, varTorque, 0.1, 'rk4_rkmk4', true, 1.0);

% Analytical precession frequency: dOmega_p = (I3 - I1)/I1 * omega3
dOmega_p = (dI3 - dI1)/dI1 * dOmega0(3);
dA = dOmega0(1);

% Visualize attitude sequence
dOrigin_Frame = zeros(3,1);

cellPlotColors = {'r', 'g', 'b'};
cellPlotNames  = {'X', 'Y', 'Z'};

objFig = figure('Renderer', 'opengl');

[~, charTextColor, ~] = DefaultPlotOpts(objFig, ...
    "charRenderer", "opengl", ...
    "bUseBlackBackground", true);

xlabel('X [-]');
ylabel('Y [-]');
zlabel('Z [-]');
title('Propagated target attitude frames');
grid off

ui32Decimation = 1;
for idAtt = 1:ui32Decimation:size(dQuatSeq, 2)

    % Convert quaternion to DCM
    dCamDCM_RenderFrameFromCam = Quat2DCM(dQuatSeq(:,idAtt), false);

    [cellCameraAxes] = PlotFrameFromDCM(dOrigin_Frame, ...
                                    dCamDCM_RenderFrameFromCam, ...
                                    cellPlotColors, ...
                                    cellPlotNames, ...
                                    objFig, ...
                                    "dAxisScale", 1.5, ...
                                    "bShowArrowHead", true);
    axis([-2 2 -2 2 -2 2]);
    view(45,45);
    drawnow
    pause(0.05)
    for idB = 1:length(cellCameraAxes)
        delete(cellCameraAxes{idB})
    end
end
return
%% test_TorqueFreeMotion_integrated


