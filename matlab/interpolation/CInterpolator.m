classdef (Abstract) CInterpolator < handle
    %% DESCRIPTION
    % Abstract base class for Interpolation classes. 
    % -------------------------------------------------------------------------------------------------------------
    %% CHANGELOG
    % 23-11-2024        Pietro Califano      First version adapted to wrap previous validated functions (codegen ready)      
    % 07-11-2024        Pietro Califano      Fixed error in method to fix discontinuities of signals
    % -------------------------------------------------------------------------------------------------------------
    %% METHODS
    % Method1: Description
    % -------------------------------------------------------------------------------------------------------------
    %% PROPERTIES
    % Property1: Description, dtype, nominal size
    % -------------------------------------------------------------------------------------------------------------
    %% DEPENDENCIES
    % [-]
    % -------------------------------------------------------------------------------------------------------------
properties (SetAccess = protected, GetAccess = public)

        % Settings
        enumInterpType;
        bEnableAutoFitCheck;
        bEnableErrorThrow  (1,1) logical {isscalar} = true
        dPercRelErrorTol   (1,1) double {isscalar} = 0.1
        ui8PolyDeg;
        
        % Data
        dDomainBounds;
        dInterpDomain = [];
        dScaledInterpDomain = [];
        dInterpCoeffsBuffer = [];

        % Auxiliary
        strFitStats;
        i32OutputVectorSize;

        % Quaternion specific
        bIsSignSwitched
        ui8HowManySwitches
        bSignSwitchDetectionMask
        dSwitchIntervals

    end

    methods (Access = public)
        %% CONSTRUCTOR
        function self = CInterpolator(dInterpDomain, ...
                                    ui8PolyDeg, ...
                                    enumInterpType, ...
                                    bEnableAutoFitCheck, ...
                                    dDomainBounds, ...
                                    i32OutputVectorSize, ...
                                    bEnableErrorThrow, ...
                                    dPercRelErrorTol )
            arguments
                dInterpDomain       (1,:) double {mustBeNumeric, isvector}
                ui8PolyDeg          (1,1) uint8 {mustBeNumeric, isscalar} = 15
                enumInterpType      (1,1) {isa(enumInterpType, 'EnumInterpType')} = EnumInterpType.VECTOR
                bEnableAutoFitCheck (1,1) logical {islogical} = true
                dDomainBounds       (1,2) double {mustBeNumeric, isvector} = zeros(1,2)
                i32OutputVectorSize (1,1) int32 {mustBeNumeric, isscalar} = -1 % Expected output size for checks
                bEnableErrorThrow   (1,1) logical {islogical, isscalar} = true
                dPercRelErrorTol    (1,1) double {mustBeNumeric, isscalar}  = 0.1
            end
            
            assert(ui8PolyDeg > 1, "MATLAB:assert:failed", 'Interpolant degree must be > 1.')
            assert(length(dInterpDomain) > ui8PolyDeg, 'Size of interpolation domanin and data must be greater than the requested ui8PolyDeg.')

            % Save class data members
            self.enumInterpType         = enumInterpType;
            self.ui8PolyDeg             = ui8PolyDeg;
            self.dInterpDomain          = dInterpDomain;
            self.bEnableAutoFitCheck    = bEnableAutoFitCheck;
            self.i32OutputVectorSize    = i32OutputVectorSize;
            self.bEnableErrorThrow      = bEnableErrorThrow;
            self.dPercRelErrorTol       = dPercRelErrorTol;

            if i32OutputVectorSize == -1
                warning('Expected output size not provided. Interpolator will get it from data matrix, and will not perform validation checks.')
            end

            % If not provided as input, compute LB and UB from dInterpDomain
            if all(dDomainBounds == [0,0]) == true
                dDomainBounds = [min(dInterpDomain, [], 'all'), max(dInterpDomain, [], 'all')];
            end

            assert(dDomainBounds(2) > dDomainBounds(1), 'Interpolation domain UB must be > than LB!')
            self.dDomainBounds = dDomainBounds;

        end

        % GETTERS

        % SETTERS

        % METHODS

    end

    methods (Access = protected)
        function [self, dModifiedDataMatrix, dSwitchIntervals, bIsSignSwitched, ...
                ui8HowManySwitches ] = fixQuatSignDiscontinuity(self, dQuatMatrix_fromAtoB) %#codegen
            arguments
                self
                dQuatMatrix_fromAtoB (:, 4) double {mustBeNumeric, ismatrix}
            end

            % Sign discontinuity detection and fix
            self.ui8HowManySwitches = uint8(0);

            % Initialize output variables
            ui32NumOfTimes      = uint32(size(dQuatMatrix_fromAtoB, 1));

            self.bSignSwitchDetectionMask = false(ui32NumOfTimes, 1);
            self.bIsSignSwitched          = false(ui32NumOfTimes, 1);

            bIsNewJump = true;
            for idT = 1:ui32NumOfTimes

                if idT > 1

                    % Check sign change of the max component of the quaternion
                    [~, ui32MaxCompIdx] = max( abs( dQuatMatrix_fromAtoB(idT,:) ) );

                    % Check if product of previous timestamp and current has negative sign
                    if dQuatMatrix_fromAtoB(idT, ui32MaxCompIdx) * dQuatMatrix_fromAtoB(idT-1,ui32MaxCompIdx) < 0
                        % Jump detected: switch sign
                        dQuatMatrix_fromAtoB(idT,:) = - dQuatMatrix_fromAtoB(idT,:);

                        self.bIsSignSwitched(idT) = true;

                        if bIsNewJump
                            self.bSignSwitchDetectionMask(idT) = true;
                            self.ui8HowManySwitches = self.ui8HowManySwitches + 1;
                            bIsNewJump = false;
                        end
                    else
                        bIsNewJump = true;
                    end

                end

            end

            % Return fixed data matrix
            dModifiedDataMatrix = dQuatMatrix_fromAtoB';

            assert(sum(all(ischange(sign(dQuatMatrix_fromAtoB)), 2) == true) == 0, 'Something may have gone wrong in fixing the discontinuity!')
            assert(self.ui8HowManySwitches <= 255, 'SAFETY STOP: possible overflow in self.ui8howManySwitches due to presence of >255 switches!')

            % Determine sign switch intervals and store in data members
            self.dSwitchIntervals = zeros(self.ui8HowManySwitches, 2, 'double');

            if self.ui8HowManySwitches > 0
                dSwitchesIDs = find(self.bSignSwitchDetectionMask, self.ui8HowManySwitches);

                self.dSwitchIntervals(:, 1) = dSwitchesIDs;
                ui32HowManySamples = length(self.bIsSignSwitched);

                for idC = 1:self.ui8HowManySwitches
                    idTmp = self.dSwitchIntervals(idC);
                    while idTmp < ui32HowManySamples && self.bIsSignSwitched(idTmp) == 1
                        idTmp = idTmp + 1;
                    end
                    self.dSwitchIntervals(idC, 2) = idTmp;
                end
            end

            dSwitchIntervals   = self.dSwitchIntervals;
            bIsSignSwitched    = self.bIsSignSwitched;
            ui8HowManySwitches = self.ui8HowManySwitches;
        end
    end

    methods (Access = public)
        % Fitting check method
        function [self, strFitStats] = checkFitPoly(self, ...
                                                   dEvalDomain, ...
                                                   dDataMatrix, ...
                                                   bEnableErrorThrow, ...
                                                   dPercRelErrorTol)
            arguments
                self,
                dEvalDomain (1,:) double {mustBeNumeric}
                dDataMatrix (:,:) double {ismatrix, mustBeNumeric}
                bEnableErrorThrow   (1,1) logical {islogical, isscalar} = self.bEnableErrorThrow
                dPercRelErrorTol    (1,1) double {mustBeNumeric, isscalar}  = self.dPercRelErrorTol
            end
            
            % HARDCODED OPTIONS
            ui32Npoints = min(1000, size(dDataMatrix,2)); % Select number of points based on input size
            ui32Npoints = ui32Npoints - 2;
            
            %% Function code
            bIS_ATT_QUAT = false;
            if self.enumInterpType == EnumInterpType.QUAT
                bIS_ATT_QUAT = true;
            end

            % Determine 
            if bIS_ATT_QUAT
                i32OutputSize = int32(4); % HARDCODED for specialization
                assert( size(dDataMatrix, 1) == i32OutputSize );
            elseif self.i32OutputVectorSize ~= -1
                % Use output size set at instantiation
                i32OutputSize = self.i32OutputVectorSize;
            else
                i32OutputSize = int32(size(dDataMatrix, 1));
            end

            % Case checks
            assert(i32OutputSize > 0)
            assert( size(dDataMatrix, 2) == length(dEvalDomain) );

            % DEVNOTE TODO: Check if evaluation point (in dEvalDomain fall outside self.dInterpDomain)
            % warning('TODO: add check on dEvalDomain vs self.dInterpDomain')

            % Evaluation at test points
            ui32testpointsIDs = uint32(sort( randi( length(dEvalDomain), ui32Npoints, 1 ), 'ascend' ));
            
            % Get points
            dTestPoints_Time = [dEvalDomain(1); dEvalDomain(ui32testpointsIDs)'; dEvalDomain(end)];
            dTestPoints_Labels = [dDataMatrix(:, 1),...
                                dDataMatrix(:, ui32testpointsIDs), ...
                                dDataMatrix(:, end)];
                            
            % Evaluate interpolant
            dChbvInterpVector = zeros(i32OutputSize, length(dTestPoints_Time));
            evalRunTime = zeros(length(dTestPoints_Time), 1);

            for idP = 1:length(dTestPoints_Time)
                dEvalPoint = dTestPoints_Time(idP);
                
                % Evaluate interpolant
                tic;

                [self, dChbvInterpVector(:, idP)] = self.evalInterpolant(dEvalPoint, true);
                
                evalRunTime(idP) = toc;
            end
            fprintf("\nAverage interpolant evaluation time: %4.4g [s]\n", mean(evalRunTime))

            % Error evaluation
            strFitStats = struct();

            if self.enumInterpType == EnumInterpType.QUAT
                strFitStats.dAbsErrVec = abs(abs(dot(dChbvInterpVector, dTestPoints_Labels, 1)) - 1);

                strFitStats.dMaxAbsErr = max(strFitStats.dAbsErrVec, [], 'all');
                strFitStats.dAvgAbsErr = mean(strFitStats.dAbsErrVec, 2);

                fprintf('Max absolute difference of (q1-dot-q2 - 1): %4.4g [-]\n', strFitStats.dMaxAbsErr);
                fprintf('Average absolute difference of (q1-dot-q2 - 1): %4.4g [-]\n', strFitStats.dAvgAbsErr);

                if bEnableErrorThrow
                    assert( all([strFitStats.dAvgAbsErr, strFitStats.dMaxAbsErr] < 1E-3), ...
                        ['ERROR: fitting validation failed to meet tolerances for attitude quaternion. ' ...
                        'Found distance greater than 1E-3 at sampling nodes.']);
                end

            else
                strFitStats.dAbsErrVec = abs(dChbvInterpVector - dTestPoints_Labels);
                strFitStats.dRelErrVec = strFitStats.dAbsErrVec./vecnorm(dTestPoints_Labels, 2, 1);

                strFitStats.dMaxAbsErr = max(strFitStats.dAbsErrVec, [], 'all');
                strFitStats.dAvgAbsErr = mean(strFitStats.dAbsErrVec, 2);

                strFitStats.dMaxRelErr = 100*max(strFitStats.dRelErrVec, [], 'all');
                strFitStats.dAvgRelErr = 100*mean(strFitStats.dRelErrVec, 2);

                % Printing
                fprintf('Max absolute error: %4.4g [-]\n', strFitStats.dMaxAbsErr);
                fprintf('Average absolute error: %4.4g, %4.4g, %4.4g [-]\n', strFitStats.dAvgAbsErr(1), ...
                    strFitStats.dAvgAbsErr(2), strFitStats.dAvgAbsErr(3));

                fprintf('\nMax relative error: %4.4g [%%]\n', strFitStats.dMaxRelErr);
                fprintf('Average relative error: %4.4g, %4.4g, %4.4g [%%]\n', strFitStats.dAvgRelErr(1), ...
                    strFitStats.dAvgRelErr(2), strFitStats.dAvgRelErr(3));

                if bEnableErrorThrow
                    assert( all([strFitStats.dMaxRelErr, max(strFitStats.dAvgRelErr)] <= dPercRelErrorTol), ...
                        ['ERROR: fitting validation failed to meet relative tolerances. ' ...
                        'Found relative error greater than 0.1% at sampling nodes.']);
                end
            end

        end

end

    % Abstract methods
    methods (Abstract, Access = public)
              
        % Interpolant evaluation method
        [self, dInterpVector] = evalInterpolant(self, dEvalPoint, i16LimitDegree);

        % Interpolant terms evaluation
        [self, dPolyTermsValues] = evalPoly(self, dEvalPoint, i16LimitDegree);

        % Data matrix fitting method
        [self, dInterpCoeffsMatrix, strFitStats] = fitDataMatrix(self, dDataMatrix);
        
    end

end
