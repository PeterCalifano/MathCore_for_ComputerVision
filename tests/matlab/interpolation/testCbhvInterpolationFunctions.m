classdef testCbhvInterpolationFunctions < matlab.unittest.TestCase
    %% DESCRIPTION
    % Unit test class for Chebyshev polynomial interpolation functions.
    % Tests fitting, evaluation, and quaternion-specific functionality.
    % -------------------------------------------------------------------------------------------------------------
    %% CHANGELOG
    % 18-07-2025   Claude Sonnet 4, Pietro Califano   Comprehensive test suite for Chebyshev interpolation functions
    % -------------------------------------------------------------------------------------------------------------
    
    properties (Access = private)
        charSavedPath   % char vector of original path
        charPathFile
    end

    methods (TestClassSetup)
        function setupPath(testCase)
            % Save current path (memory + file, if you want a record)
            testCase.charSavedPath = path;
            testCase.charPathFile  = fullfile(fileparts(mfilename("fullpath")), 'path_before_tests.mat');
            charOrigPath = testCase.charSavedPath;
            save(testCase.charPathFile, 'charOrigPath');

            restoredefaultpath;

            % Add the folder you want to test (plus subs if needed)
            charThisFilePath = fileparts(mfilename("fullpath"));
            charSrcPath = fullfile(charThisFilePath, '..', '..', '..','matlab');
            addpath(genpath(charSrcPath), '-begin');             % or addpath(genpath(p), '-begin');
            testCase.addTeardown(@() teardownPath(testCase));
        end
    end

    methods (Access = private)
        function teardownPath(testCase)
            % Restore original path
            if ~isempty(testCase.charSavedPath)
                path(testCase.charSavedPath);
            end

            % Clean up temp file
            if ~isempty(testCase.charPathFile) && isfile(testCase.charPathFile)
                delete(testCase.charPathFile);
            end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% TEST
    methods (Test)
        
        function testFitChbvPolynomialsConstant(testCase)
            %% Test fitting of constant function f(x) = 2
            dInterpDomain = linspace(0, 1, 10)';
            dDataMatrix = 2 * ones(1, length(dInterpDomain));
            ui32PolyDeg = uint32(3);
            dDomainLB = 0;
            dDomainUB = 1;
            
            [dChbvCoeffs, dScaledInterpDomain, strFitStats] = ...
                fitChbvPolynomials(ui32PolyDeg, dInterpDomain, dDataMatrix, ...
                                 dDomainLB, dDomainUB, true);
            
            % Check outputs
            testCase.verifySize(dChbvCoeffs, [1*(ui32PolyDeg+1), 1]);
            testCase.verifySize(dScaledInterpDomain, size(dInterpDomain));
            testCase.verifyTrue(isstruct(strFitStats));
            
            % For constant function, first coefficient should be ~2, others ~0
            coeffsMatrix = reshape(dChbvCoeffs, [ui32PolyDeg+1, 1]);
            testCase.verifyEqual(coeffsMatrix(1), 2, 'AbsTol', 1e-10, ...
                'First coefficient should equal constant value');
            testCase.verifyEqual(coeffsMatrix(2:end), zeros(ui32PolyDeg, 1), ...
                'AbsTol', 1e-10, 'Higher order coefficients should be zero');
        end
        
        function testFitChbvPolynomialsLinear(testCase)
            %% Test fitting of linear function f(x) = 3*x + 1
            dInterpDomain = linspace(0, 2, 15)';
            dDataMatrix = 3 * dInterpDomain' + 1;
            ui32PolyDeg = uint32(2);
            dDomainLB = 0;
            dDomainUB = 2;
            
            [dChbvCoeffs, ~, ~] = fitChbvPolynomials(ui32PolyDeg, dInterpDomain, ...
                                                   dDataMatrix, dDomainLB, dDomainUB, false);
            
            % Verify interpolation accuracy at test points
            testPoints = [0.5, 1.0, 1.5];
            for i = 1:length(testPoints)
                expectedValue = 3 * testPoints(i) + 1;
                interpolatedValue = evalChbvPolyWithCoeffs(ui32PolyDeg, uint32(1), ...
                    testPoints(i), dChbvCoeffs, dDomainLB, dDomainUB, ...
                    uint32(length(dChbvCoeffs)), ui32PolyDeg);
                
                testCase.verifyEqual(interpolatedValue, expectedValue, 'AbsTol', 1e-12, ...
                    sprintf('Linear interpolation failed at x=%.1f', testPoints(i)));
            end
        end
        
        function testFitChbvPolynomialsMultiDimensional(testCase)
            %% Test fitting of multi-dimensional vector function
            dInterpDomain = linspace(-1, 1, 20)';
            % 3D vector function: [sin(x), cos(x), x^2]
            dDataMatrix = [sin(dInterpDomain)'; cos(dInterpDomain)'; dInterpDomain'.^2];
            ui32PolyDeg = uint32(6);
            dDomainLB = -1;
            dDomainUB = 1;
            
            [dChbvCoeffs, ~, strFitStats] = ...
                fitChbvPolynomials(ui32PolyDeg, dInterpDomain, dDataMatrix, ...
                                 dDomainLB, dDomainUB, true);
            
            % Check dimensions
            expectedCoeffsSize = size(dDataMatrix, 1) * (ui32PolyDeg + 1);
            testCase.verifySize(dChbvCoeffs, [expectedCoeffsSize, 1]);
            testCase.verifyTrue(isstruct(strFitStats));
            
            % Test interpolation at a few points
            testPoints = [-0.5, 0, 0.7];
            for i = 1:length(testPoints)
                x = testPoints(i);
                expectedVector = [sin(x); cos(x); x^2];
                
                interpolatedVector = evalChbvPolyWithCoeffs(ui32PolyDeg, uint32(3), ...
                    x, dChbvCoeffs, dDomainLB, dDomainUB, ...
                    uint32(length(dChbvCoeffs)), ui32PolyDeg);
                
                testCase.verifyEqual(interpolatedVector, expectedVector, 'AbsTol', 1e-3, ...
                    sprintf('Multi-dimensional interpolation failed at x=%.1f', x));
            end
        end
        
        function testEvalChbvPolyWithCoeffsEdgeCases(testCase)
            %% Test evaluation at domain boundaries and scaling
            dInterpDomain = linspace(2, 8, 10)';
            dDataMatrix = dInterpDomain' + 5; % f(x) = x + 5
            ui32PolyDeg = uint32(2);
            dDomainLB = 2;
            dDomainUB = 8;
            
            [dChbvCoeffs, ~, ~] = fitChbvPolynomials(ui32PolyDeg, dInterpDomain, ...
                                                   dDataMatrix, dDomainLB, dDomainUB, false);
            
            % Test at domain boundaries
            testPoints = [dDomainLB, dDomainUB];
            for i = 1:length(testPoints)
                x = testPoints(i);
                expectedValue = x + 5;
                
                interpolatedValue = evalChbvPolyWithCoeffs(ui32PolyDeg, uint32(1), ...
                    x, dChbvCoeffs, dDomainLB, dDomainUB, ...
                    uint32(length(dChbvCoeffs)), ui32PolyDeg);
                
                testCase.verifyEqual(interpolatedValue, expectedValue, 'AbsTol', 1e-10, ...
                    sprintf('Boundary evaluation failed at x=%.1f', x));
            end
        end
        
        function testFixQuatSignDiscontinuityNoSwitch(testCase)
            %% Test quaternion discontinuity fix with no switches
            dQuatSequence = [ones(50, 1), 0.5*ones(50, 1), zeros(50, 2)];
            % Normalize quaternions
            for i = 1:size(dQuatSequence, 1)
                dQuatSequence(i, :) = dQuatSequence(i, :) / norm(dQuatSequence(i, :));
            end
            
            [dDataMatrix, bIsSignSwitched, ui8HowManySwitches, bSignSwitchMask] = ...
                fixQuatSignDiscontinuity(dQuatSequence);
            
            % Should be transposed but unchanged
            testCase.verifyEqual(dDataMatrix, dQuatSequence', 'AbsTol', 1e-15);
            testCase.verifyEqual(ui8HowManySwitches, uint8(0));
            testCase.verifyFalse(any(bIsSignSwitched));
            testCase.verifyFalse(any(bSignSwitchMask));
        end
        
        function testFixQuatSignDiscontinuityWithSwitch(testCase)
            %% Test quaternion discontinuity fix with sign changes
            % Create quaternion sequence with deliberate sign flip
            t = linspace(0, 1, 100);
            axis = [1; 0; 0];
            angle = pi/4 * t;
            
            dQuatSequence = zeros(length(t), 4);
            for i = 1:length(t)
                dQuatSequence(i, :) = [cos(angle(i)/2), sin(angle(i)/2)*axis'];
            end
            
            % Introduce sign flip in the middle
            switchIdx = 50;
            dQuatSequence(switchIdx:end, :) = -dQuatSequence(switchIdx:end, :);
            
            [dDataMatrix, bIsSignSwitched, ui8HowManySwitches, ~] = ...
                fixQuatSignDiscontinuity(dQuatSequence);
            
            testCase.verifyEqual(ui8HowManySwitches, uint8(1), ...
                'Should detect exactly one sign switch');
            testCase.verifyTrue(any(bIsSignSwitched), ...
                'Should flag some quaternions as sign-switched');
            
            % Verify no more discontinuities remain
            for i = 2:size(dDataMatrix, 2)
                dotProduct = abs(dot(dDataMatrix(:, i), dDataMatrix(:, i-1)));
                testCase.verifyGreaterThan(dotProduct, 0.9, ...
                    'Adjacent quaternions should have high correlation after fixing');
            end
        end
        
        function testEvalAttQuatChbvPolyWithCoeffs(testCase)
            %% Test attitude quaternion-specific polynomial evaluation
            % Create smooth quaternion trajectory
            t = linspace(0, 1, 20)';
            axis = [0; 0; 1]; % Rotation about z-axis
            angles = pi/2 * t; % 90 degree rotation over time
            
            dQuatSequence = zeros(length(t), 4);
            for i = 1:length(t)
                dQuatSequence(i, :) = [cos(angles(i)/2), sin(angles(i)/2)*axis'];
            end
            
            ui32PolyDeg = uint32(5);
            dDomainLB = 0;
            dDomainUB = 1;
            
            % Fit quaternion polynomials
            [dChbvCoeffs, ~, dSwitchIntervals, ~, ~] = ...
                fitAttQuatChbvPolynmials(ui32PolyDeg, t, dQuatSequence', ...
                                       dDomainLB, dDomainUB, false, ui32PolyDeg);
            
            % Test evaluation
            testPoint = 0.5;
            dInterpQuat = evalAttQuatChbvPolyWithCoeffs(ui32PolyDeg, uint32(4), ...
                testPoint, dChbvCoeffs, dSwitchIntervals, dDomainLB, dDomainUB, ui32PolyDeg);
            
            % Verify quaternion properties
            testCase.verifySize(dInterpQuat, [4, 1]);
            quatNorm = norm(dInterpQuat);
            testCase.verifyEqual(quatNorm, 1.0, 'AbsTol', 1e-6, ...
                'Interpolated quaternion should be unit norm');
            
            % Verify interpolation accuracy
            expectedAngle = pi/2 * testPoint;
            expectedQuat = [cos(expectedAngle/2); sin(expectedAngle/2)*axis];
            
            % Compare using quaternion dot product (accounts for double cover)
            dotProduct = abs(dot(dInterpQuat, expectedQuat));
            testCase.verifyGreaterThan(dotProduct, 0.99, ...
                'Interpolated quaternion should be close to expected');
        end
        
        function testEvalRecursiveChbv(testCase)
            %% Test Chebyshev polynomial evaluation function
            ui32PolyDeg = uint32(4);
            ui32PolyMaxDeg = uint32(6);
            
            % Test at specific points where Chebyshev polynomials have known values
            testPoints = [-1, 0, 1];
            
            for i = 1:length(testPoints)
                x = testPoints(i);
                dChbvPoly = EvalRecursiveChbv(ui32PolyDeg, x, ui32PolyMaxDeg);
                
                testCase.verifySize(dChbvPoly, [ui32PolyMaxDeg+1, 1]);
                
                % Check known Chebyshev polynomial values
                % T0(x) = 1, T1(x) = x, T2(x) = 2x^2 - 1, T3(x) = 4x^3 - 3x, T4(x) = 8x^4 - 8x^2 + 1
                expectedValues = [1; x; 2*x^2 - 1; 4*x^3 - 3*x; 8*x^4 - 8*x^2 + 1];
                
                testCase.verifyEqual(dChbvPoly(1:ui32PolyDeg+1), expectedValues, ...
                    'AbsTol', 1e-14, sprintf('Chebyshev evaluation failed at x=%.1f', x));
            end
        end
        
        function testPolynomialDegreeValidation(testCase)
            %% Test validation of polynomial degree constraints
            dInterpDomain = linspace(0, 1, 5)';
            dDataMatrix = ones(1, 5);
            
            % Test with insufficient data points
            ui32PolyDeg = uint32(6); % More coefficients than data points
            
            testCase.verifyError(@() fitChbvPolynomials(ui32PolyDeg, dInterpDomain, ...
                dDataMatrix, 0, 1, false), 'MATLAB:assertion:failed');
        end
        
        function testInterpolationAccuracyHighOrder(testCase)
            %% Test high-order polynomial interpolation accuracy
            % Use a function that should be well-approximated by Chebyshev polynomials
            dInterpDomain = linspace(-1, 1, 100)';
            % Runge function: f(x) = 1/(1 + 25*x^2) - known to be challenging
            dDataMatrix = 1 ./ (1 + 25 * dInterpDomain'.^2);
            
            ui32PolyDeg = uint32(45);
            dDomainLB = -1;
            dDomainUB = 1;
            
            [dChbvCoeffs, ~, ~] = fitChbvPolynomials(ui32PolyDeg, dInterpDomain, dDataMatrix, ...
                                                 dDomainLB, dDomainUB, true);
            
            % Test interpolation at intermediate points
            testPoints = [-0.7, -0.3, 0.2, 0.6];
            for i = 1:length(testPoints)
                x = testPoints(i);
                expectedValue = 1 / (1 + 25 * x^2);
                
                interpolatedValue = evalChbvPolyWithCoeffs(ui32PolyDeg, uint32(1), ...
                    x, dChbvCoeffs, dDomainLB, dDomainUB, ...
                    uint32(length(dChbvCoeffs)), ui32PolyDeg);
                
                relativeError = abs(interpolatedValue - expectedValue) / abs(expectedValue);
                testCase.verifyLessThan(relativeError, 0.1, ...
                    sprintf('High-order interpolation error too large at x=%.1f', x));
            end
        end
        
    end
end
