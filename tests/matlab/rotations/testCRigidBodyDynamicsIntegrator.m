classdef testCRigidBodyDynamicsIntegrator < matlab.unittest.TestCase
    methods (Test)
        function testDefaultConstructor(testCase)
            objIntegrator = CRigidBodyDynamicsIntegrator();
            objIntegrator.dInertiaMatrix = eye(3);

            % Default inertia should be identity
            testCase.verifyEqual(objIntegrator.dInertiaMatrix, eye(3));
            % Default quaternion integrator should be a QuaternionIntegrator
            testCase.verifyClass(objIntegrator.objQuatKinIntegrator, 'CQuatKinematicsIntegrator');
        end

        function testCustomConstructor(testCase)
            % Custom inertia and integrator
            dInertia = diag([2,3,4]);
            objQuatInt = CQuatKinematicsIntegrator();
            objIntegrator = CRigidBodyDynamicsIntegrator(dInertia, objQuatInt);
            
            testCase.verifyEqual(objIntegrator.dInertiaMatrix, dInertia);
            testCase.verifySameHandle(objIntegrator.objQuatKinIntegrator, objQuatInt);
        end

        function testEvalRHSAngAccel(testCase)

            % Test static RHS for angular acceleration (single RHS call)
            dInertia        = eye(3);
            dOmegaAngVel    = [0.1;0.2;0.3];
            dExtTorque      = [1;0;0];

            dAngAccel = CRigidBodyDynamicsIntegrator.EvalRHS_AngAccel_(dInertia, dOmegaAngVel, dExtTorque);
            dExpected = dInertia \ (dExtTorque - cross(dOmegaAngVel, dInertia*dOmegaAngVel));

            testCase.verifyEqual(dAngAccel, dExpected, 'AbsTol',1e-12);
        end

        function testIntegrStepRK4Signature(testCase)
            % Ensure RK4 step returns correct sizes
            dQuat0      = [0;1;0;0];
            dOmega0     = zeros(3,1);
            fcnEval     = @(t,w,q) zeros(3,1);

            dTimestamp  = 0.0;
            dDeltaTime  = 0.1;
            [dNextOmega, dNextQuat, dOmegaStagesVals, dStagesTimes] = CRigidBodyDynamicsIntegrator.IntegrStep_RK4_RKMK4(dOmega0, ...
                                                                                                    dQuat0, ...
                                                                                                    fcnEval, ...
                                                                                                    dTimestamp, ...
                                                                                                    dDeltaTime);
            
            testCase.verifyEqual(dOmegaStagesVals, repmat(dNextOmega, 1, 4)); % Mocked dynamics RHS = zero --> omega should remain constant
            testCase.verifySize(dNextOmega, [3 1]);
            testCase.verifySize(dNextQuat,  [4,1]);         % 4 stages × 3 components
            testCase.verifyEqual(dNextQuat, dQuat0, 'AbsTol', 1e-12);
            testCase.verifySize(dOmegaStagesVals, [3, 4]);  % 4 stages × 3 components
            testCase.verifySize(dStagesTimes, [1 4]);
        end

        function testIntegrateStepInterfaceAllConstant(testCase)
            % Test integrateStep_ interface and output sizes
            dInertia    = eye(3);
            dQuat0      = [1;0;0;0];
            dOmega0     = zeros(3,1);
            dTimestamp  = 0.0;
            dDeltaT     = 0.05;
            varTorque   = zeros(3,1);
            % Method signature
            % [dQuatOut, dOmegaOut] = integrateStep_(dTimestamp, ...
            %                                        dInertiaMatrix, ...
            %                                        dQuat0, ...
            %                                        dOmega0, ...
            %                                        dDeltaT, ...
            %                                        varTorque, ...
            %                                        bDecoupledKinoDynamics)

            [dQuatOut, dOmegaOut] = CRigidBodyDynamicsIntegrator.integrateStep_(dTimestamp, ...
                                                                                dInertia, ...
                                                                                dQuat0, ...
                                                                                dOmega0, ...
                                                                                dDeltaT);

            testCase.verifySize(dQuatOut,  [4 1]);
            testCase.verifySize(dOmegaOut, [3 1]);
            testCase.verifyEqual(dQuatOut, dQuat0); 
            testCase.verifyEqual(dOmegaOut, dOmega0);
        end

        function testIntegrateInterfaceAllConstant(testCase)
            % Test high-level integrate method interface
            objIntegrator = CRigidBodyDynamicsIntegrator();

            dTimegrid    = linspace(0,1,6);
            dQuat0       = [1;0;0;0];
            dOmega0      = zeros(3,1);
            varTorque    = zeros(3,1);
            dDeltaT      = 0.1;

            objIntegrator.dInertiaMatrix = eye(3);

            [dQuatSeq, dOmegaSeq] = objIntegrator.integrate(dTimegrid, dQuat0, dOmega0, varTorque, dDeltaT);
            
            testCase.verifySize(dQuatSeq,  [4 numel(dTimegrid)]);
            testCase.verifySize(dOmegaSeq, [3 numel(dTimegrid)]);
        end

        function testIntegrateUsesActualPartialInternalStep(testCase)
            objIntegrator = CRigidBodyDynamicsIntegrator(eye(3), CQuatKinematicsIntegrator());

            dTimegrid = [0.0, 0.25, 0.5];
            dQuat0 = [1;0;0;0];
            dOmega0 = zeros(3,1);
            varTorque = @(dTstamp, dOmega, dQuat0) [1;0;0];

            [~, dOmegaSeq] = objIntegrator.integrate(dTimegrid, ...
                                                     dQuat0, ...
                                                     dOmega0, ...
                                                     varTorque, ...
                                                     0.2, ...
                                                     'rk4_rkmk4', ...
                                                     true, ...
                                                     1.0);

            expectedOmega = [0.0, 0.25, 0.5;
                             0.0, 0.0,  0.0;
                             0.0, 0.0,  0.0];
            testCase.verifyEqual(dOmegaSeq, expectedOmega, 'AbsTol', 1e-10);
        end

        function testIntegrateDeducesInternalStepFromDefaultLimit(testCase)
            objIntegrator = CRigidBodyDynamicsIntegrator(eye(3), CQuatKinematicsIntegrator());

            dTimegrid = [0.0, 0.25, 0.5];
            dQuat0 = [1;0;0;0];
            dOmega0 = zeros(3,1);
            varTorque = @(dTstamp, dOmega, dQuat0) [1;0;0];

            [~, dOmegaSeq] = objIntegrator.integrate(dTimegrid, ...
                                                     dQuat0, ...
                                                     dOmega0, ...
                                                     varTorque, ...
                                                     0.0, ...
                                                     'rk4_rkmk4', ...
                                                     true, ...
                                                     0.1);

            expectedOmega = [0.0, 0.25, 0.5;
                             0.0, 0.0,  0.0;
                             0.0, 0.0,  0.0];
            testCase.verifyEqual(dOmegaSeq, expectedOmega, 'AbsTol', 1e-10);
        end


        function testZeroTorqueNonZeroInitialRate(testCase)

            % With zero torque, angular rate should vary keeping angular momentum constant
            objIntegrator = CRigidBodyDynamicsIntegrator();
            objIntegrator.dInertiaMatrix = eye(3);
            objIntegrator.dInertiaMatrix(1,1) = 5;

            dTimegrid = 0:0.1:1.0;
            dQuat0    = [1;0;0;0];
            dOmega0   = [0.5; -0.3; 0.2];
            varTorque = zeros(3,1);

            [dQuatSeq, dOmegaSeq] = objIntegrator.integrate(dTimegrid, ...
                                                        dQuat0, ...
                                                        dOmega0, ...
                                                        varTorque, ...
                                                        0.01, ...
                                                        'rk4_rkmk4', ...
                                                        true, ...
                                                        1.0);

            testCase.verifyEqual(dOmegaSeq(:,1), dOmega0, 'AbsTol',1e-12);
            for idk = 2:size(dOmegaSeq, 2)
                testCase.verifyNotEqual(dOmegaSeq(:,idk), dOmega0);
            end

            % Compute angular momentum at each time
            dAngMom = objIntegrator.dInertiaMatrix * dOmegaSeq;
            dNormAngMom = vecnorm(dAngMom, 2, 1);
            testCase.verifyLessThanOrEqual( abs(dNormAngMom - dNormAngMom(1))./dNormAngMom, 1e-6)

        end

        function testConstantTorqueNonZeroInitialRate(testCase)

            % TODO

            % With constant torque, angular momentum must not be constant
            I = diag([2,2,2]);
            objIntegrator = CRigidBodyDynamicsIntegrator(I, CQuatKinematicsIntegrator());
            
            dTimegrid = 0:0.2:1.0;
            dQuat0    = [1;0;0;0];
            dOmega0   = [0.1; 0.0; -0.1];
            dTorque   = [0.2; 0.0; 0.0];
            varTorque = @(dTstamp, dOmega, dQuat0) dTorque;

            [~, dOmegaSeq] = objIntegrator.integrate(dTimegrid, ...
                dQuat0, dOmega0, varTorque, 0.2, 'rk4_rkmk4', true, 1.0);

            % Analytical angular acceleration: alpha = I^{-1} * torque = [0.1;0;0]
            alpha = I \ dTorque;

            for k=1:size(dOmegaSeq,2)
                t = dTimegrid(k);
                expectedOmega = dOmega0 + alpha * t;
                testCase.verifyEqual(dOmegaSeq(:,k), expectedOmega, 'AbsTol',1e-6);
            end

        end

        function testTorqueFreeMotionAnalytical(testCase)

            % Validate quaternion and omega against analytical torque-free solution
            objIntegrator = CRigidBodyDynamicsIntegrator();
            
            objIntegrator.dInertiaMatrix = eye(3);
            objIntegrator.dInertiaMatrix(3,3) = 10;

            dTimegrid = 0:0.1:1.0;
            omega0 = [0.0; 0.0; pi]; % rotation about z-axis at pi rad/s
            dQuat0 = [1; 0; 0; 0];
            varTorque = zeros(3,1);
            [dQuatSeq, dOmegaSeq] = objIntegrator.integrate(dTimegrid, dQuat0, omega0, varTorque, 0.1, 'rk4_rkmk4', true, 1.0);
            
            % Analytical: constant omega and quaternion q(t) = [cos(theta/2); 0; 0; sin(theta/2)]
            for k = 1:length(dTimegrid)
            
                t = dTimegrid(k);
                theta = norm(omega0) * t;

                qExp = [cos(theta/2); 0; 0; sin(theta/2)];
                qAct = dQuatSeq(:,k) / norm(dQuatSeq(:,k));
                
                % Account for sign ambiguity
                err1 = norm(qAct - qExp);
                err2 = norm(qAct + qExp);
                testCase.verifyLessThan(min(err1, err2), 1e-6);
                
                % Omega remains constant
                testCase.verifyEqual(dOmegaSeq(:,k), omega0, 'AbsTol',1e-12);
            end
        end

        function testTorqueFreeSymmetricTop(testCase)

            % Torque-free motion for symmetric top (I1=I2 != I3)
            dI1 = 2; 
            dI2 = 2; 
            dI3 = 1;
            
            dI = diag([dI1, dI2, dI3]);
            objIntegrator = CRigidBodyDynamicsIntegrator(dI, CQuatKinematicsIntegrator());

            % Time grid
            dTimegrid = 0:0.1:2.0;

            % Initial angular velocity: small transverse and dominant spin about symmetry axis
            dOmega0 = [0.1; 0.0; 1.0];
            dQuat0 = [1; 0; 0; 0];
            
            varTorque = zeros(3,1);
            [dQuatSeq, dOmegaSeq] = objIntegrator.integrate(dTimegrid, ...
                    dQuat0, dOmega0, varTorque, 0.1, 'rk4_rkmk4', true, 1.0);
            
            % Analytical precession frequency: dOmega_p = (I3 - I1)/I1 * omega3
            dOmega_p = (dI3 - dI1)/dI1 * dOmega0(3);
            dA = dOmega0(1);
            
            for k = 1:length(dTimegrid)

                t = dTimegrid(k);
                dExp1 = dA * cos(dOmega_p * t);
                dExp2 = dA * sin(dOmega_p * t);
                
                % Compare transverse components
                testCase.verifyEqual(dOmegaSeq(1,k), dExp1, 'AbsTol',1e-4);
                testCase.verifyEqual(dOmegaSeq(2,k), dExp2, 'AbsTol',1e-4);
                
                % Spin component remains constant
                testCase.verifyEqual(dOmegaSeq(3,k), dOmega0(3), 'AbsTol',1e-12);
            end
        end

    end
end
