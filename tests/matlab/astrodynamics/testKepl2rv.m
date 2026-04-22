classdef testKepl2rv < matlab.unittest.TestCase
    %% DESCRIPTION
    % Unit tests for kepl2rv (Keplerian elements to Cartesian state conversion).
    % Validates output dimensions, circular orbit properties, round-trip consistency
    % with rv2kepl, energy and angular momentum conservation, periapsis/apoapsis positions,
    % degree input option, and random orbit coverage.
    % -------------------------------------------------------------------------------------------------------------

    properties (Constant)
        dMuEarth = 398600.4418;   % [km^3/s^2]
    end

    methods (Test)

        function testOutputDimension(testCase)
            % Output must always be [6x1]
            dxKepl = [7000; 0.01; deg2rad(45); deg2rad(30); deg2rad(60); deg2rad(0)];
            dxCart = kepl2rv(dxKepl, testCase.dMuEarth);
            testCase.verifySize(dxCart, [6, 1]);
        end

        function testCircularOrbitSpeed(testCase)
            % For a circular orbit (e=0), speed must equal sqrt(mu/r) at every position
            dMu  = testCase.dMuEarth;
            dSma = 7000;
            rng('default');

            for idTrial = 1:10
                dTA = 2*pi * rand();
                dxKepl = [dSma; 0.0; deg2rad(45); deg2rad(30); deg2rad(60); dTA];
                dxCart = kepl2rv(dxKepl, dMu);

                dSpeed     = norm(dxCart(4:6));
                dExpSpeed  = sqrt(dMu / dSma);
                testCase.verifyEqual(dSpeed, dExpSpeed, 'RelTol', 1e-10, ...
                    sprintf('Circular orbit speed incorrect at trial %d (TA=%.3f rad)', idTrial, dTA));
            end
        end

        function testCircularOrbitRadiusConstant(testCase)
            % Circular orbit: radius must equal SMA at every true anomaly
            dMu  = testCase.dMuEarth;
            dSma = 8500;
            rng(3);

            for idTrial = 1:10
                dTA = 2*pi * rand();
                dxKepl = [dSma; 0.0; deg2rad(28.5); deg2rad(0); deg2rad(0); dTA];
                dxCart = kepl2rv(dxKepl, dMu);

                dRadius = norm(dxCart(1:3));
                testCase.verifyEqual(dRadius, dSma, 'RelTol', 1e-10, ...
                    sprintf('Circular orbit radius must equal SMA (trial %d)', idTrial));
            end
        end

        function testPeriapsisAndApoapsis(testCase)
            % At TA=0 (periapsis): r = a*(1-e). At TA=pi (apoapsis): r = a*(1+e)
            dMu  = testCase.dMuEarth;
            dSma = 10000;
            dEcc = 0.3;

            dxPeri  = kepl2rv([dSma; dEcc; deg2rad(45); 0; 0; 0],   dMu);
            dxApo   = kepl2rv([dSma; dEcc; deg2rad(45); 0; 0; pi],  dMu);

            testCase.verifyEqual(norm(dxPeri(1:3)), dSma*(1-dEcc), 'RelTol', 1e-10, ...
                'Periapsis radius must be a*(1-e)');
            testCase.verifyEqual(norm(dxApo(1:3)),  dSma*(1+dEcc), 'RelTol', 1e-10, ...
                'Apoapsis radius must be a*(1+e)');
        end

        function testSpecificOrbitalEnergy(testCase)
            % Specific energy must equal -mu/(2a) for any position on the orbit
            dMu  = testCase.dMuEarth;
            rng(7);

            for idTrial = 1:15
                dSma = 7000 + 20000 * rand();
                dEcc = 0.5 * rand();
                dxKepl = [dSma; dEcc; pi*rand(); 2*pi*rand(); 2*pi*rand(); 2*pi*rand()];

                dxCart = kepl2rv(dxKepl, dMu);
                dRadius = norm(dxCart(1:3));
                dSpeed  = norm(dxCart(4:6));

                dEnergyActual   = 0.5*dSpeed^2 - dMu/dRadius;
                dEnergyExpected = -dMu / (2*dSma);

                testCase.verifyEqual(dEnergyActual, dEnergyExpected, 'RelTol', 1e-9, ...
                    sprintf('Specific energy mismatch at trial %d', idTrial));
            end
        end

        function testAngularMomentumMagnitude(testCase)
            % Angular momentum |h| = sqrt(mu * a * (1-e^2)) for any point on the orbit
            dMu  = testCase.dMuEarth;
            rng(11);

            for idTrial = 1:15
                dSma = 7000 + 15000 * rand();
                dEcc = 0.4 * rand();
                dxKepl = [dSma; dEcc; pi*rand(); 2*pi*rand(); 2*pi*rand(); 2*pi*rand()];

                dxCart = kepl2rv(dxKepl, dMu);
                dHvec  = cross(dxCart(1:3), dxCart(4:6));
                dHnorm = norm(dHvec);

                dHexpected = sqrt(dMu * dSma * (1 - dEcc^2));
                testCase.verifyEqual(dHnorm, dHexpected, 'RelTol', 1e-9, ...
                    sprintf('Angular momentum magnitude mismatch at trial %d', idTrial));
            end
        end

        function testRoundTripWithRv2kepl(testCase)
            % kepl2rv --> rv2kepl must recover original elements (within numerical tolerance)
            dMu = testCase.dMuEarth;
            rng('default');

            for idTrial = 1:20
                dSma  = 7000 + 25000 * rand();
                dEcc  = 0.05 + 0.4 * rand();   % Avoid near-circular/near-equatorial singularities
                dIncl = 0.1 + (pi - 0.2) * rand();
                dRaan = 2*pi * rand();
                dArgP = 2*pi * rand();
                dTA   = 2*pi * rand();
                dxKepl0 = [dSma; dEcc; dIncl; dRaan; dArgP; dTA];

                dxCart  = kepl2rv(dxKepl0, dMu);
                dxKepl1 = rv2kepl(dxCart, dMu);

                % SMA and eccentricity
                testCase.verifyEqual(dxKepl1(1), dxKepl0(1), 'RelTol', 1e-8, ...
                    sprintf('SMA round-trip failure at trial %d', idTrial));
                testCase.verifyEqual(dxKepl1(2), dxKepl0(2), 'AbsTol', 1e-10, ...
                    sprintf('Eccentricity round-trip failure at trial %d', idTrial));

                % Angles: use wrapped difference
                for idElem = 3:6
                    dDiff = abs(atan2(sin(dxKepl1(idElem) - dxKepl0(idElem)), ...
                                     cos(dxKepl1(idElem) - dxKepl0(idElem))));
                    testCase.verifyLessThan(dDiff, 1e-8, ...
                        sprintf('Angle element %d round-trip failure at trial %d', idElem, idTrial));
                end
            end
        end

        function testDegreeInputMatchesRadInput(testCase)
            % charUnit="deg" must produce the same output as charUnit="rad" with converted angles
            dMu  = testCase.dMuEarth;
            dSma = 9000;
            dEcc = 0.15;
            dInclDeg = 51.6;
            dRaanDeg = 120.0;
            dArgPDeg = 45.0;
            dTADeg   = 200.0;

            dxKepl_deg = [dSma; dEcc; dInclDeg; dRaanDeg; dArgPDeg; dTADeg];
            dxKepl_rad = [dSma; dEcc; deg2rad(dInclDeg); deg2rad(dRaanDeg); deg2rad(dArgPDeg); deg2rad(dTADeg)];

            dxCart_deg = kepl2rv(dxKepl_deg, dMu, "deg");
            dxCart_rad = kepl2rv(dxKepl_rad, dMu, "rad");

            testCase.verifyEqual(dxCart_deg, dxCart_rad, 'AbsTol', 1e-10, ...
                'deg and rad inputs must produce identical Cartesian states');
        end

        function testVelocityPerpendicularToPositionForCircularOrbit(testCase)
            % Circular orbit: velocity must be perpendicular to position at every point (r·v = 0)
            dMu  = testCase.dMuEarth;
            dSma = 6800;
            rng(42);

            for idTrial = 1:10
                dTA    = 2*pi * rand();
                dxKepl = [dSma; 0.0; deg2rad(35); deg2rad(80); deg2rad(20); dTA];
                dxCart = kepl2rv(dxKepl, dMu);

                dRdotV = dot(dxCart(1:3), dxCart(4:6));
                testCase.verifyEqual(dRdotV, 0.0, 'AbsTol', 1e-8, ...
                    sprintf('r·v must be zero for circular orbit (trial %d)', idTrial));
            end
        end

        function testEccentricityVectorDirection(testCase)
            % The eccentricity vector must point toward periapsis (TA=0 position unit vector)
            dMu  = testCase.dMuEarth;
            dSma = 12000;
            dEcc = 0.25;
            dIncl = deg2rad(30);
            dRaan = deg2rad(45);
            dArgP = deg2rad(90);

            % Position at TA=0 (periapsis)
            dxCartPeri = kepl2rv([dSma; dEcc; dIncl; dRaan; dArgP; 0], dMu);
            dUnitPeri  = dxCartPeri(1:3) / norm(dxCartPeri(1:3));

            % Eccentricity vector from a general position
            dxCartGen = kepl2rv([dSma; dEcc; dIncl; dRaan; dArgP; deg2rad(60)], dMu);
            dRgen = dxCartGen(1:3);
            dVgen = dxCartGen(4:6);
            dEvec = (norm(dVgen)^2/dMu - 1/norm(dRgen)) * dRgen - (dot(dRgen, dVgen)/dMu) * dVgen;
            dUnitE = dEvec / norm(dEvec);

            testCase.verifyEqual(dUnitE, dUnitPeri, 'AbsTol', 1e-9, ...
                'Eccentricity vector must point to periapsis');
        end

    end
end
