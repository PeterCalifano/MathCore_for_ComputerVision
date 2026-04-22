function [fig, o_dDCM_Target2Fixed] = PlotAttitudeQuat(i_dQuatSeq, ...
    i_dOriginPos, ...
    i_bIS_JPL_CONV, ...
    i_bPlotFrame, ...
    i_bAnimateFlag, ...
    i_dPauseTime)
%% PROTOTYPE
% [fig, o_dDCM_Target2Fixed] = PlotAttitudeQuat(i_dQuatSeq, ...
%     i_dOriginPos, ...
%     i_bIS_JPL_CONV, ...
%     i_bPlotFrame, ...
%     i_bAnimateFlag, ...
%     i_dPauseTime)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function plotting the target reference frame defined by input quaternion
% sequence with respect to a FIXED generic frame. 
% Only the tips of the frame axis are plotted if i_bPlotFrame is FALSE.
% NOTE: Input quaternion must be from TARGET to IN (q_TARGETwrtIN)
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% i_dQuatSeq:       [4, Nq]    Attitude quaternion sequence
% i_dOriginPos:     [3, Nq]    Frame origin position vector or sequence
% i_bIS_JPL_CONV:   [1]        Boolean flag for quaternion convention
% i_bPlotFrame:     [1]        Boolean flag to enable reference frame plot
% i_bAnimateFlag:   [1]        Boolean flag to enable animation of sequence
% i_dPauseTime:     [1]        Animation pause step in [s]
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% fig:                  [1]         Figure handle created by the function
% o_dDCM_Target2Fixed:  [3, 3, Nq]  Rotation matrix sequence from quat.
%                                   converting from Target to Fixed frame
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 17-10-2023    Pietro Califano     Adapted from Attitude Path Planner research codes (by Massi). Verified.
% 22-10-2023    Pietro Califano     Added option to disable animation.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% Quat2DCM()
% DefaultPlotOpts()
% -------------------------------------------------------------------------------------------------------------
%% Function code
Nq = size(i_dQuatSeq, 2);

if nargout > 2
    o_dDCM_Target2Fixed = zeros(3, 3, Nq);
end

% Default values and options
if not(exist('i_dOriginPos', 'var'))
    i_dOriginPos = zeros(3, 1);
end

if not(exist('i_bIS_JPL_CONV', 'var'))
    i_bIS_JPL_CONV = true; % JPL convention in Quat2DCM
end

if not(exist('i_bPlotFrame', 'var'))
    i_bPlotFrame = false; % Do not plot frame if false
end

if not(exist('i_dPauseTime', 'var'))
    i_dPauseTime = 0.0001; % Use default pause time
end

if not(exist('i_bAnimateFlag', 'var'))
    i_bAnimateFlag = false; % Do not animate by default
end

%% Determine execution path
posNorm = norm(i_dOriginPos(:, 1));

if i_bPlotFrame == false
    disp('Wireframe DISABLED: Animation is automatically disabled.')
    i_bAnimateFlag = false;
end

if i_bAnimateFlag == false && posNorm > 0
    if i_bPlotFrame == false
        disp('Wireframe DISABLED, but Origin DISPLACED: Frame plotting enforced.')
    else
        disp('Animation DISABLED and Origin DISPLACED: Frame is automatically plotted.')
    end
    i_bPlotFrame = true;
    i_dLineWidth = 0.25;
elseif posNorm > 0
    i_bPlotFrame = true;
    i_dLineWidth = 1.1;
end

mark0 = 'o';
markSeq = '.';
markEnd = 'x';

arrowLength = 1.2;
lineWidth = 1.1;
markerSize = 8;
scaleView = 1.5;
posScaleView = 1;
Nlim = 250;

% PLOTTING ROUTINE
if i_bAnimateFlag == true
    fig = figure('WindowState', 'maximized', 'Name', 'Attitude motion animation');
else
    fig = figure('Name', 'Attitude motion static plot');
end
hold on;

% Draw FIXED frame "Inertial"
quiver3(0, 0, 0, arrowLength + posNorm/2, 0, 0,'k', ...
    'Linewidth', 1.1); text(arrowLength + posNorm/2, 0, 0, 'X');

quiver3(0, 0, 0, 0, arrowLength + posNorm/2, 0,'k', ...
    'Linewidth', 1.1); text(0, arrowLength + posNorm/2, 0, 'Y');

quiver3(0, 0, 0, 0, 0, arrowLength + posNorm/2,'k', ...
    'Linewidth', 1.1); text(0, 0, arrowLength + posNorm/2, 'Z');

% Set axis limits and scalings
axis equal
xlim([-(scaleView*arrowLength + posNorm * posScaleView), (scaleView*arrowLength + posNorm * posScaleView)])
ylim([-(scaleView*arrowLength + posNorm * posScaleView), (scaleView*arrowLength + posNorm * posScaleView)])
zlim([-(scaleView*arrowLength + posNorm * posScaleView), (scaleView*arrowLength + posNorm * posScaleView)])

% Convert first quaternion to DCM (initial orientation)
dAq0 = Quat2DCM(i_dQuatSeq(:, 1), i_bIS_JPL_CONV); 

if nargout > 1
    o_dDCM_Target2Fixed(:, :, 1) = dAq0;
end

% Initial Triad axis in FIXED frame
u0 = dAq0 * [1;0;0];
v0 = dAq0 * [0;1;0];
w0 = dAq0 * [0;0;1];

% Plot initial triad
if posNorm == 0
    vec0(1) = plot3(u0(1), u0(2), u0(3), ['r', mark0],'Linewidth', lineWidth, 'MarkerSize', markerSize);
    vec0(2) = plot3(v0(1), v0(2), v0(3), ['g', mark0],'Linewidth', lineWidth, 'MarkerSize', markerSize);
    vec0(3) = plot3(w0(1), w0(2), w0(3), ['b', mark0],'Linewidth', lineWidth, 'MarkerSize', markerSize);
end

if i_bPlotFrame
    plotRefFrame(u0, v0, w0, arrowLength, i_dOriginPos(:, 1), i_dPauseTime, i_bAnimateFlag, i_dLineWidth);
end

if exist('vec0', 'var')
    delete(vec0);
end

if Nq > Nlim
    seqStep = floor(Nq/500);
    idSeq = 2:seqStep:Nq;

    if idSeq(end) ~= Nq
        idSeq(end) = Nq;
    end
else
    idSeq = 2:Nq;
end

DefaultPlotOpts();
hold on;

% Plot trajectory
if size(i_dOriginPos, 2) > 1
    plot3(i_dOriginPos(1, :), i_dOriginPos(2, :), i_dOriginPos(3, :), 'k--', 'LineWidth', 1.1)
end

% Plot attitude evolution
for ii = idSeq % Loop through quaternion sequence

    dAq = Quat2DCM(i_dQuatSeq(:, ii), i_bIS_JPL_CONV); % exctraction of the DCM at each time step

    if size(i_dOriginPos, 2) > 1
        displOrigin = i_dOriginPos(:, ii);
    else
        displOrigin = i_dOriginPos(:, 1);
    end

    u = dAq * [1;0;0]; % first rotating axis in fixed frame
    v = dAq * [0;1;0]; % second rotating axis in fixed frame
    w = dAq * [0;0;1]; % third rotating axis in fixed frame

    if ii == Nq
        markPlot = markEnd;
    else
        markPlot = markSeq;
    end

    if nargout > 1
        o_dDCM_Target2Fixed(:, :, ii) = dAq;
    end

    % Set axis limits and scalings
    %         view([90, 0])
    axis equal
    xlim([-(scaleView*arrowLength + posNorm * posScaleView), (scaleView*arrowLength + posNorm * posScaleView)]);
    ylim([-(scaleView*arrowLength + posNorm * posScaleView), (scaleView*arrowLength + posNorm * posScaleView)]);
    zlim([-(scaleView*arrowLength + posNorm * posScaleView), (scaleView*arrowLength + posNorm * posScaleView)]);

    if i_bPlotFrame 

        plotRefFrame(u, v, w, arrowLength, displOrigin, i_dPauseTime, i_bAnimateFlag, i_dLineWidth);

    elseif i_bPlotFrame == false && i_bAnimateFlag == false

        posNormTmp2 = norm(displOrigin);

        u = (arrowLength + posNormTmp2/2) * u + displOrigin;
        v = (arrowLength + posNormTmp2/2) * v + displOrigin;
        w = (arrowLength + posNormTmp2/2) * w + displOrigin;

        arrowTip(1) = plot3(u(1), u(2), u(3), ['r', markPlot], 'Linewidth', lineWidth, 'MarkerSize', markerSize);
        arrowTip(2) = plot3(v(1), v(2), v(3), ['g', markPlot], 'Linewidth', lineWidth, 'MarkerSize', markerSize);
        arrowTip(3) = plot3(w(1), w(2), w(3), ['b', markPlot], 'Linewidth', lineWidth, 'MarkerSize', markerSize);

        %         if i_bAnimateFlag == true
        %             % Pause and delete previous quiver3 to create animation effect
        %             pause(i_dPauseTime);
        %             delete(arrowTip);
        %         end

    end

    %     campos('auto')
    %     if i_bPlotFrame
    %         quiverObj = findobj('Type', 'quiver');
    %         delete(quiverObj(end:end-2));
    %     end

end

%% LOCAL FUNCTION
    function [vec] = plotRefFrame(i_duTarget, i_dvTarget, i_dwTarget, i_dArrowLength, i_dOriginPos, i_dPauseTime, i_bAnimateFlag, i_dLineWidth)

        if not(exist('i_dPauseTime', 'var'))
            i_dPauseTime = 0;
        end
        if not(exist('i_dAnimateFlag', 'var'))
            i_bAnimateFlag = false;
        end

        if not(exist('i_dLineWidth', 'var'))
            i_dLineWidth = 1.1;
        end

        posNormTmp = norm(i_dOriginPos);

        % Define vectors to plot
        Xvec = (i_dArrowLength + posNormTmp/2) * i_duTarget;
        Yvec = (i_dArrowLength + posNormTmp/2) * i_dvTarget;
        Zvec = (i_dArrowLength + posNormTmp/2) * i_dwTarget;

        % Plot vectors
        vec(1) = quiver3(i_dOriginPos(1), i_dOriginPos(2), i_dOriginPos(3), ...
            Xvec(1), Xvec(2), Xvec(3), 'r', ...
            'Linewidth', i_dLineWidth*0.75, 'AutoScale', 'on');

        vec(2) = quiver3(i_dOriginPos(1), i_dOriginPos(2), i_dOriginPos(3), ...
            Yvec(1), Yvec(2), Yvec(3), 'g', ...
            'Linewidth', i_dLineWidth*0.75, 'AutoScale', 'on');

        vec(3) = quiver3(i_dOriginPos(1), i_dOriginPos(2), i_dOriginPos(3), ...
            Zvec(1), Zvec(2), Zvec(3), 'b', ...
            'Linewidth', i_dLineWidth, 'AutoScale', 'on');

        if i_bAnimateFlag == true
            pause(i_dPauseTime);
            delete(vec);
        end
    end
end
