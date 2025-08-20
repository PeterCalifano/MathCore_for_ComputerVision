function dFiniteDiffJac = ComputeFiniteDiffJacobian(objFcnHandle, dX0diff, dEps, ui32OutputID, ui32NumOfOuts, varargin)
arguments
    objFcnHandle  (1,1) {mustBeA(objFcnHandle, ["function_handle", "string", "char"])}
    dX0diff       (:,1) double {isvector, isnumeric} 
    dEps          (1,1) double {isscalar} = 1e-5
    ui32OutputID  (1,1) uint32 {isscalar} = 1
    ui32NumOfOuts (1,1) uint32 {isscalar} = 1
end
arguments (Repeating)
    varargin % Varargin passed as additional args to function
end

% Compute expansion point map
dF0  = getSelectedHandleOuput(@() objFcnHandle(dX0diff), ui32OutputID, ui32NumOfOuts, varargin);
dFiniteDiffJac = zeros(numel(dF0), length(dX0diff));

d2ndOrdFiniteDiffDenom = 1/(2*dEps);

% Loop over input states to evaluate function 
for ii = 1:length(dX0diff)
   
    dEpsVec = zeros(size(dX0diff));
    dEpsVec(ii) = dEps;
    
    % Compute plus-minus function vector values
    dfPlus  = getSelectedHandleOuput(@() objFcnHandle(dX0diff + dEpsVec), ui32OutputID, ui32NumOfOuts, varargin);
    dfMinus = getSelectedHandleOuput(@() objFcnHandle(dX0diff - dEpsVec), ui32OutputID, ui32NumOfOuts, varargin);

    if isvector(dfPlus)
        dfPlus = dfPlus(:);
    end

    if isvector(dfMinus)
        dfMinus = dfMinus(:);
    end

    dDiffVec = reshape(dfPlus - dfMinus, [], 1);
    dFiniteDiffJac(:, ii) = d2ndOrdFiniteDiffDenom * (dDiffVec);
end

end


%% LOCAL FUNCTION
function outVariable = getSelectedHandleOuput(objFcnHandle, ui32OutputID, ui32NumOfOuts, varargin)
arguments
    objFcnHandle  (1,1) {mustBeA(objFcnHandle, ["function_handle", "string", "char"])}
    ui32OutputID  (1,1) uint32 {isscalar} = 1
    ui32NumOfOuts (1,1) uint32 {isscalar} = 1
end
arguments (Repeating)
    varargin % Varargin passed as additional args to function
end

cellAllOuts = cell(1, max(ui32OutputID, ui32NumOfOuts));

% Call function handle
[cellAllOuts{:}] = objFcnHandle();

% Check index
outVariable = cellAllOuts{ui32OutputID};

end


