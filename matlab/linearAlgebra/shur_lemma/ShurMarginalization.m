function [dMarginalizedCov11, dMarginalizedCov22] = ShurMarginalization(dCovMatrix, ...
                                              ui16FirstMargStateIdx, ...
                                              ui16CovarianceSize) %#codegen
arguments
    dCovMatrix              (:,:) double
    ui16FirstMargStateIdx   (1,1) uint16
    ui16CovarianceSize      (1,1) uint16 = size(dCovMatrix, 1)
end
%% PROTOTYPE
% [dMarginalCov] = ShurMarginalization(dCovMatrix, ui16FirstMargStateIdx) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Computation of conditional covariance matrix from a Joint PDF covariance using Shur Complement. Correlations 
% are removed from the prior covariance of the marginalized states. The algorithm can work only on the
% bottom portion of the state vector.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dCovMatrix              (:,:) double
% ui16FirstMargStateIdx   (1,1) uint15
% ui16CovarianceSize      (1,1) uint16 = size(dCovMatrix, 1)
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dMarginalCov
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 21-04-2024        Pietro Califano         First simple version coded.
% 03-03-2025        Pietro Califano         Update of function to compute both conditional subblocks.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code
if coder.target('MATLAB') || coder.target('MEX')
    assert(size(dCovMatrix, 1) == size(dCovMatrix, 2));
    assert(ui16FirstMargStateIdx <= size(dCovMatrix, 1));
    assert(isscalar(ui16FirstMargStateIdx));
end

% Extract prior covariance blocks
dLAMBDA22 = dCovMatrix(ui16FirstMargStateIdx:end, ui16FirstMargStateIdx:end);
dLAMBDA12 = dCovMatrix(1:ui16FirstMargStateIdx-1, ui16FirstMargStateIdx:end);
dLAMBDA11 = dCovMatrix(1:ui16FirstMargStateIdx-1, 1:ui16FirstMargStateIdx-1);

% Compute stable inverse using Cholesky
dCholLAMBDA22 = chol(dLAMBDA22, 'lower'); 
dX = dCholLAMBDA22 \ dLAMBDA12';

% Compute Schur complement
dMarginalizedCov11 = zeros( ui16CovarianceSize, ui16CovarianceSize );
dMarginalizedCov11(1:ui16FirstMargStateIdx-1, 1:ui16FirstMargStateIdx-1) = dLAMBDA11 - (dX' * dX);

% Compute LAMBDA22 marginal
if nargout > 1
    dMarginalizedCov22 = zeros( ui16CovarianceSize, ui16CovarianceSize );
    dMarginalizedCov22( ui16FirstMargStateIdx:end, ui16FirstMargStateIdx:end ) = ...
                    dLAMBDA22 - transpose(dLAMBDA12) * ( dLAMBDA11 \ dLAMBDA12 );
end

end
