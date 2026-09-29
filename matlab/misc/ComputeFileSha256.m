function charDigest = ComputeFileSha256(charFilePath)
%% SIGNATURE
% charDigest = ComputeFileSha256(charFilePath)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Compute SHA-256 over raw file bytes through the bundled DataHash package.
% Preserve the file-integrity interface used by simulation asset loaders;
% delegate bounded file reading and digest calculation to DataHash.
% Use this host-side utility with the MATLAB JVM enabled.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% charFilePath    Path to an existing local file.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% charDigest      64-character lowercase hexadecimal SHA-256 digest.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 27-09-2026  Pietro Califano, Codex gpt-6  Add file integrity check for selected kernel bundles.
% 29-09-2026  Pietro Califano, Codex        Reuse DataHash file mode from MathCore.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% DataHash (bundled under matlab/misc/DataHash), MATLAB JVM.
% -------------------------------------------------------------------------------------------------------------

arguments (Input)
    charFilePath (1, :) char {mustBeFile}
end

arguments (Output)
    charDigest (1, 64) char
end

% Hash the file contents without adding MATLAB array metadata to the digest.
charDigest = DataHash(charFilePath, 'SHA-256', 'file', 'hex');

end
