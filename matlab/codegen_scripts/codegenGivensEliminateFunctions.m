clear
close all
clc

%% CODEGEN SCRIPT for Givens elimination functions
% 26-04-2026    Pietro Califano     Codegen the row, column, and QR Givens elimination entrypoints.

charScriptDir = fileparts(mfilename('fullpath'));
charMatlabRoot = fileparts(charScriptDir);
addpath(genpath(charMatlabRoot), '-begin');

dMatrix = coder.typeof(0, [8, 8], [1, 1]);
ui32TargetSubscript = coder.typeof(uint32(0), [1, 2], [0, 0]);
ui32AuxIndex = coder.typeof(uint32(0), [1, 1], [0, 0]);
ui16ValidPtr = coder.typeof(uint16(0), [1, 1], [0, 0]);

cfg = coder.config('mex');
cfg.GenerateReport = false;
cfg.RowMajor = false;

codegen('-config', cfg, ...
        'GivensEliminateRow', ...
        '-args', {dMatrix, ui32TargetSubscript, ui32AuxIndex}, ...
        '-o', 'GivensEliminateRow_MEX');

codegen('-config', cfg, ...
        'GivensEliminateColumn', ...
        '-args', {dMatrix, ui32TargetSubscript, ui32AuxIndex}, ...
        '-o', 'GivensEliminateColumn_MEX');

codegen('-config', cfg, ...
        'GivensEliminateQR', ...
        '-args', {dMatrix, ui16ValidPtr, ui16ValidPtr}, ...
        '-o', 'GivensEliminateQR_MEX');
