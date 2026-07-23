% =========================================================================
% Main Example Script
% =========================================================================
% Paper: A. Selivanov and E. Fridman, "Sample-and-hold control of the 
%        semilinear heat equation via the L2 residue separation," 
%        IEEE Conference on Decision and Control (CDC), 2026.
%
% Requirements: 
%   - YALMIP parser (https://yalmip.github.io/)
%   - Compatible LMI Solver (e.g., SeDuMi, MOSEK, SDPT3, or CSDP)
% =========================================================================

clc; clear; close all;

%% 1. System Parameters
q     = 0.4;   % Reaction gain (> 0 leads to open-loop instability)
sigma = 0.1;   % Maximum admissible Lipschitz constant
N     = 8;     % Number of dominant modes for control design

fprintf('=== Continuous-time Controller Design ===\n');
fprintf('Parameters: q = %.2f, sigma = %.2f, N = %d\n\n', q, sigma, N);

%% 2. Continuous-time Control Design (ARE Feasibility)
[feas, K, L, X] = ARE_th1(q, sigma, N); 

if ~feas
    error('AREs are not feasible for the given parameters.');
else 
    disp('ARE feasibility confirmed.');
    disp('Control gain K =');
    disp(K);
    disp('Observer gain L =');
    disp(L);
end

%% 3. Check LMI Feasibility for a Given Sampling Period
fprintf('=== Sampled-Data Feasibility Check ===\n');

q = .4;    % Slightly smaller to distance from the critical case 
h = .1;    % Sampling period to test [s]

feas = LMI_th1(q, sigma, N, X, K, L, h); 

if feas
    fprintf('SUCCESS: The LMIs are FEASIBLE for h = %.4f s\n\n', h);
else 
    fprintf('FAILURE: The LMIs are NOT FEASIBLE for h = %.4f s\n\n', h);
end