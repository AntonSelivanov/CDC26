function [feas, K, L, X] = ARE_th1(q, sigma, N)
% ARE_TH1 Solves the Algebraic Riccati Equations (10)-(11) for controller
% synthesis under L2 residue separation.
%
% Paper Reference:
%   A. Selivanov and E. Fridman, "Sample-and-hold control of the semilinear 
%   heat equation via the L2 residue separation," IEEE CDC, 2026.
%
% Inputs:
%   q     - Reaction gain of the PDE (scalar, q > 0)
%   sigma - Lipschitz constant of nonlinearity (scalar, sigma > 0)
%   N     - Number of dominant modes for control design (integer, N >= 1)
%
% Outputs:
%   feas  - Logical flag: true (1) if AREs are feasible, false (0) otherwise
%   K     - Controller gain vector (1 x N), empty if infeasible
%   L     - Observer gain vector (N x 1), empty if infeasible
%   X     - Positive-definite solution matrix of control ARE (N x N)

%% 0. Default Outputs & Parameter Checks
feas = false; 
K    = []; 
L    = []; 
X    = [];

% Ensure residual modes are open-loop stable: lambda_N = N^2 > q + sigma
if (N^2 - q - sigma) <= 0
    warning('ARE_th1:ResidualInstability', ...
        'N = %d is too small for q = %.2f, sigma = %.2f. Require N^2 > q + sigma.', N, q, sigma);
    return;
end

positivityGap = 1e-10;

%% 1. System Matrices of Dominant Modes (Eq. 7)
mode_indices = (0:N-1)';
A = q * eye(N) - diag(mode_indices.^2);

% Mode projections for boundary actuation and boundary measurement
B = [1/sqrt(pi); sqrt(2/pi) * (-1).^(1:N-1)'];
C = [1/sqrt(pi), sqrt(2/pi) * ones(1, N-1)];

%% 2. L2 Gain of Residual Modes (Eq. 9)
% Compute tail mode infinite sum approximation (j >= N)
j_tail = N:1e5;
gamma  = (2/pi) * sum(1 ./ (j_tail.^2 - q - sigma));

if gamma <= 0
    return;
end

%% 3. Solve Algebraic Riccati Equations (Eq. 10 and 11)
Q = (sigma / gamma) * eye(N);
G = (sigma * gamma) * eye(N);

try
    X = icare(A, B, Q, [], [], [], G);    % Control ARE (10a)
    Z = icare(A', C', Q, [], [], [], G);  % Observer ARE (10b)
catch
    return; % Return infeasible if ARE solver fails to converge
end

% Validate solver outputs
if isempty(X) || isempty(Z) || any(isnan(X(:))) || any(isnan(Z(:)))
    X = [];
    return;
end

% Enforce numerical symmetry
X = (X + X') / 2;
Z = (Z + Z') / 2;

%% 4. Verify Positivity & Spectral Radius Condition (Eq. 11)
eigX = real(eig(X));
eigZ = real(eig(Z));

if min(eigX) > -positivityGap && min(eigZ) > -positivityGap
    % Check spectral condition: lambda_max(X * Z) < 1 / gamma^2
    eigXZ = real(eig(X * Z));

    if max(eigXZ) < (1 / gamma^2)
        feas = true;

        % Compute Feedback and Observer Gains (Eq. 12)
        K = B' * X;
        L = (Z / (eye(N) - (gamma^2) * (X * Z))) * C';
    end
end

end