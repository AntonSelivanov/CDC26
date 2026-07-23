function res = LMI_th1(q, sigma, N, X, K, L, h)
% LMI_TH1 Evaluates the feasibility of LMIs in Theorem 1 for sampled-data 
% boundary control of the semilinear heat equation.
%
% Paper Reference:
%   A. Selivanov and E. Fridman, "Sample-and-hold control of the semilinear 
%   heat equation via the L2 residue separation," IEEE CDC, 2026.
%
% Requirements:
%   - YALMIP Parser (https://yalmip.github.io/)
%   - Compatible LMI Solver (e.g., SeDuMi, MOSEK, SDPT3, or CSDP)
%
% Inputs:
%   q     - Reaction gain (scalar, q > 0)
%   sigma - Lipschitz constant of nonlinearity (scalar, sigma > 0)
%   N     - Number of dominant modes (integer, N >= 1)
%   X     - Solution matrix from ARE (10a) (N x N)
%   K     - Controller gain vector (1 x N)
%   L     - Observer gain vector (N x 1)
%   h     - Sampling period in seconds (scalar, h > 0)
%
% Output:
%   res   - Logical flag: true (1) if LMIs are feasible, false (0) otherwise

%% 0. Default Output & Parameter Validation
res = false;

if h <= 0
    return; % Sampling period must be strictly positive
end

%% 1. System Matrices & Residual Tail Gain (Eqs. 7, 9, 15)
j_tail = N:1e5;
gamma  = (2/pi) * sum(1 ./ (j_tail.^2 - q - sigma)); % L2 gain of residual modes (9)

mode_indices = (0:N-1)';
A    = q * eye(N) - diag(mode_indices.^2);
B    = [1/sqrt(pi); sqrt(2/pi) * (-1).^(1:N-1)'];
C    = [1/sqrt(pi), sqrt(2/pi) * ones(1, N-1)];
Abar = A - L*C + gamma*sigma*X; % Closed-loop observer matrix from (15)

%% 2. Decision Variables Definition
Pz = sdpvar(N, N, 'symmetric');
Pe = sdpvar(N, N, 'symmetric');
Wz = sdpvar(N, N, 'symmetric');
We = sdpvar(N, N, 'symmetric');

%% 3. Assembly of Main LMI Matrix (Psi)
Psi = blkvar;

% Row 1
Psi(1,1) = Pz*(A - B*K) + (A - B*K)'*Pz + K'*K + (sigma/gamma)*eye(N); 
Psi(1,2) = -Pz*B*K + gamma*sigma*X*Pe + K'*K; 
Psi(1,3) = gamma*sigma*Pz; 
Psi(1,5) = K'*K - Pz*B*K; 
Psi(1,6) = K'*K - Pz*B*K; 
Psi(1,7) = h*(A - B*K)'*Wz; 
Psi(1,8) = h*gamma*sigma*X*We; 

% Row 2
Psi(2,2) = Pe*Abar + Abar'*Pe + K'*K; 
Psi(2,3) = -gamma*sigma*Pe; 
Psi(2,4) = Pe*L; 
Psi(2,5) = K'*K; 
Psi(2,6) = K'*K; 
Psi(2,7) = -h*K'*B'*Wz; 
Psi(2,8) = h*Abar'*We; 

% Row 3
Psi(3,3) = -gamma*sigma*eye(N); 
Psi(3,7) = gamma*sigma*h*Wz; 
Psi(3,8) = -gamma*sigma*h*We; 

% Row 4
Psi(4,4) = -1 / (gamma^2); 
Psi(4,8) = h*L'*We; 

% Row 5
Psi(5,5) = K'*K - ((pi^2)/4)*Wz; 
Psi(5,6) = K'*K; 
Psi(5,7) = -h*K'*B'*Wz; 

% Row 6
Psi(6,6) = K'*K - ((pi^2)/4)*We; 
Psi(6,7) = -h*K'*B'*Wz; 

% Row 7 & 8
Psi(7,7) = -Wz; 
Psi(8,8) = -We;  

Psi = sdpvar(Psi); % Convert block variable to standard SDP matrix

%% 4. Define LMI Constraints
tol = 1e-7; % Small strict positivity margin for numerical stability
LMIs = [Pz >= tol*eye(N), ...
        Pe >= tol*eye(N), ...
        Wz >= tol*eye(N), ...
        We >= tol*eye(N), ...
        Psi <= -tol*eye(size(Psi, 1))]; 

%% 5. Solve LMIs
% Try SeDuMi first, with fallback to any installed SDP solver (e.g., MOSEK, SDPT3)
options = sdpsettings('solver', 'sedumi,mosek,sdpt3,csdp', 'verbose', 0); 
sol = optimize(LMIs, [], options); 

%% 6. Verify Solution Feasibility
% problem == 0 (successfully solved), problem == 4 (numerical issues, needs residual verification)
if sol.problem == 0 || sol.problem == 4
    [primal, dual] = check(LMIs); 
    % Allow a small negative residual (-1e-6) to account for solver precision noise
    res = (min(primal) >= -1e-6) && (min(dual) >= -1e-6); 
end

end