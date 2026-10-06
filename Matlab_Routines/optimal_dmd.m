function Xdmd_opt = optimal_dmd(X, r)
%OPTIMAL_DMD Computes the Optimal Dynamic Mode Decomposition (DMD) solution
%
%   Xdmd_opt = OPTIMAL_DMD(X, r) returns the optimal DMD reconstruction of 
%   the snapshot matrix X using a reduced rank r. The method computes 
%   DMD eigenvalues and modes, then finds the optimal amplitude vector 
%   using a least-squares formulation based on a Vandermonde matrix.
%
%   Inputs:
%       X  - Snapshot matrix of size (n x m), where each column is the 
%            system state at a different time step
%       r  - Target rank for truncated DMD
%
%   Output:
%       Xdmd_opt - Optimal DMD reconstruction of the snapshot matrix X
%                  (size n x m)
%
%   Requirements:
%       - External functions:
%           * dmd: function handle to compute [lambda, Phi] = dmd(X1, X2, r)
%           * vanderm: constructs a Vandermonde matrix for a set of eigenvalues
%
tic
% Split snapshots into time-shifted matrices
X1 = X(:, 1:end-1);
X2 = X(:, 2:end);

% Compute DMD eigenvalues and modes
[lambda, Phi] = dmd(X1, X2, r);

% Construct Vandermonde matrix from eigenvalues
Vand = vanderm(lambda, size(X, 2));

% Construct system for optimal amplitudes
P = (Phi' * Phi) .* conj(Vand * Vand');
q = conj(diag((Vand * X') * Phi));  % Note: X is used as D here

% Check whether P is rank-deficient
if rank(P) < size(P, 1)
    Xdmd_opt = zeros(size(X));
    return
end

% Solve for optimal amplitude vector b using Cholesky factorization
Pl = chol(P, 'lower');
b_opt = Pl' \ (Pl \ q);

tic
% Compute time evolution
Psi = diag(b_opt) * Vand;

% Reconstruct the optimal DMD solution
Xdmd_opt = Phi * Psi;
Xdmd_opt = real(Xdmd_opt);  % Discard imaginary part due to numerical noise
toc
end
