function [U, Udot, Uddot] = NewmarkIntegrationLinear(M, K, F, dt, sp, tsim, t_index, z, U0, Udot0, beta, gamma)

% This function implements the implicit Newmark-beta time integration
% method for a *linear* dynamic system with zero-displacement boundary conditions.
%
% Inputs:
%    - M           : Mass matrix (size: sp x sp)
%    - K           : Stiffness matrix (size: sp x sp)
%    - F           : Force matrix (size: sp x N_steps), each column is the force at a time step
%    - dt          : Time step size (scalar)
%    - sp          : Number of degrees of freedom
%    - tsim        : Total simulation time (scalar)
%    - t_index     : Starting time index for reading force input
%    - z           : Indices of degrees of freedom with zero displacement
%    - U0          : Initial displacement (size: sp x 1)
%    - Udot0       : Initial velocity (size: sp x 1)
%    - beta        : Newmark-beta parameter (commonly 0.25 for unconditional stability)
%    - gamma       : Newmark-gamma parameter (commonly 0.5)
%
% Outputs:
%    - U           : Displacement matrix (size: sp x N_steps+1)
%    - Udot        : Velocity matrix (size: sp x N_steps+1)
%    - Uddot       : Acceleration matrix (size: sp x N_steps+1)
%
% The method solves the linear equation:
%       (K + a1) * U_{k+1} = F_hat
% with no nonlinear/Newton iterations.


% Coefficient matrices for Newmark scheme
a1 = M / (beta * dt^2);
a2 = M / (beta * dt);
a3 = (1 - 2 * beta) * M / (2 * beta);

% Initialize matrices for displacement, velocity, and acceleration
num_steps = nearest(tsim / dt);
U = zeros(sp, num_steps);
Udot = zeros(sp, num_steps);
Uddot = zeros(sp, num_steps);

% Initialize based on initial conditions
ukp1 = U0;
udkp1 = Udot0;
F(:, size(F, 2) + 1) = zeros(sp, 1);
uddkp1 = M \ (F(:, t_index) - K * ukp1);


% Store initial conditions
U(:, 1) = ukp1;
Udot(:, 1) = udkp1;
Uddot(:, 1) = uddkp1;

% Precompute constant matrix
K_plus_a1 = K + a1;

% Time-stepping loop
for i = 1:num_steps
    uk = ukp1;
    udk = udkp1;
    uddk = uddkp1;

    % Effective force (linear system)
    F_hat = F(:, i + t_index) + a1 * uk + a2 * udk + a3 * uddk;

    % --------- LINEAR SOLVE (no nonlinear iteration) ---------
    ukp1 = K_plus_a1 \ F_hat;

    % Compute velocity and acceleration using Newmark formulas
    udkp1 = gamma / (beta * dt) * (ukp1 - uk) + ...
            (1 - gamma / beta) * udk + ...
            (1 - gamma / (2 * beta)) * dt * uddk;

    uddkp1 = 1 / (beta * dt^2) * (ukp1 - uk) - ...
             1 / (beta * dt) * udk - ...
            (1 - 2 * beta) / (2 * beta) * uddk;

    % Enforce zero displacement, velocity, and acceleration BCs
    if ~isempty(z)
        ukp1(z) = 0;
        udkp1(z) = 0;
        uddkp1(z) = 0;
    end

    % Store the results
    U(:, i + 1) = ukp1;
    Udot(:, i + 1) = udkp1;
    Uddot(:, i + 1) = uddkp1;
end

end
