function [dU_dt, d2U_dt2] = finiteDifference8thOrder(U, dt)
    % U: Matrix with each column representing a time instant, each row a degree of freedom
    % dt: Time step size
    % dU_dt: First derivative approximation
    % d2U_dt2: Second derivative approximation
    
    [n_dof, n_time] = size(U);

    % Preallocate derivative matrices
    dU_dt = zeros(n_dof, n_time);
    d2U_dt2 = zeros(n_dof, n_time);

    % Coefficients for the central finite difference scheme (8th-order)
    c = [1/280, -4/105, 1/5, -4/5, 4/5, -1/5, 4/105, -1/280];
    a = [-1/560, 8/315, -1/5, 8/5, -205/72, 8/5, -1/5, 8/315, -1/560];

    % Central finite difference for time steps 5 to n_time - 4
    for k = 5:(n_time - 4)
        dU_dt(:, k) = (c(1) * U(:, k-4) + c(2) * U(:, k-3) + c(3) * U(:, k-2) + ...
                       c(4) * U(:, k-1) + c(5) * U(:, k+1) + c(6) * U(:, k+2) + ...
                       c(7) * U(:, k+3) + c(8) * U(:, k+4)) / dt;

        d2U_dt2(:, k) = (a(1) * U(:, k-4) + a(2) * U(:, k-3) + a(3) * U(:, k-2) + ...
                         a(4) * U(:, k-1) + a(5) * U(:, k)   + a(6) * U(:, k+1) + ...
                         a(7) * U(:, k+2) + a(8) * U(:, k+3) + a(9) * U(:, k+4)) / dt^2;
    end
    
    % Forward finite difference for the first four time steps
    for k = 1:4
        dU_dt(:, k) = (-25/12 * U(:, k) + 4 * U(:, k+1) - 3 * U(:, k+2) + ...
                        4/3 * U(:, k+3) - 1/4 * U(:, k+4)) / dt;
        
        d2U_dt2(:, k) = (35/12 * U(:, k) - 26/3 * U(:, k+1) + 19/2 * U(:, k+2) - ...
                         14/3 * U(:, k+3) + 11/12 * U(:, k+4)) / dt^2;
    end
    
    % Backward finite difference for the last four time steps
    for k = (n_time-3):n_time
        dU_dt(:, k) = (25/12 * U(:, k) - 4 * U(:, k-1) + 3 * U(:, k-2) - ...
                        4/3 * U(:, k-3) + 1/4 * U(:, k-4)) / dt;
        
        d2U_dt2(:, k) = (35/12 * U(:, k) - 26/3 * U(:, k-1) + 19/2 * U(:, k-2) - ...
                         14/3 * U(:, k-3) + 11/12 * U(:, k-4)) / dt^2;
    end
end
