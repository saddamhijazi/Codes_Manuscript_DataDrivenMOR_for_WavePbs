function f = fcoeffunction_custom(location,state,pars)
% fcoeffunction_custom Evaluates a spatially localized, time-dependent
%                       forcing function with user-defined parameters.
%                       The forcing consists of a Gaussian spatial
%                       distribution and an exponentially damped
%                       sinusoidal temporal excitation.
%
% Inputs:
%   location - Spatial coordinates at which the forcing is evaluated.
%              The x- and y-coordinates are accessed through location.x
%              and location.y.
%   state    - PDE solution state containing the current simulation time,
%              accessed through state.time.
%   pars     - Forcing parameters:
%              pars(1) - Amplitude of the excitation.
%              pars(2) - x-coordinate of the Gaussian center.
%              pars(3) - y-coordinate of the Gaussian center.
%              pars(4) - Spatial width parameter of the Gaussian.
%              pars(5) - Exponential decay coefficient.
%              pars(6) - Frequency of the sinusoidal excitation.
%
% Outputs:
%   f        - Value of the forcing function at the specified spatial
%              location and time.

    % Parameters
    A     = pars(1);

    x0    = pars(2);
    y0    = pars(3);

    sigma = pars(4);

    alpha = pars(5);

    f0    = pars(6);
    omega = 2*pi*f0;

    % Spatial Gaussian
    g = exp( -((location.x-x0).^2 + (location.y-y0).^2) ...
              / sigma^2 );

    % Temporal excitation
    h = exp(-alpha*state.time) ...
            .* sin(omega*state.time);

    % Final forcing
    f = A * g .* h;

end