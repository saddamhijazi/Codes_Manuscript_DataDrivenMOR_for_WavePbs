function f = fcoeffunction(location,state)
% fcoeffunction Evaluates a spatially localized, time-dependent forcing
%               function defined by a Gaussian spatial distribution and
%               an exponentially damped sinusoidal temporal excitation.
%
% Inputs:
%   location - Spatial coordinates at which the forcing is evaluated.
%              The x- and y-coordinates are accessed through location.x
%              and location.y.
%   state    - PDE solution state containing the current simulation time,
%              accessed through state.time.
%
% Outputs:
%   f        - Value of the forcing function at the specified spatial
%              location and time.

    % Parameters
    A     = 1.0;

    x0    = 2.5;
    y0    = 2.5;

    sigma = 2.0;

    alpha = 0.1;

    f0    = 0.2;
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