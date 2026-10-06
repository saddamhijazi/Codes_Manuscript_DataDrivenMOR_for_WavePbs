function u = excitation(t,x,y)
% excitation Evaluates a spatially localized harmonic excitation centered
%            at (2.5,2.5) with a radius of 0.3.
%
% Inputs:
%   t - Time at which the excitation is evaluated.
%   x - x-coordinate of the evaluation point.
%   y - y-coordinate of the evaluation point.
%
% Outputs:
%   u - Excitation value. The excitation is nonzero only inside a circular
%       region of radius 0.3 centered at (2.5,2.5), with time dependence
%       given by 10 + sin(2*pi*450000*t) and a spatially varying amplitude.

    if sqrt((x-2.5)^2+(y-2.5)^2) < 0.3
        u = (10+sin(2*pi*450000*t))*exp(1+ 0.09/((x-2.5)^2+(y-2.5)^2-0.09) );
    else 
        u = 0;
    end 
end
