function u = excitiationCustom(t,x,y,f,xc,yc)
% excitiationCustom Evaluates a spatially localized harmonic excitation
%                   with user-defined frequency and center coordinates.
%
% Inputs:
%   t  - Time at which the excitation is evaluated.
%   x  - x-coordinate of the evaluation point.
%   y  - y-coordinate of the evaluation point.
%   f  - Frequency of the harmonic excitation.
%   xc - x-coordinate of the excitation center.
%   yc - y-coordinate of the excitation center.
%
% Outputs:
%   u  - Excitation value. The excitation is nonzero only inside a
%        circular region of radius 0.3 centered at (xc,yc), with time
%        dependence given by 10 + sin(2*pi*f*t) and a spatially varying
%        amplitude.

    if sqrt((x-xc)^2+(y-yc)^2) < 0.3
        u = (10+sin(2*pi*f*t))*exp(1+ 0.09/((x-xc)^2+(y-yc)^2-0.09) );
    else 
        u = 0;
    end 
end
