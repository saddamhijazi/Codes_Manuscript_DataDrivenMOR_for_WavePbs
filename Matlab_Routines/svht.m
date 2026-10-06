function tau = svht(X, sv)
% svht for sigma unknown

sorted_dims = sort(size(X));
m = sorted_dims(1);
n = sorted_dims(2);
beta = m/n; % ratio between 0 and 1

if nargin == 1
    sv = svd(X);
end
omega_approx = 0.56 * beta^3 - 0.95 * beta^2 + 1.82 * beta + 1.43;

tau = median(sv) * omega_approx;

end