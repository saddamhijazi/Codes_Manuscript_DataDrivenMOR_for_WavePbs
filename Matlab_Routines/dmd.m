function [lambda,Phi] = dmd(X,Y,r)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

[U,S,V]=svd(X,'econ');

if nargin == 2
    r = size(X,2);
end

Ur=U(:,1:r); 
Sr=S(1:r,1:r);
Vr=V(:,1:r);

Atilde=Ur'*Y*Vr/Sr;

[W,D] = eig(Atilde);

Phi = Y*Vr/Sr*W;
Phi = complex(Phi);
lambda = diag(D);

end