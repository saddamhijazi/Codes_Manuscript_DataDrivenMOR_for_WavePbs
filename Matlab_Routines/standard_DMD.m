function Y = standard_DMD(Snapshots,r)

X = Snapshots; 
Y1 = X(:,1:end-1); 
Y2 =  X(:,2:end); 

[U,S,V] = svd(Y1,'econ');

Ur = U(:,1:r); 
Sr = S(1:r,1:r);
Vr = V(:,1:r);

Atilde = Ur'*Y2*Vr/Sr;

[m,n] = size(X);
Yr = zeros(r,n);

tic
Yr(:,1) = Ur'*X(:,1);

for i = 2:n
    Yr(:,i) = Atilde*Yr(:,i-1);
end

Y = zeros(m,n);
Y = Ur*Yr;

end

       
