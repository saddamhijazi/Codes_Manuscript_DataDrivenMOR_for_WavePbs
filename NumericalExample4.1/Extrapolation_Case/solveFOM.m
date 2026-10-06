function [U,M,K,F] = solveFOM(params, nodes, mesh, model, tau, T)
% solveFOM Solves the full order model for the specified parameters and
%          spatial and temporal discretizations.
%
% Inputs:
%   params - Model damage parameters defining the location of the damage,
%   the magnitude of the damage and the width of the damage.
%   nodes  - Coordinates of the spatial discretization nodes.
%   mesh   - Mesh.
%   model  - PDE model.
%   tau    - Time step size.
%   T      - Final simulation time.
%
% Outputs:
%   U      - FOM solution/snapshot matrix.
%   M      - Mass matrix of the spatially discretized model.
%   K      - Stiffness matrix of the spatially discretized model.
%   F      - Force/right-hand-side matrix.

% Parameter Values 
xD = params(1);
yD = params(2);
cD = params(3);
dD = params(4); 

nmb_nodes = size(nodes,2);
nE = findNodes(mesh,"region","Edge",[1 2 3 4]);
nmb_eNodes = size(nE,2);

c = @(location,state) 0.25-(0.25-cD)*exp(-((location.x-xD).^2 + (location.y-yD).^2)/dD);
fcoeff = @(location,state) fcoeffunction(location,state);
specifyCoefficients(model,"m",1,"d",0,"c",c,"a",0,"f",fcoeff);
applyBoundaryCondition(model,"dirichlet","Edge",[1 2 3 4],"u",0);

% Solve with the Theta-Scheme
nmb_timestep = floor(T/tau);
theta = 0.25;

% Extract Matrices
state.time = tau; 
FEMn = assembleFEMatrices(model,"nullspace",state);
K = FEMn.Kc; M = FEMn.M;
U = zeros(nmb_nodes-nmb_eNodes,nmb_timestep);

F = zeros(nmb_nodes-nmb_eNodes,nmb_timestep);
F(:,2) = FEMn.Fc;

%Theta-Scheme
for i = 2:nmb_timestep-1
    state.time = i*tau;
    FEMn = assembleFEMatrices(model,"nullspace",state);
    F(:,i+1) = FEMn.Fc;
    U(:,i+1) = (M + tau^2*theta*K) \ (M*(2*U(:,i)-U(:,i-1)) ...
                + tau^2*K*((2*theta-1)*U(:,i) - theta*U(:,i-1)) ...
                + tau^2*(theta*F(:,i+1) + (1-2*theta)*F(:,i) + theta*F(:,i-1)) );
end
end
