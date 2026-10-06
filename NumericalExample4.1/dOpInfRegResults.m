close all; clear; clc;
% Damage Parameters
xD = 3;  % Position 
yD = 3; 
cD = 0.05; % Wave speed in the damaged area 
dD = 0.06; % Width damage variable
params = [xD,yD,cD,dD];

T = 5; 
tau = 0.02; 
dx = 0.05;

nmb_timestep = T/tau; % Number of time steps
theta = 0.25;
[nodes,mesh,model] = gitter(dx);
nmb_nodes = size(nodes,2);
nE = findNodes(mesh,"region","Edge",[1 2 3 4]);
nmb_eNodes = size(nE,2);

tic
% Solve the FOM
[U,M,K] = solveFOM(params,nodes,mesh,model,tau,T);
toc

tsim = T;
nSkip = 0.1/tau;
dt = tau;
t = nSkip*dt:dt:tsim-dt;
load('ScalingVals.mat')

U = U / s_U;

U = U(:,nSkip+1:end);

N = size(U,1);
Nt = size(U,2);
%% Solve q different full order model setups in order to regulairze the
% operator inference problem, the x-center of the wave excitiation is
% varied in order to result in different input/excititaion signals
tstart = tic;
q = 6;
U_multiple = zeros(N,q*Nt);
xc = linspace(2.5,3.5,q);
for j=1:q
    [U_s,~,~] = solveFOM_Custom(params,nodes,mesh,model, ...
        tau,T,450000,xc(j),2.5);
    U_s = U_s(:,nSkip+1:end);
    U_multiple(:,(j-1)*Nt+1:j*Nt) = U_s;
end
U_multiple = U_multiple / norm(U_multiple,'fro');
tEnd_multiTraj = toc(tstart);
r = 40;

[V,~,~] = svd(U,'econ');
Phi = V(:,1:r);

% Compute the continuous intrusive POD Operators
M_r = Phi' * M * Phi;
K_r = Phi' * K * Phi;

tstart = tic;
U_hat = Phi' * U_multiple;

Uk_hat = zeros(r,(Nt-2)*q);
Ukp1_hat = zeros(r,(Nt-2)*q);
Ukp2_hat = zeros(r,(Nt-2)*q);

% Assemble the snapshots matrices corresponding to different input
% trajectories 
for j=1:q
    U_hat_s = U_hat(:,(j-1)*Nt+1:j*Nt);
    Uk_hat(:,(j-1)*(Nt-2)+1:j*(Nt-2)) = U_hat_s(:,1:end-2);
    Ukp1_hat(:,(j-1)*(Nt-2)+1:j*(Nt-2)) = U_hat_s(:,2:end-1);
    Ukp2_hat(:,(j-1)*(Nt-2)+1:j*(Nt-2)) = U_hat_s(:,3:end);
end
%%
% Regularization coefficient
lambda = 1e-8;

% Assemble the matrix of the data in the least squares problem
D_hat = [Ukp1_hat; Uk_hat];

disp(['The rank of the matrix D_hat is ' num2str(rank(D_hat))])
disp(['The condition number of the matrix D_hat is ' ...
    num2str(cond(D_hat))])
disp(['The condition number of the matrix D_hat D_hat^T is ' ...
    num2str(cond(D_hat*D_hat'))])

P_hat = zeros(r,2*r);
rhs = D_hat * Ukp2_hat';

% Tikhonov-regularized normal equations
A_reg = D_hat*D_hat' + lambda*eye(2*r);

% Solve the linear systems of normal equations
for m=1:r
    P_hat(m,:) = A_reg \ rhs(:,m);
end

H_hat = P_hat(:,1:r);
R_hat = P_hat(:,r+1:end);

tEnd_AssemblyAndSolving1 = toc(tstart);
%%
% Compute the re-projected trajectories corresponding to different input
% signals
tstart = tic;
U_bar_multiple = zeros(r,size(U_multiple,2));
for j=1:q
    U_s = zeros(size(U,1),size(U,2));
    U_s(:,1:2) = U_multiple(:,(j-1)*Nt+1:(j-1)*Nt+2);
    U_bar_multiple(:,(j-1)*Nt+1:(j-1)*Nt+2) = Phi' * U_s(:,1:2);
    for i=2:size(U,2)-1
        U_s = Phi * (Phi' * U_s);
        U_s(:,i+1) = (M + tau^2*theta*K) \ (M*(2*U_s(:,i)-U_s(:,i-1)) ...
                    + tau^2*K*((2*theta-1)*U_s(:,i) - theta*U_s(:,i-1)) );
        U_bar_multiple(:,(j-1)*Nt+i+1) = Phi' * U_s(:,i+1);
    end
end

U_bar_multiple = U_bar_multiple / norm(U_bar_multiple, 'fro');

% Assemble the snapshots matrices corresponding to different input
% trajectories in the case of re-projection
for j=1:q
    U_bar_s = U_bar_multiple(:,(j-1)*Nt+1:j*Nt);
    Uk_bar(:,(j-1)*(Nt-2)+1:j*(Nt-2)) = U_bar_s(:,1:end-2);
    Ukp1_bar(:,(j-1)*(Nt-2)+1:j*(Nt-2)) = U_bar_s(:,2:end-1);
    Ukp2_bar(:,(j-1)*(Nt-2)+1:j*(Nt-2)) = U_bar_s(:,3:end);
end

% Assemble the matrix of the data in the least squares problem of the 
% re-projection dOpInf
D_bar = [Ukp1_bar; Uk_bar];


disp(['The rank of the matrix D_bar is ' num2str(rank(D_bar))])
disp(['The condition number of the matrix D_bar is ' ...
    num2str(cond(D_bar))])
disp(['The condition number of the matrix D_bar D_bar^T  is ' ...
    num2str(cond(D_bar*D_bar'))])

P_bar = zeros(r,2*r);
rhs2 = D_bar * Ukp2_bar';

% Tikhonov-regularized normal equations
A_reg = D_bar*D_bar' + lambda*eye(2*r);

% Solve the linear systems of normal equations
for m=1:r
    P_bar(m,:) = A_reg \ rhs2(:,m);
end

H_bar = P_bar(:,1:r);
R_bar = P_bar(:,r+1:end);
tEnd_AssemblyAndSolving2 = toc(tstart);
%% Compute the intrusive discrete POD operators
H_r = (M_r + tau^2*theta*K_r) \ (2*M_r - tau^2*(1-2*theta)*K_r);
R_r = -1*eye(r);

%% Compute the relative error of the learned operators 
% without re-projection 
H_hat_error = vecnorm(H_hat(:) - H_r(:),2)/vecnorm(H_r(:),2) * 100;
R_hat_error = vecnorm(R_hat(:) - R_r(:),2)/vecnorm(R_r(:),2) * 100;

fprintf(['Inferred Opeartors Errors:\nH_hat_error: %d, ' ...
    'R_hat_error: %d\n'], H_hat_error, R_hat_error);
%% Compute the relative error of the learned operators 
% with re-projection 
H_bar_error = vecnorm(H_bar(:) - H_r(:),2)/vecnorm(H_r(:),2) * 100;
R_bar_error = vecnorm(R_bar(:) - R_r(:),2)/vecnorm(R_r(:),2) * 100;

fprintf(['Inferred Opeartors Errors:\nH_bar_error: %d, ' ...
    'R_bar_error: %d\n'], H_bar_error, R_bar_error);
%% Compute the POD-G-ROM solution
tstart = tic;
u_r = zeros(r,size(U,2));
u_r(:,1:2) = Phi' * U(:,1:2);

for i=2:size(U,2)-1
    u_r(:,i+1) = ( H_r*u_r(:,i) ...
            + R_r*u_r(:,i-1));
end

U_POD = Phi*u_r;
tEnd_Online1 = toc(tstart);
%% Compute the dOpInf-ROM solution (without re-projection)
tstart = tic;
u_hat = zeros(r,size(U,2));
u_hat(:,1:2) = Phi' * U(:,1:2);

for i=2:size(U,2)-1
    u_hat(:,i+1) = ( H_hat*u_hat(:,i) ...
            + R_hat*u_hat(:,i-1));
end

U_OpInf = Phi*u_hat;
tEnd_Online2 = toc(tstart);
%% Compute the dOpInf-ROM solution (with re-projection)
tstart = tic;
u_bar = zeros(r,size(U,2));
u_bar(:,1:2) = Phi' * U(:,1:2);

for i=2:size(U,2)-1
    u_bar(:,i+1) = ( H_bar*u_bar(:,i) ...
            + R_bar*u_bar(:,i-1));
end

U_OpInf_ReProj = Phi*u_bar;
tEnd_Online3 = toc(tstart);
%%
sensor = 9739 - 400;
% sensor = 100;

U1=U(sensor,:);
U2=U_POD(sensor,:);
U3=U_OpInf(sensor,:);
U4=U_OpInf_ReProj(sensor,:);

U_s_FOM = s_U * U1;
U_s_POD = s_U * U2;
U_s_OpInf = s_U * U3;
U_s_OpInf_ReProj = s_U * U4;


error_s_POD = vecnorm(U_s_POD-U_s_FOM,2)/vecnorm(U_s_FOM,2)*100;
error_s_OpInf = vecnorm(U_s_OpInf-U_s_FOM,2)/vecnorm(U_s_FOM,2)*100;
error_s_OpInf_Reproj = vecnorm(U_s_OpInf_ReProj-U_s_FOM,2)/vecnorm(U_s_FOM,2)*100;

disp(['The relative error in displacment POD-G-ROM approximation' ...
    ' at the sensor is ' num2str(error_s_POD)])
disp(['The relative error in displacment dOpInf-ROM (w/o reprojection) approximation' ...
    ' at the sensor is ' num2str(error_s_OpInf)])
disp(['The relative error in displacment dOpInf-ROM (w reprojection) approximation' ...
    ' at the sensor is ' num2str(error_s_OpInf_Reproj)])



color1 = "#469990"; % Teal
color2 = "#911EB4"; % Strong Purple


figure('units','normalized','outerposition',[0 0 1 1])
plot(t,U_s_FOM,'r');
grid;
hold on
plot(t,U_s_POD,'--b');
hold on
plot(t,U_s_OpInf,'Marker','o','MarkerSize',15,Color=color1,LineStyle="-.");
hold on
plot(t,U_s_OpInf_ReProj,'Marker','*','MarkerSize',15,Color=color2,LineStyle="none");
xlim([0,t(end)]);
legend('FOM',['POD-ROM $r=' num2str(r) '$'],['dOpInf $r=' num2str(r) '$'], ...
    ['dOpInf ReProj $r=' num2str(r) '$'],'Location','best');
change_font_sizes(1,15,40);
change_xaxis(80,1,'$t$ in [s]');
change_yaxis(80,1,'Displacement in [m]');
change_legend(80)
change_line_width(5)
set(gcf,'color','w')
export_fig(['FOM_vs_ROM_r_BasisFromMultiInputData' num2str(r)],'-pdf');

error_u_POD = vecnorm(U_POD - U) ./ vecnorm(U) * 100;
error_u_OpInf = vecnorm(U_OpInf - U) ./ vecnorm(U) * 100;
error_u_OpInf_ReProj = vecnorm(U_OpInf_ReProj - U) ./ vecnorm(U) * 100;


figure('units','normalized','outerposition',[0 0 1 1])
plot(t,error_u_POD,'--b');
grid;
hold on
plot(t,error_u_OpInf,'Marker','o','MarkerSize',15,Color=color1,LineStyle="-.");
hold on
plot(t,error_u_OpInf_ReProj,'Marker','*','MarkerSize',15,Color=color2,LineStyle="none");
xlim([0,t(end)]);
legend(['POD-ROM $r=' num2str(r) '$'],['dOpInf $r=' num2str(r) '$'], ...
    ['dOpInf ReProj $r=' num2str(r) '$'],'Location','best');
change_font_sizes(1,15,40);
change_xaxis(80,1,'$t$ in [s]');
change_yaxis(80,1,'$\epsilon_u(t)$ in $\%$');
change_legend(80)
change_line_width(5)
set(gcf,'color','w')
export_fig(['dOpInf_Reg_' num2str(r)],'-pdf');

%%
% fileName0 = ['FOM_vs_PODGROM_Results_r' num2str(r) '.mat'];
% save(fileName0, "U_s_FOM" ,"U_s_POD","t","error_u_POD");
% %%
% fileName00 = ['FOM_vs_dOpInf_Results_r' num2str(r) '.mat'];
% save(fileName00,"U_s_OpInf","U_s_OpInf_ReProj","error_u_OpInf","error_u_OpInf_ReProj");
% %%
% fileName1 = 'dOpInf_results_different_r.mat';
% 
% if isfile(fileName1)
%     % Load existing data
%     load(fileName1);
% 
%     % Check if the current r value already exists
%     if ~ismember(r, r_vals)
%         % Append and sort
%         r_vals = [r_vals, r];
%         epsilon_s_vals = [epsilon_s_vals, error_s_OpInf];
%         H_hat_error_vals = [H_hat_error_vals, H_hat_error];
%         R_hat_error_vals = [R_hat_error_vals, R_hat_error];
%         cond_DDT_vals = [cond_DDT_vals, cond(D_hat*D_hat')];
% 
%         % Sort by r_vals and reorder epsilon_s_vals accordingly
%         [r_vals, sortIdx] = sort(r_vals);
%         epsilon_s_vals = epsilon_s_vals(sortIdx);
%         H_hat_error_vals = H_hat_error_vals(sortIdx);
%         R_hat_error_vals = R_hat_error_vals(sortIdx);
%         cond_DDT_vals = cond_DDT_vals(sortIdx);
% 
%         % Save updated data
%         save(fileName1, "r_vals", "epsilon_s_vals","H_hat_error_vals", ...
%             "R_hat_error_vals", "cond_DDT_vals");
%     else
%         disp(['Value of r already exists in the OpInf-ROM errors file.' ...
%             ' No update performed.']);
%     end
% 
% else
%     % Create new file if it doesn't exist
%     r_vals = r;
%     epsilon_s_vals = error_s_OpInf;
%     H_hat_error_vals = H_hat_error;
%     R_hat_error_vals = R_hat_error;
%     cond_DDT_vals = cond(D_hat*D_hat');
%     save(fileName1, "r_vals", "epsilon_s_vals","H_hat_error_vals", ...
%             "R_hat_error_vals", "cond_DDT_vals");
% end
% %%
% fileName2 = 'dOpInfReproj_results_different_r.mat';
% 
% if isfile(fileName2)
%     % Load existing data
%     load(fileName2);
% 
%     % Check if the current r value already exists
%     if ~ismember(r, r_vals)
%         % Append and sort
%         r_vals = [r_vals, r];
%         epsilon_s_vals = [epsilon_s_vals, error_s_OpInf_Reproj];
%         H_bar_error_vals = [H_bar_error_vals, H_bar_error];
%         R_bar_error_vals = [R_bar_error_vals, R_bar_error];
%         cond_DDT_vals = [cond_DDT_vals, cond(D_bar*D_bar')];
% 
%         % Sort by r_vals and reorder epsilon_s_vals accordingly
%         [r_vals, sortIdx] = sort(r_vals);
%         epsilon_s_vals = epsilon_s_vals(sortIdx);
%         H_bar_error_vals = H_bar_error_vals(sortIdx);
%         R_bar_error_vals = R_bar_error_vals(sortIdx);
%         cond_DDT_vals = cond_DDT_vals(sortIdx);
% 
%         % Save updated data
%         save(fileName2, "r_vals", "epsilon_s_vals","H_bar_error_vals", ...
%             "R_bar_error_vals", "cond_DDT_vals");
%     else
%         disp(['Value of r already exists in the OpInf-ROM errors file.' ...
%             ' No update performed.']);
%     end
% 
% else
%     % Create new file if it doesn't exist
%     r_vals = r;
%     epsilon_s_vals = error_s_OpInf_Reproj;
%     H_bar_error_vals = H_bar_error;
%     R_bar_error_vals = R_bar_error;
%     cond_DDT_vals = cond(D_bar*D_bar');
%     save(fileName2, "r_vals", "epsilon_s_vals","H_bar_error_vals", ...
%             "R_bar_error_vals", "cond_DDT_vals");
% end