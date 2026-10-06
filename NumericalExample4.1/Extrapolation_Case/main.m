clc; clear; close all;

% Damage Parameters
xD = 3;  % Position 
yD = 3; 
cD = 0.05; % Wave speed in the damage area
dD = 0.06; % Damage width 
params = [xD,yD,cD,dD];

T = 40; 
tau = 0.2; 
dx = 0.05;

nmb_timestep = T/tau; % Number of time steps
theta = 0.25;
[nodes,mesh,model] = gitter(dx);
nmb_nodes = size(nodes,2);
nE = findNodes(mesh,"region","Edge",[1 2 3 4]);
nmb_eNodes = size(nE,2);

compute_sol = 1;

if(compute_sol)
    tstart = tic;
    % Solve the full order model
    [U,M,K,F] = solveFOM(params,nodes,mesh,model,tau,T);
    tEnd = toc(tstart);
    save('Snapshots.mat','U','M','K','F');
else
    load('Snapshots.mat');
end


L = tau^2*(theta*F(:,3:end) + ...
    (1-2*theta)*F(:,2:end-1) + theta*F(:,1:end-2));
s_L = norm(L,'fro');
L = L / s_L;

t = 0:tau:T-tau;

s_U = norm(U,'fro');

U = U / s_U;

N = size(U,1);
Nt = size(U,2);

% Define several positions of the actuators to simulate the FOM for
% different trajectories in order to assemble a regression matrix which has
% modest condition number
xy_c = [0 5; 5 0; 0 0; 5 5; 2.5 2.5];
q = size(xy_c,1);
U_multiple = zeros(N,q*Nt);
F_multiple = zeros(N,q*Nt);

pars_c = zeros(q,4);
pars_c(:,1) = 1.0;
pars_c(:,4) = 2.0;
pars_c(:,5) = 0.1;
pars_c(:,6) = 0.2;


[U_POD, S, V] = svd(U, 'econ'); % Perform SVD on snapshot matrix u
r = 80;
Phi = U_POD(:, 1:r); % Reduced basis

Uk_hat = zeros(r,(Nt-2)*q);
Ukp1_hat = zeros(r,(Nt-2)*q);
Ukp2_hat = zeros(r,(Nt-2)*q);
L_r_mult = zeros(r,(Nt-2)*q);

tstart = tic;
for j=1:q
    pars_c(j,2:3) = xy_c(j,:);

    % Folder and filename
    folderName = sprintf('T%d', T);
    
    if ~exist(folderName, 'dir')
        mkdir(folderName);
    end
    
    fileName = sprintf('U_s_%ds_xc%.1f_yc%.1f.mat', ...
    T, xy_c(j,1), xy_c(j,2));
    
    filePath = fullfile(folderName, fileName);
    
    % Load if already computed
    if exist(filePath, 'file')
        data = load(filePath, 'U_s', 'F_s');
        U_s = data.U_s;
        F_s = data.F_s;
    % Otherwise, compute the FOM corresponding to the actuator parameters
    else
        [U_s,~,~,F_s] = solveFOM_Custom( ...
            params, nodes, mesh, model, ...
            tau, T, pars_c(j,:));
    
        save(filePath, 'U_s', 'F_s');
    end
    U_multiple(:,(j-1)*Nt+1:j*Nt) = U_s;
    F_multiple(:,(j-1)*Nt+1:j*Nt) = F_s;
    Uk_hat(:,(j-1)*(Nt-2)+1:j*(Nt-2)) = Phi' * U_s(:,1:end-2);
    Ukp1_hat(:,(j-1)*(Nt-2)+1:j*(Nt-2)) = Phi' * U_s(:,2:end-1);
    Ukp2_hat(:,(j-1)*(Nt-2)+1:j*(Nt-2)) = Phi' * U_s(:,3:end);
    L_s = tau^2*(theta*F_s(:,3:end) + ...
    (1-2*theta)*F_s(:,2:end-1) + theta*F_s(:,1:end-2));
    L_r_mult(:,(j-1)*(Nt-2)+1:j*(Nt-2)) = Phi' * L_s;
end
tEnd_multiTraj = toc(tstart);

s_U_multi = norm(U_multiple,'fro');

U_multiple = U_multiple / s_U_multi;
Uk_hat = Uk_hat / s_U_multi;
Ukp1_hat = Ukp1_hat / s_U_multi;
Ukp2_hat = Ukp2_hat / s_U_multi;
L_r_mult = L_r_mult / s_L;
%%
M_r = Phi' * M * Phi;
K_r = Phi' * K * Phi;
L_r = Phi' * L;

A_r = (M_r + tau^2*theta*K_r);
B_r = (2*M_r - tau^2*(1-2*theta)*K_r);
C_r = -1*A_r;
%% Solve the dOpInf regression problem without the re-projection method
D_hat = [Ukp2_hat; Ukp1_hat; Uk_hat];

disp(['The rank of the matrix D_hat is ' num2str(rank(D_hat))])
disp(['The condition number of the matrix D_hat is ' ...
    num2str(cond(D_hat))])
disp(['The condition number of the matrix D_hat D_hat^T  is ' ...
    num2str(cond(D_hat*D_hat'))])

P_hat = zeros(r,3*r);
rhs = D_hat * L_r_mult';

for m=1:r
    P_hat(m,:) = (D_hat*D_hat') \ (rhs(:,m));
end

A_hat = P_hat(:,1:r);
B_hat = -1 * P_hat(:,r+1:2*r);
C_hat = -1 * P_hat(:,2*r+1:end);
%% Solve the dOpInf regression problem with the re-projection method
tstart = tic;
U_bar_multiple = zeros(r,size(U_multiple,2));
for j=1:q
    U_s = zeros(size(U,1),size(U,2));
    U_s(:,1:2) = U_multiple(:,(j-1)*Nt+1:(j-1)*Nt+2);
    F_s = F_multiple(:,(j-1)*Nt+1:j*Nt);
    U_bar_multiple(:,(j-1)*Nt+1:(j-1)*Nt+2) = Phi' * U_s(:,1:2);
    for i=2:size(U,2)-1
        U_s = Phi * (Phi' * U_s);
        U_s(:,i+1) = (M + tau^2*theta*K) \ (M*(2*U_s(:,i)-U_s(:,i-1)) ...
                + tau^2*K*((2*theta-1)*U_s(:,i) - theta*U_s(:,i-1)) ...
                + tau^2*(theta*F_s(:,i+1) + (1-2*theta)*F_s(:,i) + theta*F_s(:,i-1)) );
        U_bar_multiple(:,(j-1)*Nt+i+1) = Phi' * U_s(:,i+1);
    end
end

s_U_mulit_t = norm(U_bar_multiple, 'fro');

U_bar_multiple = U_bar_multiple / s_U_mulit_t;

for j=1:q
    U_bar_s = U_bar_multiple(:,(j-1)*Nt+1:j*Nt);
    Uk_bar(:,(j-1)*(Nt-2)+1:j*(Nt-2)) = U_bar_s(:,1:end-2);
    Ukp1_bar(:,(j-1)*(Nt-2)+1:j*(Nt-2)) = U_bar_s(:,2:end-1);
    Ukp2_bar(:,(j-1)*(Nt-2)+1:j*(Nt-2)) = U_bar_s(:,3:end);
end

D_bar = [Ukp2_bar; Ukp1_bar; Uk_bar];

disp(['The rank of the matrix D_bar is ' num2str(rank(D_bar))])
disp(['The condition number of the matrix D_bar is ' ...
    num2str(cond(D_bar))])
disp(['The condition number of the matrix D_bar D_bar^T  is ' ...
    num2str(cond(D_bar*D_bar'))])


P_bar = zeros(r,3*r);
rhs2 = D_bar * L_r_mult';

for m=1:r
    P_bar(m,:) = (D_bar*D_bar') \ (rhs2(:,m));
end

A_bar = P_bar(:,1:r);
B_bar = -1 * P_bar(:,r+1:2*r);
C_bar = -1 * P_bar(:,2*r+1:end);
tEnd_AssemblyAndSolving2 = toc(tstart);
%%
A_error = vecnorm(B_bar(:) - A_r(:),2)/vecnorm(A_r(:),2) * 100;
B_error = vecnorm(B_bar(:) - B_r(:),2)/vecnorm(B_r(:),2) * 100;
C_error = vecnorm(C_bar(:) - C_r(:),2)/vecnorm(C_r(:),2) * 100;

fprintf(['Inferred Opeartors Errors:\nA_error: %d, ' ...
    'B_error: %d, C_error: %d\n'], A_error, B_error, C_error); 
%% Solve the POD-ROM
tstart = tic;
u_r = zeros(r,size(U,2));
u_r(:,1:2) = Phi' * U(:,1:2);
Frop_f = Phi' * F;

for i=2:size(U,2)-1
    u_r(:,i+1) = A_r \ ( B_r*u_r(:,i) ...
            + C_r*u_r(:,i-1) ...
            + L_r(:,i-1) / s_U * s_L );
end

U_POD = Phi*u_r;
tEnd_Online1 = toc(tstart);
%% Solve the dOpInf ROM w/o re-projection
tstart = tic;
u_hat = zeros(r,size(U,2));
u_hat(:,1:2) = Phi' * U(:,1:2);

for i=2:size(U,2)-1
    u_hat(:,i+1) = A_hat \ ( B_hat*u_hat(:,i) ...
            + C_hat*u_hat(:,i-1) ...
            + L_r(:,i-1) );
end

U_OpInf = Phi*u_hat*s_U_multi/s_U;
tEnd_Online2 = toc(tstart);
%% Solve the dOpInf ROM w re-projection
tstart = tic;
u_bar = zeros(r,size(U,2));
u_bar(:,1:2) = Phi' * U(:,1:2);

for i=2:size(U,2)-1  
    u_bar(:,i+1) = A_bar \ ( B_bar*u_bar(:,i) ...
            + C_bar*u_bar(:,i-1) ...
            + L_r(:,i-1) );
end

U_OpInf_ReProj = Phi*u_bar*s_U_mulit_t/s_U;
tEnd_Online3 = toc(tstart);
%%

sensor = 2343;

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

color2 = "#911EB4"; % Strong Purple

%% Compute the solution for the extrapolation window
tsim = 120;
t_long = 0:tau:tsim-tau;

% Folder and filename
    folderName = sprintf('T%d', tsim);
    
    if ~exist(folderName, 'dir')
        mkdir(folderName);
    end
    
    fileName = sprintf('LongTimeData_T%ds.mat', tsim);
    
    filePath = fullfile(folderName, fileName);
    
    % Load if already computed
    if exist(filePath, 'file')
        data = load(filePath, 'U_longTime', 'Fhist_long');
        U_longTime = data.U_longTime;
        Fhist_long = data.Fhist_long;
    else
        t_start = tic;
        % Solve the FOM
        [U_longTime,~,~,Fhist_long] = solveFOM(params,nodes,mesh,model,tau,tsim);
        t_End_long = toc(t_start);
    
        save(filePath, 'U_longTime', 'Fhist_long','t_End_long');
    end

U_longTime = U_longTime / s_U;

L = tau^2*(theta*Fhist_long(:,3:end) + ...
    (1-2*theta)*Fhist_long(:,2:end-1) + theta*Fhist_long(:,1:end-2));
L_r_long = Phi' * L;
L_r_long = L_r_long / s_L;
%% Solve the POD-ROM
u_POD = zeros(r,size(U_longTime,2));
u_POD(:,1:2) = Phi'*U_longTime(:,1:2);

for i=2:size(U_longTime,2)-1
    u_POD(:,i+1) = A_r \ ( B_r*u_POD(:,i) ...
            + C_r*u_POD(:,i-1) ...
            + L_r_long(:,i-1) / s_U * s_L );
end

U_POD = Phi*u_POD;
%% Solve the dOpInf ROM w/o re-projection
u_hat = zeros(r,size(U_longTime,2));
u_hat(:,1:2) = Phi' * U_longTime(:,1:2);

for i=2:size(U_longTime,2)-1
    u_hat(:,i+1) = A_hat \ ( B_hat*u_hat(:,i) ...
            + C_hat*u_hat(:,i-1) ...
            + L_r_long(:,i-1) );
end

U_OpInf = Phi*u_hat*s_U_multi/s_U;
%% Solve the dOpInf ROM w re-projection
u_bar = zeros(r,size(U_longTime,2));
u_bar(:,1:2) = Phi' * U_longTime(:,1:2);

for i=2:size(U_longTime,2)-1  
    u_bar(:,i+1) = A_bar \ ( B_bar*u_bar(:,i) ...
            + C_bar*u_bar(:,i-1) ...
            + L_r_long(:,i-1) );
end

U_OpInf_ReProj = Phi*u_bar*s_U_mulit_t/s_U;
%% Plot the results
sensor = 2343;

U1_long=U_longTime(sensor,:);
U2_long=U_POD(sensor,:);
U3_long=U_OpInf(sensor,:);
U4_long=U_OpInf_ReProj(sensor,:);

U_s_FOM = s_U * U1_long;
U_s_POD = s_U * U2_long;
U_s_OpInf = s_U * U3_long;
U_s_OpInf_ReProj = s_U * U4_long;

error = vecnorm(U2_long-U1_long,2)/vecnorm(U1_long,2)*100;
disp(['The relative error in displacment ROM approximation at the sensor is ' num2str(error)])

figure('units','normalized','outerposition',[0 0 1 1])
plot(t_long,U_s_FOM,'r','MarkerSize',15);
grid;
hold on
plot(t_long,U_s_POD,'--b','MarkerSize',15);
hold on
plot(t_long,U_s_OpInf_ReProj,'Marker','*','MarkerSize',15,Color=color2,LineStyle="none");
line = xline(T,'-',{'Training Window'});
line.LabelHorizontalAlignment = 'left';
line.LineWidth = 5;
line.FontSize = 80;
xlim([0,t_long(end)+tau]);
legend('FOM', ['POD-ROM $r=' num2str(r) '$'] , ['dOpInf ReProj $r=' num2str(r) '$'] ...
    ,'Location','best');
change_font_sizes(1,15,40);
change_xaxis(80,1,'$t$ in [s]');
change_yaxis(80,1,'Displacement in [m]');
change_legend(80)
change_line_width(5)
set(gcf,'color','w')
export_fig(['FOM_vs_ROMs_Extrapolation_ComparisonAtSpatialPoint' num2str(sensor)],'-pdf');

% Compute the L^2 spatial relative error as function of time for all ROMs
error_u_POD_long = vecnorm(U_POD-U_longTime,2)./vecnorm(U_longTime,2) * 100;
error_u_OpInf_long = vecnorm(U_OpInf - U_longTime) ./ vecnorm(U_longTime) * 100;
error_u_OpInf_ReProj_long = vecnorm(U_OpInf_ReProj - U_longTime) ./ vecnorm(U_longTime) * 100;

pause(1);

figure('units','normalized','outerposition',[0 0 1 1])
plot(t_long,error_u_POD_long,'--b');
grid;
hold on
plot(t_long,error_u_OpInf_ReProj_long,'Marker','*','MarkerSize',15,Color=color2,LineStyle="none");
line = xline(T,'-',{'Training Window'});
line.LabelHorizontalAlignment = 'left';
line.LineWidth = 5;
line.FontSize = 80;
xlim([0,t_long(end) + tau]);
lgd = legend(['POD-ROM $r=' num2str(r) '$'] , ['dOpInf ReProj $r=' num2str(r) '$'] ...
    ,'Location','best');
lgd.Position = [0.55 0.70 0.20 0.15];
change_font_sizes(1,15,40);
change_xaxis(80,1,'$t$ in [s]');
change_yaxis(80,1,'$\epsilon_u(t)$ in $\%$');
change_legend(80)
change_line_width(5)
set(gcf,'color','w')
export_fig('FOM_vs_ROMs_Extrapolation_Spatial_Error','-pdf');



extrapolation_mean_error = mean(error_u_POD_long(size(U,2)+1:end));

%% Save the CPU times for the record
if(~isfile("CPU_Times.mat"))
    save('CPU_Times','tEnd_multiTraj','tEnd_AssemblyAndSolving2', ...
        'tEnd','tEnd_Online1','tEnd_Online2','tEnd_Online3')
end