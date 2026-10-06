clear; close all; clc;

dataFile = 'Comsol_data_1001Snapshots.mat';
load(dataFile);

U_COM(1:2,:) = [];
Udot_COM(1:2,:) = [];

N_T = size(U_COM,2);

fex = 120000; T=1/fex; n=5; tE = 5*T; tsim=5*n*T;
Vamp = 10e-9;
UU = zeros(1,N_T);
t = linspace(0,tsim,N_T);
dt=t(2)-t(1);
index_tE = find(t>=tE,1);

for i=1:index_tE
    tt=(i-1)*dt;
    UU(i)=Vamp*sin(2*pi*fex*tt)*sin(pi*fex/n*tt)^2;
end

L0 = zeros(size(U_COM,1),N_T);
L0(2,:) = UU;
L0(966,:) = UU;

[Uddot_COM, ~] = finiteDifference8thOrder(Udot_COM,dt);

tic
% Step 1: Perform POD on displacement data
r = 100;
[U, S, V] = svd(U_COM, 'econ'); % Perform SVD on snapshot matrix u
toc
Phi = U(:, 1:r); % Reduced basis

t(1:4) = [];
t(end-3:end) = [];
U_COM(:,1:4) = [];
U_COM(:,end-3:end) = [];
Udot_COM(:,1:4) = [];
Udot_COM(:,end-3:end) = [];
Uddot_COM(:,1:4) = [];
Uddot_COM(:,end-3:end) = [];
s_U = norm(U_COM,'fro');
s_Udot = norm(Udot_COM,'fro');
s_Uddot = norm(Uddot_COM,'fro');
s_R = norm(UU,'fro');
U_COM = U_COM/s_U;
Udot_COM = Udot_COM/s_Udot;
Uddot_COM = Uddot_COM/s_Uddot;
UU = UU/s_R;

UU(1:4) = [];
UU(end-3:end) = [];

M = eye(r,r);
u_r_0 = Phi'*U_COM(:,1);
udot_r_0 = Phi'*(Udot_COM(:,1));
uddot_r_0 = Phi'*(Uddot_COM(:,1));
%%
tic
U_hat = Phi'*U_COM;
Uddot_hat = Phi'*Uddot_COM;
D_hat = [U_hat; UU];

disp(['The condition number of the matrix D_hat is ' num2str(cond(D_hat))])
disp(['The condition number of the matrix D_hat D_hat^T  is ' ...
num2str(cond(D_hat*D_hat'))])
disp(['The rank of the matrix D_hat is ' num2str(rank(D_hat))])

P_hat = (Uddot_hat*D_hat') / (D_hat*D_hat') ;

K_M_hat = -1*P_hat(:,1:r);
B_M_hat = P_hat(:,r+1:end);
toc

tic
F_M_hat = B_M_hat*UU;

beta = 1/6;
gamma = 0.5;
[u_opinf, ~ , uddot_opinf] = NewmarkIntegrationLinear(M/s_Uddot,...
K_M_hat/s_U,F_M_hat,dt,r,tsim-8*dt,1,[]...
,u_r_0*s_U,udot_r_0*s_Udot,beta,gamma);

U_OpInf = Phi*u_opinf/s_U;
toc

%%
sensor = 15648;
U1 = U_COM(sensor,:);
U2 = U_OpInf(sensor,:);
error = vecnorm(U2-U1,2)/vecnorm(U1,2)*100;
disp(['The relative error in displacment ROM approximation at the sensor is ' num2str(error)])

figure('units','normalized','outerposition',[0 0 1 1])
plot(t,U1*s_U,'r');
grid;
hold on
plot(t,U2*s_U,'--b');
% title(['Model with damage: FOM and ROM results at the spatial point number ' num2str(sensor)]);
line = xline(tE,'-',{'Excitation Time'});
line.LabelHorizontalAlignment = 'left';
line.LineWidth = 5;
line.FontSize = 80;
xlim([0 tsim])
xticks([5e-5 1e-4 1.5e-4 2e-4 tsim]);
legend('FOM',['OpInf ROM $r=' num2str(r) '$'],'Location','best');
change_font_sizes(1,15,40);
change_xaxis(80,1,'$t$ in [s]');
change_yaxis(80,1,'Out-of-plane displacement in [m]');
change_legend(80)
change_line_width(5)
set(gcf,'color','w')
export_fig(['opinf_romVScomsol_r' num2str(r)],'-pdf');

pause(1);

error_Opinf_Direct = vecnorm(U_OpInf-U_COM,2)./vecnorm(U_COM,2)*100;

beta = 1/4;
[u_opinf2, ~ , uddot_opinf2] = NewmarkIntegrationLinear(M/s_Uddot,...
K_M_hat/s_U,F_M_hat,dt,r,tsim-8*dt,1,[]...
,u_r_0*s_U,udot_r_0*s_Udot,beta,gamma);


U_OpInf2 = Phi*u_opinf2/s_U;
error_Opinf_Direct2 = vecnorm(U_OpInf2-U_COM,2)./vecnorm(U_COM,2)*100;

figure('units','normalized','outerposition',[0 0 1 1])
plot(t,error_Opinf_Direct,'r');
grid;
hold on
plot(t,error_Opinf_Direct2,'--b');
legend('Linear accelearation','Constant accelearation','Location','best');
xlim([tsim/30 tsim])
change_font_sizes(1,15,40);
change_xaxis(80,1,'$t$ in [s]');
change_yaxis(80,1,'$\epsilon_u$ in $\%$');
change_legend(80)
change_line_width(5)
set(gcf,'color','w')
export_fig(['error_u_r' num2str(r) '_opinf'],'-pdf');

%%
fileName1 = 'OpInf_results_different_r.mat';

if isfile(fileName1)
    % Load existing data
    load(fileName1);
    
    % Check if the current r value already exists
    if ~ismember(r, r_vals)
        % Append and sort
        r_vals = [r_vals, r];
        epsilon_s_vals = [epsilon_s_vals, error];
      
        % Sort by r_vals and reorder epsilon_s_vals accordingly
        [r_vals, sortIdx] = sort(r_vals);
        epsilon_s_vals = epsilon_s_vals(sortIdx);
        
        % Save updated data
        save(fileName1, "r_vals", "epsilon_s_vals");
    else
        disp(['Value of r already exists in the OpInf-ROM errors file.' ...
            ' No update performed.']);
    end
    
else
    % Create new file if it doesn't exist
    r_vals = r;
    epsilon_s_vals = error;
    save(fileName1, "r_vals", "epsilon_s_vals");
end
