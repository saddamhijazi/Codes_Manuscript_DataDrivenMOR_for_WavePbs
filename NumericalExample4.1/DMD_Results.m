clear; close all; clc;

load('Snapshots.mat')

tau = 0.02;
T = 5;

t = 0:tau:T-tau;
t(1:floor(0.1/tau)) = [];
U(:,1:floor(0.1/tau)) = [];
r = 100;

Y = standard_DMD(U,r);

t_start_DMD = tic;
Xdmd_opt = optimal_dmd(U, r);
t_end_DMD = toc(t_start_DMD);

%%
sensor = 9339;
U_s_FOM = U(sensor,:);
U_s_DMD = Y(sensor,:);
U_s_DMD_Opt = Xdmd_opt(sensor,:);

color1 = "#7E2F8E";

figure('units','normalized','outerposition',[0 0 1 1])
plot(t,U(sensor,:),'r');
grid;
hold on
plot(t,Y(sensor,:),'--b');
hold on
plot(t,Xdmd_opt(sensor,:),'Marker','o',Color=color1,LineStyle="-");
xlim([0,t(end)]);
legend('FOM','DMD','DMD Optimal','Location','best');
change_font_sizes(1,15,40);
change_xaxis(60,1,'$t$ in [s]');
change_yaxis(60,1,'The displacement in [m]');
change_legend(60)
change_line_width(2.5)
set(gcf,'color','w')
% export_fig(['DMD_Comparison_Classical_vs_OptimalDMD_r' num2str(r)],'-pdf');

error_s_DMD = vecnorm( U_s_FOM - U_s_DMD , 2 ) / vecnorm(U_s_FOM,2) * 100;
disp(['The relative error in displacment DMD approximation at the sensor is ' num2str(error_s_DMD)])

error_s_DMDOpt = vecnorm( U_s_FOM - U_s_DMD_Opt , 2 ) / vecnorm(U_s_FOM,2) * 100;
disp(['The relative error in displacment DMD Optimal approximation at the sensor is ' num2str(error_s_DMDOpt)])

epsilon_u_DMD = vecnorm(Y-U,2)./vecnorm(U,2)*100;
epsilon_u_DMDOpt = vecnorm(Xdmd_opt-U,2)./vecnorm(U,2)*100;

fileName = ['FOM_vs_DMD_Results_r' num2str(r) '.mat'];
save(fileName, "U_s_FOM" ,"U_s_DMD","U_s_DMD_Opt", ...
    "t","epsilon_u_DMD","epsilon_u_DMDOpt");

fileName1 = 'DMD_Epsilon_s_vs_r.mat';

if isfile(fileName1)
    % Load existing data
    load(fileName1, "r_vals", "epsilon_s_vals");
    
    % Check if the current r value already exists
    if ~ismember(r, r_vals)
        % Append and sort
        r_vals = [r_vals, r];
        epsilon_s_vals = [epsilon_s_vals, error_s_DMD];
        
        % Sort by r_vals and reorder epsilon_s_vals accordingly
        [r_vals, sortIdx] = sort(r_vals);
        epsilon_s_vals = epsilon_s_vals(sortIdx);
        
        % Save updated data
        save(fileName1, "r_vals", "epsilon_s_vals");
    else
        disp(['Value of r already exists in the DMD errors file.' ...
            ' No update performed.']);
    end
    
else
    % Create new file if it doesn't exist
    r_vals = r;
    epsilon_s_vals = error_s_DMD;
    save(fileName1, "r_vals", "epsilon_s_vals");
end


fileName2 = 'DMD_Opt_Epsilon_s_vs_r.mat';

if isfile(fileName2)
    % Load existing data
    load(fileName2, "r_vals", "epsilon_s_vals");
    
    % Check if the current r value already exists
    if ~ismember(r, r_vals) && sum(Xdmd_opt(:)) ~=0
        % Append and sort
        r_vals = [r_vals, r];
        epsilon_s_vals = [epsilon_s_vals, error_s_DMDOpt];
        
        % Sort by r_vals and reorder epsilon_s_vals accordingly
        [r_vals, sortIdx] = sort(r_vals);
        epsilon_s_vals = epsilon_s_vals(sortIdx);
        
        % Save updated data
        save(fileName2, "r_vals", "epsilon_s_vals");
    else
        disp(['Value of r already exists in DMD-Opt errors file or ill-conditioned problem encountered.' ...
            ' No update performed.']);
    end
    
else
    % Create new file if it doesn't exist
    r_vals = r;
    epsilon_s_vals = error_s_DMDOpt;
    save(fileName2, "r_vals", "epsilon_s_vals");
end
