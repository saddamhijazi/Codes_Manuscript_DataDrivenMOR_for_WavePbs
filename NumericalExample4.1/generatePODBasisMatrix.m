clear; clc; close all;

% Damage Parameters
xD = 3;  % Position 
yD = 3; 
cD = 0.05; % Wave speed in the damage area
dD = 0.06; % Damage width 
params = [xD,yD,cD,dD];

T = 5; 
tau = 0.02; 
dx = 0.05;

nmb_timestep = T/tau;
theta = 0.25;
[nodes,mesh,model] = gitter(dx);
nmb_nodes = size(nodes,2);
nE = findNodes(mesh,"region","Edge",[1 2 3 4]);
nmb_eNodes = size(nE,2);
B = constructB(nmb_nodes,nmb_eNodes); 

compute_sol = 1;

if(compute_sol)
    tic
    % Solve the FOM
    [U,M,K] = solveFOM(params,nodes,mesh,model,tau,T);
    toc
    save('Snapshots.mat','U','M','K','B');
else
    load('Snapshots.mat');
end

[~, Uddot] = finiteDifference8thOrder(U, tau);

tsim = T;
nSkip = 0.1/tau;
dt = tau;
t = nSkip*dt:dt:tsim-dt;

umax = max(max(U));
umin = min(min(U));

%%
figure('units','normalized','outerposition',[0 0 1 1])
pdeplot(mesh,XYData=B*U(:,5))%,ColorMap=summer)
axis([0 5 0 5])
ax = gca; ax.XTick = [0 1 2 3 4 5]; ax.YTick = [0 1 2 3 4 5]; ax.FontSize = 24;
grid on 
grid minor
cb = colorbar;
cb.TickLabelInterpreter = 'latex';
change_font_sizes(1,15,80);
change_xaxis(80,1,'$x_1$');
change_yaxis(80,1,'$x_2$');
set(gcf,'color','w')
ax = gca;
ax.Units = 'normalized';
% ax.Position = [0.13 0.11 0.775 0.815];  % standard MATLAB axes
exportgraphics(gcf, 'u0.pdf');

V = 1 / tau * ( U(:,2:end) - U(:,1:end-1) );

figure('units','normalized','outerposition',[0 0 1 1])
pdeplot(mesh,XYData=B*V(:,4))%,ColorMap=summer)
axis([0 5 0 5])
ax = gca; ax.XTick = [0 1 2 3 4 5]; ax.YTick = [0 1 2 3 4 5]; ax.FontSize = 24;
grid on 
grid minor
cb = colorbar;
cb.TickLabelInterpreter = 'latex';
change_font_sizes(1,15,80);
change_xaxis(80,1,'$x_1$');
change_yaxis(80,1,'$x_2$');
set(gcf,'color','w')
ax = gca;
ax.Units = 'normalized';
% ax.Position = [0.13 0.11 0.775 0.815];  % standard MATLAB axes
exportgraphics(gcf, 'v0.pdf');
%%

figure('units','normalized','outerposition',[0 0 1 1])
pdeplot(mesh,XYData=B*U(:,1))%,ColorMap=summer)
axis([0 5 0 5 umin umax]); 
ax = gca; ax.XTick = [0 1 2 3 4 5]; ax.YTick = [0 1 2 3 4 5]; ax.FontSize = 24;
grid on 
grid minor
cb = colorbar;
cb.TickLabelInterpreter = 'latex';
change_font_sizes(1,15,80);
change_xaxis(80,1,'$x_1$');
change_yaxis(80,1,'$x_2$');
set(gcf,'color','w')
ax = gca;
ax.Units = 'normalized';
% ax.Position = [0.13 0.11 0.775 0.815];  % standard MATLAB axes
exportgraphics(gcf, 'uFOMt0.pdf');

figure('units','normalized','outerposition',[0 0 1 1])
pdeplot(mesh,XYData=B*U(:,125))%,ColorMap=summer)
axis([0 5 0 5 umin umax]); 
ax = gca; ax.XTick = [0 1 2 3 4 5]; ax.YTick = [0 1 2 3 4 5]; ax.FontSize = 24;
grid on 
grid minor
cb = colorbar;
cb.TickLabelInterpreter = 'latex';
change_font_sizes(1,15,80);
change_xaxis(80,1,'$x_1$');
change_yaxis(80,1,'$x_2$');
set(gcf,'color','w')
ax = gca;
ax.Units = 'normalized';
% ax.Position = [0.13 0.11 0.775 0.815];  % standard MATLAB axes
exportgraphics(gcf, 'uFOMt25.pdf');

figure('units','normalized','outerposition',[0 0 1 1])
pdeplot(mesh,XYData=B*U(:,end))%,ColorMap=summer)
axis([0 5 0 5 umin umax]); 
ax = gca; ax.XTick = [0 1 2 3 4 5]; ax.YTick = [0 1 2 3 4 5]; ax.FontSize = 24;
grid on 
grid minor
cb = colorbar;
cb.TickLabelInterpreter = 'latex';
change_font_sizes(1,15,80);
change_xaxis(80,1,'$x_1$');
change_yaxis(80,1,'$x_2$');
set(gcf,'color','w')
ax = gca;
ax.Units = 'normalized';
% ax.Position = [0.13 0.11 0.775 0.815];  % standard MATLAB axes
exportgraphics(gcf, 'uFOMt5.pdf');

%%

s_U = norm(U,'fro');
s_Uddot = norm(Uddot,'fro');

U = U / s_U;
U = U(:,nSkip+1:end);

save('ScalingVals.mat','s_U','s_Uddot');


r = 200;
[V,S,~] = svd(U,'econ');
Phi = V(:,1:r);
save('Phi_POD.mat','Phi');

figure('units','normalized','outerposition',[0 0 1 1])
semilogy(1:1:size(diag(S)),diag(S),'LineWidth',4);
axis([0 200 10^-20 10^0])
grid;
change_font_sizes(1,15,40);
change_xaxis(80,1,'$r$');
change_yaxis(80,1,'$\sigma_r$');
change_line_width(5)
set(gcf,'color','w')
export_fig('Singular_Values','-pdf');


figure('units','normalized','outerposition',[0 0 1 1])
plot(1:1:size(diag(S)),cumsum(diag(S))/sum(diag(S)),'LineWidth',4);
xlim([1 200])
grid;
change_font_sizes(1,15,40);
change_xaxis(80,1,'$r$');
change_yaxis(80,1,'$c_r$');
change_line_width(5)
set(gcf,'color','w')
export_fig('CumulativeSum_Values','-pdf');