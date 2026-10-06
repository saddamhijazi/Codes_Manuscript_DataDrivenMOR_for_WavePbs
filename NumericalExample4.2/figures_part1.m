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


r = 80;

s = 1;
t_start = t(s);
t = t(s:end);
U_COM = U_COM(:,s:end);

Y = standard_DMD(U_COM,r);

t_start_DMD = tic;
Xdmd_opt = optimal_dmd(U_COM, r);
t_end_DMD = toc(t_start_DMD);


t_start_mrdmd = tic;
tic
nodes = mrdmd(U_COM,0,0,0,6,1,false);
toc

tic
X_mrdmd = 0*U_COM;
for i=0:7
    [Phi,Psi] = stitch_mrdmd(nodes,i);
    X_mrdmd = X_mrdmd + Phi*Psi;
end
toc

X_mrdmd = real(X_mrdmd);
t_end_mrdmd = toc(t_start_mrdmd);

%%

figure('units','normalized','outerposition',[0 0 1 1])
plot(t,UU,'r');
grid;
xlim([0,t(end)]);
xticks([5e-5 1e-4 1.5e-4 2e-4 tsim]);
change_font_sizes(1,15,40);
change_xaxis(80,1,'$t$ in [s]');
change_yaxis(80,1,'$f(t)$ in [N]');
change_line_width(5)
set(gcf,'color','w')
export_fig('Force_5CycleHanningWindow','-pdf');

%%
level_test = 2;
[Phi2,Psi2] = stitch_mrdmd(nodes,level_test);

n = 1:size(Psi2,1);

if(length(n)>1)
    figure('units','normalized','outerposition',[0 0 1 1])
    plot_styles(t,real(Psi2),'Markersize',5);
    grid;
    xlim([0,t(end)]);
    xticks(0:2e-5:tsim);
    legend(cellstr(num2str(n', 'n=%-d')),'Location','best');
    change_font_sizes(1,15,40);
    change_xaxis(80,1,'$t$ in [s]');
    change_yaxis(80,1,'$\Psi$');
    change_legend(80)
    change_line_width(5)
    set(gcf,'color','w')
    export_fig(strcat('DMD_timeCoefficients_FML_DamagedModel_level_',num2str(level_test)),'-pdf');
end

%%

color1 = "#7E2F8E";
Spacepoint=15648;

figure('units','normalized','outerposition',[0 0 1 1])
plot(t,U_COM(Spacepoint,:),'r');
grid;
hold on
plot(t,Y(Spacepoint,:),'--b');
hold on
plot(t,Xdmd_opt(Spacepoint,:),'Marker','o','MarkerSize',15,Color=color1,LineStyle="None");
xlim([0,t(end)]);
xticks(0:2e-5:tsim);
legend('FOM','DMD','DMD Optimal','Location','best');
change_font_sizes(1,15,40);
change_xaxis(80,1,'$t$ in [s]');
change_yaxis(80,1,'Out-of-plane displacement in [m]');
change_legend(80)
change_line_width(5)
set(gcf,'color','w')
export_fig(['DMD_Comparison_Classical_vs_OptimalDMD_r' num2str(r)],'-pdf');

figure('units','normalized','outerposition',[0 0 1 1])
plot(t,U_COM(Spacepoint,:),'r');
grid;
hold on
plot(t,X_mrdmd(Spacepoint,:),'*k','MarkerSize',15);
xlim([0,t(end)]);
xticks(0:2e-5:tsim);
c=legend('FOM','mrDMD','Location','best');
change_font_sizes(1,15,40);
change_xaxis(80,1,'$t$ in [s]');
change_yaxis(80,1,'Out-of-plane displacement in [m]');
change_legend(80)
change_line_width(5)
set(gcf,'color','w')
export_fig(strcat('FOM_vs_mrDMD_FML_DamagedModel_atPoint', num2str(Spacepoint)),'-pdf');
%%
[~, S , ~ ] = svd(U_COM,'econ');

% plot the singular values in semi-log scale
figure('units','normalized','outerposition',[0 0 1 1])
semilogy(diag(S))
xlim([1 1001])
grid
change_font_sizes(1,15,40);
change_xaxis(80,1,'$r$');
change_yaxis(80,1,'$\sigma_r$');
change_line_width(5)
set(gcf,'color','w')
export_fig('decay_SV_nonPar','-pdf');



sv = diag(S);

% plot cumulative sum of singular values
figure('units','normalized','outerposition',[0 0 1 1])
plot(cumsum(sv)/sum(sv),'LineWidth',2);
xlim([1 1001])
change_font_sizes(1,15,40);
change_xaxis(80,1,'$r$');
change_yaxis(80,1,'$c_r$');
change_line_width(5)
set(gcf,'color','w')
export_fig('CumSumSV_nonPar','-pdf');
