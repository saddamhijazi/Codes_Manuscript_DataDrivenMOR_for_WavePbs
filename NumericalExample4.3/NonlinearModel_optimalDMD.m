clear; close all; clc;

data = load('NonlinearAluminiumModel_Comsol_Data.mat');

U_FOM = data.U;
fex = 100000; T=1/fex; n=5; dt=T/20; tE = n*T; tsim=5*n*T;
Vamp = 10e-9;
t = linspace(0,tsim,size(U_FOM,2));

r = 80;
Y = standard_DMD(U_FOM,r);

Xdmd_opt = optimal_dmd(U_FOM, r);

Spacepoint = 15642;

%%
color1 = "#7E2F8E";

figure('units','normalized','outerposition',[0 0 1 1])
plot(t,U_FOM(Spacepoint,:),'r');
grid;
hold on
plot(t,Y(Spacepoint,:),'--b');
hold on
plot(t,Xdmd_opt(Spacepoint,:),'Marker','o',Color=color1,LineStyle="-");
xlim([0,t(end)]);
xticks([5e-5 1e-4 1.5e-4 2e-4 tsim]);
legend('FOM','DMD','DMD Optimal','Location','best');
change_font_sizes(1,15,40);
change_xaxis(80,1,'$t$ in [s]');
change_yaxis(80,1,'Out-of-plane displacement in [m]');
change_legend(80)
change_line_width(5)
set(gcf,'color','w')
export_fig(['DMD_Comparison_Classical_vs_OptimalDMD_r' num2str(r) '_NonLinearAluModel'],'-pdf');

error = vecnorm( U_FOM(Spacepoint,:) - Xdmd_opt(Spacepoint,:) , 2 ) / vecnorm(U_FOM(Spacepoint,:),2) * 100;
disp(['The relative error in displacment DMD Optimal approximation at the sensor is ' num2str(error)])
