clear; close all; clc;

data = load('NonlinearAluminiumModel_Comsol_Data.mat');
U_FOM = data.U;
fex = 100000; T=1/fex; n=5; dt=T/20; tE = n*T; tsim=5*n*T;
Vamp = 10e-9;
t = linspace(0,tsim,size(U_FOM,2));


level_test = 1;
nodes = mrdmd(U_FOM,0,0,0,6,1,false);
[Phi1,Psi1] = stitch_mrdmd(nodes,level_test);

n = 1:size(Psi1,1);
if(length(n)>1)
    figure('units','normalized','outerposition',[0 0 1 1])
    plot_styles(t,real(Psi1),'Markersize',5);
    grid;
    xlim([t(1) t(end)])
    % title(['Nonlinear model: DMD time evolution coefficients for level ' num2str(level_test)]);
    legend(cellstr(num2str(n', 'n=%-d')),'Location','best');
    change_font_sizes(1,15,40);
    change_xaxis(80,1,'$t$ in [s]');
    change_yaxis(80,1,'$\Psi$');
    change_legend(80)
    change_line_width(5)
    set(gcf,'color','w')
    export_fig(strcat('DMD_timeCoefficients_NonlinearElasticityModel_level_',num2str(level_test)),'-pdf');
end


X_mrdmd = 0*U_FOM;
for i=0:7
    [Phi,Psi] = stitch_mrdmd(nodes,i);
    X_mrdmd = X_mrdmd + Phi*Psi;
end

X_mrdmd = real(X_mrdmd);


Spacepoint = 15642;

figure('units','normalized','outerposition',[0 0 1 1])
plot(t,U_FOM(Spacepoint,:),'r');
grid;
hold on
plot(t,X_mrdmd(Spacepoint,:),'*k','MarkerSize',15);
xlim([0 tsim]);
xticks(0:2e-5:tsim);
legend('FOM','mrDMD','Location','best');
change_font_sizes(1,15,40);
change_xaxis(80,1,'$t$ in [s]');
change_yaxis(80,1,'Out-of-plane displacement in [m]');
change_legend(80)
change_line_width(5)
set(gcf,'color','w')
export_fig(strcat('FOM_vs_mrDMD_NonlinearElasticityModel_atPoint', num2str(Spacepoint)),'-pdf');

error = vecnorm(U_FOM(Spacepoint,:)-X_mrdmd(Spacepoint,:),2)./vecnorm(U_FOM(Spacepoint,:),2)*100;
disp(['The relative error in displacment mrDMD approximation at the sensor is ' num2str(error)])

