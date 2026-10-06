clear; close all; clc;

data = load('Nonlinear_Alu_DamagedModel_SymEx_Setting1.mat');
U_FOM = data.U;
U_FOM(1:2,:) = [];
fex = 100000; T=1/fex; n=5; dt=T/20; tE = n*T; tsim=5*n*T;
Vamp = 10e-9;
t = linspace(0,tsim,size(U_FOM,2));


nodes = mrdmd(U_FOM,0,0,0,8,2,false);

X_mrdmd = 0*U_FOM;
for i=0:7
    [Phi,Psi] = stitch_mrdmd(nodes,i);
    X_mrdmd = X_mrdmd + Phi*Psi;
end

X_mrdmd = real(X_mrdmd);


Spacepoint = 10656;

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
export_fig(strcat('FOM_vs_mrDMD_NonlinearElasticityDamagedModel_atPoint', num2str(Spacepoint)),'-pdf');

error = vecnorm(U_FOM(Spacepoint,:)-X_mrdmd(Spacepoint,:),2)./vecnorm(U_FOM(Spacepoint,:),2)*100;
disp(['The relative error in displacment mrDMD approximation at the sensor is ' num2str(error)])

