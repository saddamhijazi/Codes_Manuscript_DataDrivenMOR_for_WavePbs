clc; clear; close all;

color25 = "#800000"; % Maroon
color27 = "#469990"; % Teal

r = 20;
data0 = load(['FOM_vs_PODGROM_Results_r' num2str(r) '.mat']);
data1 = load(['FOM_vs_DMD_Results_r' num2str(r) '.mat']);

t = data1.t;
U_ref = data1.U_s_FOM;
U_POD = data0.U_s_POD;
U_DMD = data1.U_s_DMD;
U_DMD_Opt = data1.U_s_DMD_Opt;
%%

figure('units','normalized','outerposition',[0 0 1 1])
plot(t,U_ref,'r','MarkerSize',15);
grid;
hold on
plot(t,U_POD,'--b','MarkerSize',15);
hold on
plot(t,U_DMD,'*k','MarkerSize',15);
hold on
plot(t,U_DMD_Opt,'o','MarkerSize',15,Color=color27);
hold on
legend('FOM', ['POD-ROM $r=' num2str(r) '$'] ,['DMD $r=' num2str(r) '$'], ...
    ['Optimal DMD $r=' num2str(r) '$'],'Location','best');
change_font_sizes(1,15,40);
change_xaxis(80,1,'$t$ in [s]');
change_yaxis(80,1,'Displacement in [m]');
change_legend(80)
change_line_width(5)
set(gcf,'color','w')
export_fig(['ComparisonDMD_WaveEquation_DamagedModel_r' num2str(r)],'-pdf');

pause(1);

r = 40;
data0 = load(['FOM_vs_PODGROM_Results_r' num2str(r) '.mat']);
data1 = load(['FOM_vs_DMD_Results_r' num2str(r) '.mat']);

t = data1.t;
U_POD = data0.U_s_POD;
U_DMD = data1.U_s_DMD;
U_DMD_Opt = data1.U_s_DMD_Opt;

figure('units','normalized','outerposition',[0 0 1 1])
plot(t,U_ref,'r','MarkerSize',15);
grid;
hold on
plot(t,U_POD,'--b','MarkerSize',15);
hold on
plot(t,U_DMD,'*k','MarkerSize',15);
hold on
plot(t,U_DMD_Opt,'o','MarkerSize',15,Color=color27);
hold on
legend('FOM', ['POD-ROM $r=' num2str(r) '$'] ,['DMD $r=' num2str(r) '$'], ...
    ['Optimal DMD $r=' num2str(r) '$'],'Location','best');
change_font_sizes(1,15,40);
change_xaxis(80,1,'$t$ in [s]');
change_yaxis(80,1,'Displacement in [m]');
change_legend(80)
change_line_width(5)
set(gcf,'color','w')
export_fig(['ComparisonDMD_WaveEquation_DamagedModel_r' num2str(r)],'-pdf');

%%
pause(2)

r = 20;
data0 = load(['FOM_vs_PODGROM_Results_r' num2str(r) '.mat']);
data2 = load(['FOM_vs_OpInf_Results_r' num2str(r) '.mat']);
data3 = load(['FOM_vs_dOpInf_Results_r' num2str(r) '.mat']);

U_POD = data0.U_s_POD;
U_OpInf = data2.U_s_OpInf;
U_dOpInf = data3.U_s_OpInf;
U_dOpInfReProj = data3.U_s_OpInf_ReProj;

figure('units','normalized','outerposition',[0 0 1 1])
plot(t,U_ref,'r','MarkerSize',15);
grid;
hold on
plot(t,U_POD,'--b','MarkerSize',15);
hold on
plot(t,U_OpInf,'*k','MarkerSize',15);
hold on
plot(t,U_dOpInf,'o','MarkerSize',15,Color=color27);
hold on
plot(t,U_dOpInfReProj,'+','MarkerSize',15,Color=color25);
legend('FOM', ['POD-ROM $r=' num2str(r) '$'] ,['OpInf $r=' num2str(r) '$'], ...
    ['dOpInf $r=' num2str(r) '$'], ['dOpInf ReProj $r=' num2str(r) '$'] ...
    ,'Location','best');
change_font_sizes(1,15,40);
change_xaxis(80,1,'$t$ in [s]');
change_yaxis(80,1,'Displacement in [m]');
change_legend(80)
change_line_width(5)
set(gcf,'color','w')
export_fig(['ComparisonOpInf_WaveEquation_DamagedModel_r' num2str(r)],'-pdf');


pause(1);

r = 40;
data0 = load(['FOM_vs_PODGROM_Results_r' num2str(r) '.mat']);
data2 = load(['FOM_vs_OpInf_Results_r' num2str(r) '.mat']);
data3 = load(['FOM_vs_dOpInf_Results_r' num2str(r) '.mat']);

U_POD = data0.U_s_POD;
U_OpInf = data2.U_s_OpInf;
U_dOpInf = data3.U_s_OpInf;
U_dOpInfReProj = data3.U_s_OpInf_ReProj;

figure('units','normalized','outerposition',[0 0 1 1])
plot(t,U_ref,'r','MarkerSize',15);
grid;
hold on
plot(t,U_POD,'--b','MarkerSize',15);
hold on
plot(t,U_OpInf,'*k','MarkerSize',15);
hold on
plot(t,U_dOpInf,'o','MarkerSize',15,Color=color27);
hold on
plot(t,U_dOpInfReProj,'+','MarkerSize',15,Color=color25);
legend('FOM', ['POD-ROM $r=' num2str(r) '$'] ,['OpInf $r=' num2str(r) '$'], ...
    ['dOpInf $r=' num2str(r) '$'], ['dOpInf ReProj $r=' num2str(r) '$'] ...
    ,'Location','best');
change_font_sizes(1,15,40);
change_xaxis(80,1,'$t$ in [s]');
change_yaxis(80,1,'Displacement in [m]');
change_legend(80)
change_line_width(5)
set(gcf,'color','w')
export_fig(['ComparisonOpInf_WaveEquation_DamagedModel_r' num2str(r)],'-pdf');

%%
A = load('PODROM_results_different_r.mat');
B = load('DMD_Epsilon_s_vs_r.mat');
C = load('DMD_Opt_Epsilon_s_vs_r.mat');
D = load('OpInf_results_different_r.mat');
E = load('dOpInf_results_different_r.mat');
G = load('dOpInfReproj_results_different_r.mat');


figure('units','normalized','outerposition',[0 0 1 1])
semilogy(E.r_vals,E.cond_DDT_vals,'o','MarkerSize',15,Color=color27);
grid;
hold on;
semilogy(G.r_vals,G.cond_DDT_vals,'+','MarkerSize',15,Color=color25);
xlim([0 100])
legend('dOpInf', 'dOpInf ReProj', ...
    'Location','best');
change_font_sizes(1,15,40);
change_xaxis(80,1,'$r$');
change_yaxis(80,1,'Condition number');
change_legend(80)
change_line_width(5)
set(gcf,'color','w')
export_fig('dOpInf_conditionNumber_vs_r','-pdf');


r_dmd_vals = A.r_vals;
epsilon_DMD = A.epsilon_s_vals;
epsilon_DMD_Opt = B.epsilon_s_vals;

figure('units','normalized','outerposition',[0 0 1 1])
plot(A.r_vals,A.epsilon_s_vals,'--b','MarkerSize',15);
grid;
hold on;
plot(B.r_vals,B.epsilon_s_vals,'*k','MarkerSize',15);
hold on;
plot(C.r_vals,C.epsilon_s_vals,'o','MarkerSize',15,Color=color27);
xlim([10 100])
ylim([0 40])
legend('POD-ROM' , 'DMD', ...
    'Optimal DMD', ...
    'Location','best');
change_font_sizes(1,15,40);
change_xaxis(80,1,'$r$');
change_yaxis(80,1,'$\epsilon_s$ in $\%$');
change_legend(80)
change_line_width(5)
set(gcf,'color','w')
export_fig('ComparisonDMD_WaveEquation_DamagedModel_SensorError_vs_r','-pdf');


figure('units','normalized','outerposition',[0 0 1 1])
plot(A.r_vals,A.epsilon_s_vals,'--b','MarkerSize',15);
grid;
hold on;
plot(D.r_vals,D.epsilon_s_vals,'*k','MarkerSize',15);
hold on;
plot(E.r_vals,E.epsilon_s_vals,'o','MarkerSize',15,Color=color27);
hold on;
plot(G.r_vals,G.epsilon_s_vals,'+','MarkerSize',15,Color=color25);
xlim([10 100])
legend('POD-ROM' , 'OpInf', ...
    'dOpInf', 'dOpInf ReProj', ...
    'Location','best');
change_font_sizes(1,15,40);
change_xaxis(80,1,'$r$');
change_yaxis(80,1,'$\epsilon_s$ in $\%$');
change_legend(80)
change_line_width(5)
set(gcf,'color','w')
export_fig('ComparisonOpInf_WaveEquation_DamagedModel_SensorError_vs_r','-pdf');

%%
epsilon_u_POD = data0.error_u_POD;
epsilon_u_DMD = data1.epsilon_u_DMD;
epsilon_u_DMD_Opt = data1.epsilon_u_DMDOpt;
epsilon_u_OpInf = data2.error_u_OpInf;
epsilon_u_dOpInf = data3.error_u_OpInf;
epsilon_u_dOpInfReProj = data3.error_u_OpInf_ReProj;

figure('units','normalized','outerposition',[0 0 1 1])
plot(t,epsilon_u_POD,'--b','MarkerSize',15);
grid;
hold on
plot(t,epsilon_u_DMD,'*k','MarkerSize',15);
hold on
plot(t,epsilon_u_DMD_Opt,'o','MarkerSize',15,Color=color27);
legend(['POD-ROM $r=' num2str(r) '$'] ,['DMD $r=' num2str(r) '$'] ...
    ,['Optimal DMD $r=' num2str(r) '$'],'Location','best');
change_font_sizes(1,15,40);
change_xaxis(80,1,'$t$ in [s]');
change_yaxis(80,1,'$\epsilon_u(t)$ in $\%$');
change_legend(80)
change_line_width(5)
set(gcf,'color','w')
export_fig(['Comparison_WaveEquation_DMD_DamagedModel_Epsilon_u_r' num2str(r)],'-pdf');


figure('units','normalized','outerposition',[0 0 1 1])
plot(t,epsilon_u_POD,'--b','MarkerSize',15);
grid;
hold on
plot(t,epsilon_u_OpInf,'*k','MarkerSize',15);
hold on
plot(t,epsilon_u_dOpInf,'o','MarkerSize',15,Color=color27);
hold on
plot(t,epsilon_u_dOpInfReProj,'+','MarkerSize',15,Color=color25);
legend(['POD-ROM $r=' num2str(r) '$'] ,['OpInf $r=' num2str(r) '$'], ...
    ['dOpInf $r=' num2str(r) '$'], ['dOpInf ReProj $r=' num2str(r) '$'] ...
    ,'Location','best');
change_font_sizes(1,15,40);
change_xaxis(80,1,'$t$ in [s]');
change_yaxis(80,1,'$\epsilon_u(t)$ in $\%$');
change_legend(80)
change_line_width(5)
set(gcf,'color','w')
export_fig(['Comparison_WaveEquation_OpInf_DamagedModel_Epsilon_u_r' num2str(r)],'-pdf');
