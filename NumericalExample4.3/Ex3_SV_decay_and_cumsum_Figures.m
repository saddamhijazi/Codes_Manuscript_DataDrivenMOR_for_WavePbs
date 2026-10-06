clear; close all; clc;

data = load('NonlinearAluminiumModel_Comsol_Data.mat');
U_FOM = data.U;

[~,S,~] = svd(U_FOM,'econ');

figure('units','normalized','outerposition',[0 0 1 1])
semilogy(diag(S))
grid;
axis([0 500 10^-24 10^-6])
yticks([10^-20 10^-15 10^-10])
xticks([0 125 250 375 500])
change_font_sizes(1,15,40);
change_xaxis(80,1,'$r$');
change_yaxis(80,1,'$\sigma_r$');
change_line_width(5)
set(gcf,'color','w')
savefig('decay_SV_nonlinear_alu_model.fig')
export_fig('decay_SV_nonlinear_alu_model','-pdf');


figure('units','normalized','outerposition',[0 0 1 1])
plot(cumsum(diag(S))/sum(diag(S)))
grid
axis([0 500 0.1 1])
xticks([0 125 250 375 500])
yticks([0.1 0.25 0.5 0.75 1])
change_font_sizes(1,15,40);
change_xaxis(80,1,'$r$');
change_yaxis(80,1,'$c_r$');
change_line_width(5)
set(gcf,'color','w')
xlim([1 501])
savefig('CumSumSV_AluModel.fig')
export_fig('CumSumSV_AluModel','-pdf');
