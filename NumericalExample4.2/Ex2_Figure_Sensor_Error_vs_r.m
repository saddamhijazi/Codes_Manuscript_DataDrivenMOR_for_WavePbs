clc; clear; close all;
load('OpInf_results_different_r.mat')

figure('units','normalized','outerposition',[0 0 1 1])
plot(r_vals,epsilon_s_vals,'r');
hold on
rectangle('Position',[60 40 35 38],...
          'EdgeColor','b','LineWidth',3)
axes('Position',[0.6 0.55 0.25 0.3]) % [left bottom width height]
box on
plot(r_vals(6:10), epsilon_s_vals(6:10), 'r','LineWidth',1.5)
change_font_sizes(1,15,40);
change_xaxis(60,1,'$r$');
change_yaxis(60,1,'$\epsilon_s$ in $\%$');
change_title(40,'zoomed picture');
change_xaxis(80,2,'$r$');
change_yaxis(80,2,'$\epsilon_s$ in $\%$');
change_legend(60)
change_line_width(5)
set(gcf,'color','w')
export_fig('OpInf_DamagedFML_SensorError_vs_r','-pdf');