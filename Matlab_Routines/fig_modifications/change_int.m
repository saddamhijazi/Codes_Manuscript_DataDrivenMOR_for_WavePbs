function change_int()

h=get(gcf,'children');
ax=findobj(h,'Type','axes');
leg=findobj(h,'Tag','legend');
for i=1:size(ax,1)
    ylhand = get(ax(i),'ylabel');
    set(ylhand,'Interpreter','latex');
        xlhand = get(ax(i),'xlabel');
    set(xlhand,'Interpreter','latex');
end
for i=1:size(leg,1)
    set(leg(i),'Interpreter','latex')
end
