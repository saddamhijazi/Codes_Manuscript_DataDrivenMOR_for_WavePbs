function change_line_width(lwidth)

h=get(gcf,'children');
lin=findobj(h,'Type','line');

for i=1:size(lin,1)
   set(lin(i),'LineWidth',lwidth);
end
