function change_linestyle(nline,linestyle,marker,color)

hline = findobj(gcf, 'type', 'line');
set(hline(nline),'LineStyle',linestyle)
set(hline(nline),'Marker',marker) 
set(hline(nline),'color',color) 