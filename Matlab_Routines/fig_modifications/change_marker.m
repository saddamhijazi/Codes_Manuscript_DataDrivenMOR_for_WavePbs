function change_marker(nline,marker,markersize)

hline = findobj(gcf, 'type', 'line');
set(hline(nline),'Marker',marker) 
set(hline(nline),'MarkerSize',markersize) 