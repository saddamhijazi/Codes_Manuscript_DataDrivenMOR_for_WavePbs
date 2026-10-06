function change_legendlabel(name,pos1,pos2,varargin)

    h=get(gcf,'children');
    ax=findobj(h,'Type','axes');
    ylhand = get(ax(pos1),'legend');
    ylhand.String{pos2} = name;
   