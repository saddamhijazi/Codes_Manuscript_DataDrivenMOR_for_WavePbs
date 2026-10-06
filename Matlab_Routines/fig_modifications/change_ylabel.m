function change_ylabel(name,pos,varargin)

    h=get(gcf,'children');
    ax=findobj(h,'Type','axes');
    ylhand = get(ax(pos),'ylabel');
    set(ylhand,'string',name,'Interpreter','latex');