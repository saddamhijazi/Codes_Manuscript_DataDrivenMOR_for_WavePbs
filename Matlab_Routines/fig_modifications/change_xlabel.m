function change_xlabel(name,pos,varargin)

    h=get(gcf,'children');
    ax=findobj(h,'Type','axes');
    ylhand = get(ax(pos),'xlabel');
    set(ylhand,'string',name,'Interpreter','latex');
   