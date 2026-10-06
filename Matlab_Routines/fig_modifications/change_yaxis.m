function change_yaxis(fsize,pos,name,varargin)

if nargin==3;
    h=get(gcf,'children');
    ax=findobj(h,'Type','axes');
    ylhand = get(ax(pos),'ylabel');
    set(ylhand,'string',name,'fontsize',fsize,'Interpreter','latex');
elseif nargin==2;    
    h=get(gcf,'children');
    ax=findobj(h,'Type','axes');
    ylhand = get(ax(pos),'ylabel');
    set(ylhand,'fontsize',fsize,'Interpreter','latex');
end
    