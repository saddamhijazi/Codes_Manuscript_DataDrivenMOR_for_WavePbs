function change_font_sizes(fontaxes,fontlegend,fontgca,filename,varargin)

if nargin==4;
    close all
    uiopen(filename,1);
    h=get(gcf,'children');
    ax=findobj(h,'Type','axes');
    leg=findobj(h,'Tag','legend');
    for i=1:size(ax,1)
        set(ax(i),'FontSize',fontgca)
        ylhand = get(ax(i),'ylabel');
        %set(ylhand,'Interpreter','latex','FontSize',fontaxes);
        set(ylhand,'FontSize',fontaxes);
        xlhand = get(ax(i),'xlabel');
        %set(xlhand,'Interpreter','latex','FontSize',fontaxes);
        set(xlhand,'FontSize',fontaxes);
        set(ax(i),'Xgrid','on')
        set(ax(i),'Ygrid','on')
    end
    for i=1:size(leg,1)
        set(leg(i),'Interpreter','latex')
        set(leg(i),'FontSize',fontlegend)
    end
    export_fig(filename);
    close all
    
    
else
    h=get(gcf,'children');
    ax=findobj(h,'Type','axes');
    leg=findobj(h,'Tag','legend');
    for i=1:size(ax,1)
        set(ax(i),'FontSize',fontgca)
        ylhand = get(ax(i),'ylabel');
        set(ylhand,'Interpreter','latex','FontSize',fontaxes);
        set(ylhand,'FontSize',fontaxes);
        xlhand = get(ax(i),'xlabel');
        set(xlhand,'Interpreter','latex','FontSize',fontaxes);
        set(xlhand,'FontSize',fontaxes);
        set(ax(i),'Xgrid','on')
        set(ax(i),'Ygrid','on')
        xdata = get(findobj(ax(1), 'type', 'line'), 'XData');
        %         set(ax(i),'XLim',[min(min(xdata{:})) max(max(xdata{:}))])
        set(ax(i),'Box','on')
    end
    for i=1:size(leg,1)
        set(leg(i),'Interpreter','latex')
        set(leg(i),'FontSize',fontlegend)
    end
end





