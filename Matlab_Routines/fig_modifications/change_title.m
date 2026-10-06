function change_title(fsize,name,varargin)

if nargin==2;
    tihand = get(gca,'title');
    set(tihand,'string',name,'fontsize',fsize,'Interpreter','latex');
else
    tihand = get(gca,'title');
    set(tihand,'fontsize',fsize);
end
    