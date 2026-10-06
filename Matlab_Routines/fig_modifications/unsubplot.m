function []=unsubplot(n,m,k,f);
% IO:
% unsubplot(n,m,k,f);
% n,m,k corresponding to subplot(n,m,k)
% f corresponing to figure(f)
% Please note don't work with subplots where
% suptitle have been used
%
% Anders Björk 2000-08-30
g=findobj;
%max(g((g-ceil(g)==0)));
% if ~any(f==g);
%    error('Non existing figure number');
% end;

h=figure(f);
hh=subplot(n,m,k);
c=copyobj(hh,get(hh,'Parent'));
f=figure;
set(c,'Parent',f);
ha=gca;
set(c,'Position','default');
hay=get(ha,'Ylabel');
hax=get(ha,'Xlabel');
NewFontSize=10;
set(hay,'Fontsize',NewFontSize);
set(hax,'Fontsize',NewFontSize);

