close all
files=dir('*.fig');
for i=1:size(files);
filen=['./' files(i).name];
uiopen(filen,1)
change_title(10,' ')
savefig(filen)
end
close all