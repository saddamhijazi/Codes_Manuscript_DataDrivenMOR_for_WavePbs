figfiles = dir('*.fig');

for i=1:size(figfiles,1)
    [~,name_fil,~] =  fileparts(figfiles(i).name);
    name = ['./',figfiles(i).name];
    export_fig(name,'-pdf');
    % export_fig(name,'-pdf','-nocrop');
    % cropcmd = ['pdfcrop ./',name_fil,'.pdf',' ./',name_fil,'.pdf'];
    % system(cropcmd);
end
