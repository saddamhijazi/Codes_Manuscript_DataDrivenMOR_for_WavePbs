close all;
files=dir('*.fig');
for i=1:size(files);
    filen=['./' files(i).name];
    uiopen(filen,1)
    [~,name_fil,~] =  fileparts(files(i).name);
    change_font_sizes(45,35,30)
    name_pdf = [name_fil,'.pdf'];
    export_fig(name_fil,'-pdf','-nocrop');
    cropcmd = ['pdfcrop ./',name_fil,'.pdf',' ./',name_fil,'.pdf'];
    system(cropcmd);
    close all;
end