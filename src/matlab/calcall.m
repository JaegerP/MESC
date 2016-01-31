%files = dir('allimg');
%files = files(3:length(files));
%parfor j = 3:length(files)
%    file=files(j);
%    disp(file.name)
%    i=imread(sprintf('allimg/%s',file.name));
%    dist=1;
%    if max(size(i)) == 400
%        dist=0.6375;
%    end
%    if max(size(i)) == 800
%        dist=1.275;
%    end
%    try
dist=1.275;
parfor i = [6:10]
    if i ~= 8
    try
        [cfdata,imgdata]=grainsize(sprintf('good/%i/MgO-Fe.bmp',i),dist);
        writetable(cfdata,sprintf('good/%i/MgO-Fe.bmp_cfdata.csv',i));
        writetable(imgdata,sprintf('good/%i/MgO-Fe.bmp_imgdata.csv',i));
        system(sprintf('for k in *.fig; do mv "$k" "good/%i/MgO-Fe_$k"; done',i) );
    catch 
        disp(sprintf('Error at %i',i))
    end   
    disp(sprintf('Done with MgO_Fe for %i', i));
    try
        [cfdata,imgdata]=grainsize(sprintf('good/%i/MgO-Fe-Pt.bmp',i),dist);
        writetable(cfdata,sprintf('good/%i/MgO-Fe-Pt.bmp_cfdata.csv',i));
        writetable(imgdata,sprintf('good/%i/MgO-Fe-Pt.bmp_imgdata.csv',i));
        system(sprintf('for k in *.fig; do mv "$k" "good/%i/MgO-Fe-Pt_$k"; done',i) );
    catch 
        disp(sprintf('Error at %i',i))
    end
    disp(sprintf('Done with MgO_Fe_Pt for %i', i));
    end
    
end