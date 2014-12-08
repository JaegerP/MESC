path = ('/home/max/Pictures/Probe1(2)');
liste=dir(path);
files={liste(~[liste.isdir]).name};
h=length(files);
for a = 1:h
    
    bild{a} = imread( fullfile(path,files{a}));
    
    bildzugeschnitten = imcrop(bild{a},[18.5 270.5 2284 2200]);
    bildzugeschnitten1{a} = bildzugeschnitten;
    bildgrau{a} = rgb2gray(bildzugeschnitten1{a});
end
for a = 1:(h)
    
       
    difference{a}=imabsdiff(bildgrau{a},bildgrau{h});
    invert{a}= imcomplement(difference{a});
    diffsw{a}= im2bw(invert{a},0.995);
    
    
    %imshow(diffsw{a});
    zwischenergebnis{a}= medfilt2(diffsw{a},[3 3]);
    
 %   hold on;
  %  imshow(bildzugeschnitten1{a});
    
   % alpha{a} = imshow(zwischenergebnis{a});
    %set(alpha{a},'AlphaData',0.20);
end
for a=1:h
figure;
imshow(bildzugeschnitten1{a});
hold on;
pixel = ones(1,2201)*zwischenergebnis{a}*ones(2285,1);
pixelanfang = ones(1,2201)*zwischenergebnis{1}*ones(2285,1);
pixelinprozent= pixel/pixelanfang;
disp(pixelinprozent);
beta = imshow(zwischenergebnis{a});
set(beta,'AlphaData',0.35);

end




%bildsw{a} = im2bw(bildgrau{a},0.012)    %Level einstellen schwarz wei??