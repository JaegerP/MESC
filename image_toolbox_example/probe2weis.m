path = ('/home/max/Pictures/Probe2');
liste=dir(path);
files={liste(~[liste.isdir]).name};
h=length(files);
for a = 1:h
    
    bild{a} = imread( fullfile(path,files{a}));
    
    bildzugeschnitten = imcrop(bild{a},[98.5 434.5 3300 1736]);
    bildzugeschnitten1{a} = bildzugeschnitten;
    bildgrau{a} = rgb2gray(bildzugeschnitten1{a});
end
for a = 1:(h)
    subplot(3,4,a)
       
    difference{a}=imabsdiff(bildgrau{a},bildgrau{h});
    invert{a}= imcomplement(difference{a});
    diffsw{a}= im2bw(invert{a},0.992);
    
    
    %imshow(diffsw{a});
    zwischenergebnis{a}= medfilt2(diffsw{a},[5 5]);
    
    hold on;
    imshow(bildzugeschnitten1{a});
    alpha{a} = imshow(zwischenergebnis{a});
    set(alpha{a},'AlphaData',0.15);
end
for a=1:h
figure;
imshow(bildzugeschnitten1{a});
hold on;

beta = imshow(zwischenergebnis{a});
set(beta,'AlphaData',0.35);

end




%bildsw{a} = im2bw(bildgrau{a},0.012)    %Level einstellen schwarz wei??