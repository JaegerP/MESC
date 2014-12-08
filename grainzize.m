%
% grainsizeanalyser
% TU Kaiserslautern, AG Magnetism
%
% authors: Philipp Jaeger
%
% input file location
m_InputPath = 'sample2.jpg';
%
grainMinSize = 5;
grainMaxSize = 200;


%m_InputImage = imcrop(imread(m_InputPath),[138.5 8.5 718 718]);
m_InputImage = imread(m_InputPath);
%m_ColorScale = 
%m_ImageGrey = rgb2gray(m_InputImage);
m_ImageHistogram = imhist(m_InputImage);
%imshow(m_InputImage);
%figure; xlim([80 130]); plot(m_ImageHistogram);
m_ImageHistogram=m_ImageHistogram/max(m_ImageHistogram);
figure; bar(1:256,m_ImageHistogram); xlim([80,120]);
m_ImageBW = im2bw(m_InputImage,104.0/255);
%imshow(m_ImageBW);
BW2 = bwpropfilt(m_ImageBW, 'ConvexArea',[grainMinSize grainMaxSize]);
figure
%subimage(BW2)
subplot(1,3,1),imshow(m_InputImage)
subplot(1,3,2),imshow(m_ImageBW)
subplot(1,3,3),imshow(BW2)

