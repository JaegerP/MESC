% determine surface properties from STM and AFM images
%
%                          _____--_____                                              
%               _____------             ------____                                                         
%              |                                  |                           
%              |                                  |                               
%              |__________________________________|                                                                                       
%                    /       /        |       |                                     
%                   /       /         |       |                                                
%                  /       /          |       |                                              
%                 /       /           |_______|                                            
%                /       /             _______                                                       
%               /       /             |       |                                     
%              /       /              |       |                                                
%             /       /               |       |                                              
%            /       /                |_______|                                                 
%           /       /                  _______                                                     
%          /       /                  |       |                                     
%         /       /                   |       |                                               
%        /       /                    |       |                                              
%       /_______/                     |_______|                                                                
%                                                             
%       TU Kaiserslautern
%       Fachbereich Physik
%       AG Magnetismus
%                                                               
%       Authors: * Philipp Jaeger
%
%================================================================
%
% This script analyses AFM/STM images and calculates the roughness
% in the lateral directions. Images have to be cropped such that
% no artifacts are visible in order to obtain reasonable data.
% Currently, only grayscale images are supported.
%
% This script is based on Philipp Fuhrmann's script written in R,   
% and on the PhD thesis of Tim Mewes, Appendix A
%
% In order to run this script, some configuration data has to be hard-coded
% in the upper code section. 
%

%% PARAMETER SETTING

% DO NOT CHANGE CODE OUTSIDE THIS SECTION, EXCEPT YOU KNOW WHAT YOU ARE 
% DOING!!!

function [cfdata,imgdata]=grainsize (m_InputPath, m_PixelDist)

% declare globals
% distance from one edge of the emage to another
global m_Distance;

% input file location
%m_InputPath='test2.bmp';

% total input image size ( before cropping !!! )
% in px
m_TotalImageSize_px = 400;

% spatial frame size of image ( pysical distance from one edge to another )
% in nm
m_TotalImageSize_nm = 100; %Probably, that's it

% xy Calibration factor (some magic factor taken from Philipp Fuhrmann's
% original script)
%m_XYCalibration = 5.1; %STM before Jan '13

% boundadies for grain area (in pixels)
% you may need to adjust this.
m_GrainAreaBoundary = [30 2000];

% boundaries for grain excentricity.
m_ExcentricityBoundary = [0 1];

% Number of bins for grainsize histogram
m_NumberOfBins = 20; % Or whatever you want.

% hight difference in which the color remains the same
%m_HightTick = total height of color scale / number of colors;
m_HightTick = 1.1191  / 256;

% estimated Layer thickness
%m_LayerInitial = 10*m_HightTick;

%% CALCULATION

%close open figures if any
close all;


m_InputImage = imread( m_InputPath );
if size(m_InputImage,3) == 3
    m_InputImage = rgb2gray( m_InputImage );
end


%m_Distance = m_TotalImageSize_nm / m_TotalImageSize_px * m_XYCalibration;
m_Distance = m_PixelDist;
%m_TotalImageSize_px = max ( size ( m_InputImage ) );

inputData = double( m_InputImage ) * m_HightTick;

[ cf_x, cf_y, avg, span, rms ] = calcCF( inputData );

% start parameters: [a c d r], function d*(1-exp(-(x/a)^r)+c)
[xfig,xdata]=calcParam( cf_x, 'x', [2 0 0 1] );%
[yfig,ydata]=calcParam( cf_y, 'y', [2 0 0 1] );
ydata=[ydata, avg, span, rms ]';
xdata=[xdata, avg, span, rms ]';

%return;

%% IMAGE ANALYSIS
%image segmentation
[binaryImage,maskedImage] = segmentImage(m_InputImage);
% Fill holes in regions.
binaryImage = imfill(binaryImage, 'holes');
% Filter image based with Image Segmenter and Image Region Analyzer.
% We are filtering on convex area to drop out connected grains and so on, 
% but finally, we need the diameter.
binaryImage = bwpropfilt(binaryImage, 'ConvexArea', m_GrainAreaBoundary);
binaryImage = bwpropfilt(binaryImage, 'Eccentricity', m_ExcentricityBoundary);
%binaryImage = bwpropfilt(binaryImage, 'EquivDiameter', [0 70]);
labeledImage = bwlabel(binaryImage);
measurements = regionprops(labeledImage, 'EquivDiameter');
%allDiameters = [measurements.EquivDiameter];
% Plot Histgram
figure('Name','Grainsize Distributuon');
[diamDistribution binDiameters] = hist(struct2array(measurements), m_NumberOfBins);
bar(binDiameters*m_Distance, diamDistribution, 'BarWidth', 1.0);
hold on;
title( sprintf('Grainsize Distributuon for %s' , m_InputPath) );
% fit function from Karoutsos, Pappaoioanou, et al., "Growth modes of 
% nanocrystalline Ni/Pt multilayers with deposition temperature", 2007
ft = fittype( 'k/s*exp(-(log(x)-log(D))^2/2/s^2)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.StartPoint=[20 1 5];
opts.Lower=[0 0 0];
opts.Upper=[100 10 100];
[xData, yData]=prepareCurveData(binDiameters*m_Distance, diamDistribution);
[fr, gof] = fit( xData, yData, ft, opts );
h=plot(fr, 'r', 'predobs');
%disp(sprintf('Fit parameters for grainsize distribution:\nMaximum: %f\nLogarithmic standard deviation: %f',fr.D,fr.s));
% make figure more fancy
set(h(1),'linewidth',2);
set(h(2),'linewidth',2);
set(h(3),'linewidth',2);
xlabel('grain diameter [nm]');
ylabel('occurences');
legend( 'measurement', 'fit function', '95% Confidence Bounds', 'Location', 'NorthEast' );
ylim([0, max(yData)+5]);
savefig('GrainsizeDist.fig');
%close all;
hold off;
%figure();
%imshow(binaryImage);
%disp(xdata);
%disp({'a','c','d','r','span','rms','avg'}');
cfdata=table({'a','c','d','r','span','rms','avg'}',xdata,ydata,'RowNames',{'a','c','d','r','span','rms','avg'}','VariableNames',{'parameter','CF_X','CF_Y'}');
imgdata=table([fr.D,fr.s]','RowNames',{'D','s'}');
end