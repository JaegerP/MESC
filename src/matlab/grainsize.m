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
% DO NOT CHANGE CODE OUTSIDE THIS SECTION, EXCEPT YOU KNOW WHAT YOU ARE 
% DOING. BE SURE TO HAVE YOUR MEASURMENT DATA BACKUPED. UNECPEXTED
% BEHAVIOUR MIGHT CAUSE LOSS OF DATA.
%

%% PARAMETER SETTING

% declare globals
% distance from one edge of the emage to another
global m_Distance;
% hight difference in which the color remains the same
global m_HightTick;
%global m_LayerInitial;

% input file location
%m_InputPath = 'testimg.j2c';
m_InputPath='1nm.png';

% total input image size ( before cropping !!! )
% in px
m_TotalImageSize_px = 1333;

% spatial frame size of image ( pysical distance from one edge to another )
% in nm
m_TotalImageSize_nm = 200;

% xy Calibration factor (some magic factor taken from Philipp Fuhrmann's
% original script)
m_XYCalibration = 5.1; %STM before Jan '13

% maximum height difference (maximum - minimum) in nm
m_TotalHeight = 10;

% boundadies for grain area (in pixels)
% you may need to adjust this.
m_GrainAreaBoundary = [50 3000];

% boundaries for grain excentricity.
m_ExcentricityBoundary = [0 1];

% Number of bins for grainsize histogram
m_NumberOfBins = 20; % Or whatever you want.

% estimated Layer thickness
%m_LayerInitial = 10*m_HightTick;

%% CALCULATION

%close open figures if any
close all;

%m_HightTick = m_TotalHeight / mean (span);
% some magic factor, too
m_HightTick = 0.61634556  / 256;

m_Distance = m_TotalImageSize_nm / m_TotalImageSize_px * m_XYCalibration;
m_InputImage = imread( m_InputPath );
if size(m_InputImage,3) == 3
    m_InputImage = rgb2gray( m_InputImage );
end
inputData = double( m_InputImage ) * m_HightTick;

[ cf_x, cf_y, avg, span, rms ] = calcCF( inputData );

% start parameters: [a d r], function d*(1-exp(-(x/a)^r))
%calcParam( cf_x, 'x', [0 0 0] );
calcParam( cf_x, 'x', [3 10 1] );
calcParam( cf_y, 'y', [3 10 1] );
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
opts.StartPoint=[20 1 1];
[xData, yData]=prepareCurveData(binDiameters*m_Distance, diamDistribution);
[fitresult, gof] = fit( xData, yData, ft, opts );
h=plot(fitresult, 'r', 'predobs');
% make figure more fancy
set(h(1),'linewidth',2);
set(h(2),'linewidth',2);
set(h(3),'linewidth',2);
xlabel('grain diameter [nm]');
ylabel('occurences');
legend( 'measurement', 'fit function', '95% Confidence Bounds', 'Location', 'NorthEast' );
ylim([0, max(yData)+5]);
savefig('GrainsizeDist.fig');
close all;