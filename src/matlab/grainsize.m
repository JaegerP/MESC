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
global m_Distance;
global m_HightTick;
%global m_LayerInitial;

% input file location
m_InputPath = 'testimg.j2c';
%m_InputPath='../no7.bmp';

% total input image size ( before cropping !!! )
% in px
m_TotalImageSize_px = 800;

% spatial frame size of image ( pysical distance from one edge to another )
% in nm
m_TotalImageSize_nm = 200;

%xy Calibration factor 
xyCalibration = 5.1; %STM before Jan '13

% maximum height difference (maximum - minimum) in nm
m_TotalHeight = 10;

% estimated Layer thickness
%m_LayerInitial = 10*m_HightTick;

%% CALCULATION

%colse open figures if any
close all;

%m_HightTick = m_TotalHeight / mean (span);
m_HightTick = 0.61634556  / 256;

m_Distance = m_TotalImageSize_nm / m_TotalImageSize_px * xyCalibration;
m_InputImage = imread( m_InputPath );
if size(m_InputImage,3) == 3
    m_InputImage = rgb2gray( m_InputImage );
end
inputData = double( m_InputImage ) * m_HightTick;

[ cf_x, cf_y, avg, span, rms ] = calcCF( inputData );

%calcParam( cf_x, 'x', [0 0 0] );
calcParam( cf_x, 'x', [5 2 1] );


labeledImage = bwlabel(binaryImage);
measurements = regionprops(labeledImage, 'EquivDiameter');
allDiameters = [measurements.EquivDiameter];
numberOfBins = 50; % Or whatever you want.
[diamDistribution binDiameters] = hist(allDiameters, numberOfBins);
bar(binDiameters, diamDistribution, 'BarWidth', 1.0);