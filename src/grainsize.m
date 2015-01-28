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

% input file location
m_InputPath = 'testimg.j2c';

% total input image size ( before cropping !!! )
% in px
m_TotalImageSize_px = 800;

% spatial frame size of image ( pysical distance from one edge to another )
% in nm
m_TotalImageSize_nm = 200;

% maximum height difference (maximum - minimum) in nm
m_TotalHeight = 10;

%% CALCULATION

%colse open figures if any
close all;

m_Distance = m_TotalImageSize_nm / m_TotalImageSize_px;
m_InputImage = imread( m_InputPath );

[ cf_x, cf_y, avg, span, rms ] = calcCF( m_InputImage );

m_HightTick = m_TotalHeight / double( span );

calcParam( cf_x, 'x', [10 1 3] );
calcParam( cf_y, 'y', [10 1 3] );
