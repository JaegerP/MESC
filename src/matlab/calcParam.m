function [ fig, data ] = calcParam( cf, label, startParam )
%Calculate characteristic parameters
%
% Arguments:
%   *   cf:         A Correlation function, represented as an row vector
%
% Returns: a list of parameters
%   *   a figure handle containing data and fit function
%   *   struct of fit parameters [a,c,d,r]

% declare globals
global m_Distance;
%global m_HightTick;

%% Fit curves
distScale = ( 1:max( size( cf ) ) ) * m_Distance;

%m_HightTick = ( 1 );
[ fig, fr, gof ] = createFit( distScale, cf, label, startParam );

roughness=fr.r;
thickness=sqrt((fr.d)/2);%*m_HightTick;
corrlen=fr.a;
%text(0,0,sprintf('Fit parameters %s-Direction:\nroughness: %.3f\nlayer thickness: %.3f\ncorrelation length: %.3f', label, roughness, thickness, corrlen));
%disp(sprintf('Fit parameters %s-Direction:\nroughness: %.3f\ninterface layer thickness: %.3f\ncorrelation length: %.3f \nwhite noise level: %.3f * 10^-3', label, roughness, thickness, corrlen, fr.c*1000));
savefig(sprintf('CorrelationFunction_%s.fig',label));
data=[fr.a,0,fr.d,fr.r];
close all;
end

