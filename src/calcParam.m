function [ roughness, thickness, corrlen ] = calcParam( cf, label, startParam )
%Calculate characteristic parameters
%
% Arguments:
%   *   cf:         A Correlation function, represented as an row vector
%
% Returns: a list of parameters
%   *   roughness
%   *   thickness
%   *   corrlen:   Correlation Length

% declare globals
global m_Distance;
global m_HightTick;

%% Fit curves
distScale = ( 1:max( size( cf ) ) ) * m_Distance;

%m_HightTick = ( 1 );

[ fig, fitresult, gof ] = createFit( distScale, cf, label, startParam );

roughness=fitresult.r;
thickness=sqrt((fitresult.d)/2);%*m_HightTick;
corrlen=fitresult.a;
%text(0,0,sprintf('Fit parameters %s-Direction:\nroughness: %.3f\nlayer thickness: %.3f\ncorrelation length: %.3f', label, roughness, thickness, corrlen));
disp(sprintf('Fit parameters %s-Direction:\nroughness: %.3f\nlayer thickness: %.3f\ncorrelation length: %.3f', label, roughness, thickness, corrlen));
end

