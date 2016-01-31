function [ cf_x, cf_y, avg, span, rms ] = calcCF( mat )
% calcCF - calculate correlation function in x and y direction of 
% 2D matrix mat, e.g. a greyscale image
%
% Returns: a list of parameters
%  *    cf_x: Correlation Function in x direction (rows)
%  *    cf_y: Correlation Function in y direction (columns)
%  *    avg:  the average value of all entries
%  *    span: the span of the image (maximal value - minimal value)
%  *    rms:  the Root Mean Sqare of all entries
% Invocation: calcCF( mat )
%
%==========================================================================
%
% Author:       Philipp Jaeger
% Copyright:    AG Magnetismus
%               TU Kaiserslautern
%
    span = max(max(mat)) - min(min(mat));
    avg = mean(mean(mat));
    s = size(mat);
    rms  = sqrt(1/s(1)/s(2)*sum(sum(mat)));

    %% calculate height correlation function in x direction
    cf_x=zeros(1,round(s(1)/2));
    for i = 1:round((s(1)-1)/2)	
        for k = 1:(s(1)-i) 
            cf_x(i) = cf_x(i) + sum((mat(i+k,:) - mat(k,:)).^2);
        end
        %normalize
        cf_x(i) = cf_x(i) / (s(2) * (s(1) - i));
        %disp(['Point', num2str(i) , 'of', num2str(s(1))]);       
    end
    
    %% calculate height correlation function in y direction
    cf_y=zeros(1,round(s(2)/2));
    for i = 1:round((s(2)-1)/2)	
        for k = 1:(s(2)-i) 
            cf_y(i) = cf_y(i) + sum((mat(:,i+k) - mat(:,k)).^2);
        end
        %normalize
        cf_y(i) = cf_y(i) / (s(1) * (s(2) - i));
        %disp(['Point', num2str(i) , 'of', num2str(s(1))]);
    end
end

