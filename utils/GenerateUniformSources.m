function [ s ] = GenerateUniformSources( N, M, I )
%GENERATESOURCES Generate N signals (uniformly distributed) with M samples in each, given intervals
%   in I.
%   The returned matrix s is arranged like that: the columns represent
%   the samples while the rows represent the different signals.

    if( nargin == 2 )
        I = zeros(N,2);
        I(:,1) =-0.5;
        I(:,2) = 0.5;
    end
    
    s = rand( N, M );
    
    for i=1:N
        s(i,:) = I(i,1) + (I(i,2) - I(i,1)) .* s(i,:);
    end
end

