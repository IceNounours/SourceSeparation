function [ s ] = GenerateGaussianSources( N, M, U, S )
%GENERATESOURCES Generate N signals (uniformly distributed) with M samples in each, given
%   means U and variances S.
%   The returned matrix s is arranged like that: the columns represent
%   the samples while the rows represent the different signals.

    s = randn( N, M );

    for i=1:N
        s(i,:) = U(i) + S(i) .* sign( s(i,:) ) .* s(i,:).^2;
    end
end

