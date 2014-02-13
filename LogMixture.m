function [ y ] = LogMixture( s )
%GENERATEMIXTURE Summary of this function goes here
%   Detailed explanation goes here

    y(1,:) = s(1,:);
    y(2,:) = log( 1 + log( (1 + s(1,:)) .* (1 + s(2,:)) ) );
end