function [ y ] = ApplyDistortionR( P, M, x )
%APPLYDISTORTION P: polynoms to use, x: signal to distort
% P(i,j,k) : i -> output signal in which it is used,
%            j -> input signal to which it is applied
%            k -> k-th coefficient to which it is applied
% M(i) : mean to add to the signals

    numSignals = size(x,1);
    numSamples = size(x,2);
    numCoeffs = size(P,3);
    
    y = zeros( numSignals, numSamples );

    for i=1:numSignals
        y(i,:) = zeros( 1, numSamples );
        for l=1:numSignals
            yy = zeros( 1, numSamples );
            z0 = x(l,:);
            for k=numCoeffs:-1:1
                yy = (yy + P(i,l,k)) .* z0;
            end
            y(i,:) = y(i, :) + yy;
        end
        
        y(i,:) = y(i, :) + M(i);
    end
end

