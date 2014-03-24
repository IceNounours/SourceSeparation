function [ y ] = ApplyDistortion( P, x )
%APPLYDISTORTION P: polynoms to use, x: signal to distort
% P(i,j,k) : i -> output signal in which it is used,
%            j -> input signal to which it is applied
%            k -> k-th coefficient to which it is applied

    numSignals = size(x,1);
    numSamples = size(x,2);
    numCoeffs = size(P,3);
    
    y = zeros( numSignals, numSamples );

    for i=1:numSignals
        yy = zeros( 1, numSamples );
        for l=1:numSignals
            z0 = x(l,:);
            z = ones(1, numSamples );
            for k=numCoeffs:-1:1
                yy = yy .* z0 + P(i,l,k) .* z;
            end
        end
        
        y(i,:) = yy;
    end
end

