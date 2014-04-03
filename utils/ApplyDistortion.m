function [ y ] = ApplyDistortion( P, x )
%APPLYDISTORTION P: polynoms to use, x: signal to distort

    numSignals = size(x,1);
    numSamples = size(x,2);
    numCoeffs = ( size(P,2) - 1 ) / numSignals;
    
    y = zeros( numSignals, numSamples );

    for i=1:numSignals
        for l=1:numSignals
            yy = zeros( 1, numSamples );
            z0 = x(l,:);
            for k=numCoeffs:-1:1
                yy = (yy + P( i, 1 + numCoeffs * (l-1) + k )) .* z0;
            end
            
            y(i,:) = y(i,:) + yy;
        end
        
        y(i,:) = y(i,:) + P( i, 1 ) * ones( 1, numSamples );
    end
end

