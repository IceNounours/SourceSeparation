function [ W ] = KurtSymmetricOrtho( z )
%KurtSymmetricOrtho find a first directing vector for whitened data z.
%   Detailed explanation goes here

    numSignals = size(z,1);
    numSamples = size(z,2);
    
    W = rand( numSignals );
    W = normc( W );
    
    prevW = zeros( numSignals );
    o = ones( numSamples, 1 ) / numSamples;
    
    maxSteps = 100;
    steps = 0;
    while( (min( abs(dot(W,prevW)) ) < 0.9999 )&& (steps < maxSteps) )
        prevW = W;

        for j = 1 : numSignals
            w = W(:,j);
            
            wz = (w')*z;
            s3 = wz .* wz .* wz;

            x = zeros( numSignals, numSamples );
            for i=1:numSignals
                x(i,:)  = z(i,:) .* s3;
            end

            m = x * o;
            w = m - 3*w;
            
            W(:,j) = w;
        end

        % symmetric orthogonalization
        W = (sqrtm(W*W')) \ W;
        
        steps = steps + 1;
    end
end

