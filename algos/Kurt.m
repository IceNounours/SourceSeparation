function [ W ] = Kurt(z)
%KURTOSIS find a first directing vector for whitened data z.
%   Detailed explanation goes here

    numSignals = size(z,1);
    numSamples = size(z,2);
    
    W = zeros( numSignals );
    
    for l=1:numSignals
        
        % Choose randomly a starting point
        w = rand(numSignals, 1);
        w = normc(w);

        o = ones( numSamples, 1 ) / numSamples;

        prevw = zeros(numSignals, 1);

        maxSteps = 100;
        steps = 0;
        while( (abs(dot(w,prevw)) < 0.9999) && (steps < maxSteps) )
            prevw = w;
            wz = (w')*z;
            s3 = wz .* wz .* wz;

            x = zeros( numSignals, numSamples );
            for i=1:numSignals
                x(i,:)  = z(i,:) .* s3;
            end

            m = x * o;
            w = m - 3*w;
            
            for k=1:(l-1)
                w = w - dot(w,W(:,k))*W(:,k);
            end
            
            w = normc(w);
            
            steps = steps + 1;
        end
    
        W(:,l) = w;
    end
end

