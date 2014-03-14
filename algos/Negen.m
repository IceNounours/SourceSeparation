function [ W ] = Negen( z )
%NEGEN Perfomrs ICA with negentropy
%   Detailed explanation goes here

    numSignals = size(z,1);
    numSamples = size(z,2);
    
    W = zeros( numSignals );
    
    o = ones( numSamples, 1 ) / numSamples;
    
    for l=1:numSignals
        
        % Choose randomly a starting point
        w = rand(numSignals, 1);
        w = normc(w);

        prevw = zeros(numSignals, 1);

        maxSteps = 100;
        steps = 0;
        while( (abs(dot(w,prevw)) < 0.9999) && (steps < maxSteps) )
            prevw = w;
            wz = (w')*z;
            
            gwz = g( wz );
            gpwz = gp( wz );
            Egpwz = mean(gpwz);

            x = zeros( numSignals, numSamples );
            for i=1:numSignals
                x(i,:)  = z(i,:) .* gwz;
            end

            m = x * o;
            w = m - Egpwz * w;
            
            for k=1:(l-1)
                w = w - dot(w,W(:,k))*W(:,k);
            end
            
            w = normc(w);
            
            steps = steps + 1;
        end
    
        W(:,l) = w;
    end
end

function [ y ] = G( x )

    y = log( cosh( x ) );
end

function [ y ] = g( x )

    y = tanh( x );
end

function [ y ] = gp( x )
    z = tanh( x );
    y = 1.0 - (z .* z );
end
