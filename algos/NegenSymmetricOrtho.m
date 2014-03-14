function [ W ] = NegenSymmetricOrtho( z )
%NEGEN Perfomrs ICA with negentropy
%   Detailed explanation goes here

    numSignals = size(z,1);
    numSamples = size(z,2);
    
    W = rand( numSignals );
    W = normc( W );
    
    prevW = zeros( numSignals );
    
    o = ones( numSamples, 1 ) / numSamples;
    
    maxSteps = 100;
    steps = 0;
     while( (min( abs(dot(W,prevW)) ) < 0.9999) && (steps < maxSteps) )
        prevW = W;
        for l=1:numSignals
            w = W(:,l);
            wz = (w')*z;

            gwz = g( wz );
            gpwz = gp( wz );
            Egpwz = mean(gpwz);

            x = zeros( numSignals, numSamples );
            for i=1:numSignals
                x(i,:) = z(i,:) .* gwz;
            end

            m = x * o;
            w = m - Egpwz * w;

            W(:,l) = w;
        end
        
        % symmetric orthogonalization
        W = (sqrtm(W*W')) \ W;
        
        steps = steps + 1;
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
