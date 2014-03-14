function [ s ] = GenerateBernouilliSources( numSignals, numSamples, des )
%GENERATEBERNOUILLISOURCES Generate bernouilli sources

    X = rand( numSignals, numSamples );
    s = zeros( numSignals, numSamples );
    
    for i = 1 : numSignals
        p = des(i);
        for j = 1 : numSamples
            if( X(i,j) < p )
                s(i,j) = 0;
            else
                s(i,j) = 1;
            end
        end
    end
end

