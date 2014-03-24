function [ P ] = CoinsU2PolysFast( x, numCoeffs )
%COINSU2POLYS Summary of this function goes here
%   Detailed explanation goes here

    numSamples = size( x, 2 );
    numSignals = size( x, 1 );

    dmax = max( dot(x,x) );
    di = 0.5 * dmax * ones( numSignals, numSignals, numCoeffs );
    
    % initialize the vectors
    P = zeros( numSignals, numSignals, numCoeffs );
    for i = 1:numSignals
        for k = 1:numSignals
            for l = 1:numCoeffs
                P(i,k,l) = di(i,k,l);
            end
        end
    end
    
    carre = CreateMesh( 32, 32);
    grid = TransformedGrid( ApplyDistortion(P, carre) );
    [ nW, sW ] = CompterPointsCalculerSurface( grid, x );
    
    de = di;
    maxSteps = 1000;
    step = 0;
    while( step < maxSteps )
        
        nWprev = nW;
        sWprev = sW;
        
        Delta = zeros(numSignals, numSignals, numCoeffs);
        Pt = zeros( numSignals, numSignals, numCoeffs );
        
        for i=1:numSignals
            for k=1:numSignals
                for l=1:numCoeffs
                    for m=0:1% a changer
                        Pt(i,k,l) = P(i,k,l) + 2.0 * (m-0.5) * de(i,k,l);

                            grid = TransformedGrid( ApplyDistortion(Pt, carre) );
                            [ nWt, sWt ] = CompterPointsCalculerSurface( grid, x );

                        if( Meilleure( nWprev, sWprev, nWt, sWt) )
                            %found a better config
                            nW = nWt;
                            sW = sWt;
                            P = Pt;
                            Delta(i,k,l) = 2.0 * (m-0.5) * de(i,k,l);
                        end
                    end
                end
            end
        end
        
        if( any(Delta) )
            subStep = 0;
            maxSubSteps = 10;
            
            while( (subStep < maxSubSteps) )
                for i=1:numSignals
                    for k=1:numSignals
                        for l=1:numCoeffs %continue along Delta
                            Pt(i,k,l) = P(i,k,l) + Delta(i,k,l);
                        end
                    end
                end
                
                grid = TransformedGrid( ApplyDistortion(Pt, carre) );
                [ nWt, sWt ] = CompterPointsCalculerSurface( grid, x );

                if( Meilleure( nWprev, sWprev, nWt, sWt) )
                    %found a better config
                    nW = nWt;
                    sW = sWt;
                    P = Pt;
                    Delta = 1.1 * Delta;% increase delta
                    de = 1.1 * de;
                else
                    break;
                end
                subStep = subStep + 1;
            end
        else
            de = 0.9 * de;%decrease de.
        end
        
        step = step + 1;
    end
    
    
end

function [ meilleure ] = Meilleure( nWprev, sWprev, nW, sW )
    l1 = nWprev/sWprev;
    l2 = nW/sW;
    meilleure = l2 > l1;
end

function [ nW, sW ] = CompterPointsCalculerSurface( tGrid, x )
    numSamples = size( x, 2 );
    
    nW = 0;
    sW = tGrid.surface;%surface
    
    for i=1:numSamples
        if( TransformedGridIsPointIn( tGrid, x(:,i) ) )
            nW = nW + 1;
        end
    end
end