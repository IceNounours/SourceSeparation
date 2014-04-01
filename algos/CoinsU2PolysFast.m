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
    
    numPointsPerNode = 32;
    gridSize = 32;
    pq = PointQuadTree( numPointsPerNode, x );
    carre = CreateMesh( gridSize, gridSize );
    tcarre = ApplyDistortion(P, carre);
    [ nW, sW ] = CompterPointsCalculerSurface( pq, tcarre );
    
    de = di;
    maxSteps = 10;
    step = 0;
    while( step < maxSteps )
        
        Delta = zeros(numSignals, numSignals, numCoeffs);
        Pt = zeros( numSignals, numSignals, numCoeffs );
        
        for i=1:numSignals
            for k=1:numSignals
                for l=1:numCoeffs
                    for m=0:1% a changer
                        Pt(i,k,l) = P(i,k,l) + 2.0 * (m-0.5) * de(i,k,l);

                        tcarre = ApplyDistortion(Pt, carre);
                        [ nWt, sWt ] = CompterPointsCalculerSurface( pq, tcarre );

                        if( Meilleure( nW, sW, nWt, sWt) )
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
                
                tcarre = ApplyDistortion(Pt, carre);
                [ nWt, sWt ] = CompterPointsCalculerSurface( pq, tcarre );

                if( Meilleure( nW, sW, nWt, sWt) )
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

function [ nW, sW ] = CompterPointsCalculerSurface( pq, ttriangles )

    nW = PointQuadTreeCompterPoints( pq, ttriangles );
    sW = 0;
    
    numTriangles = floor( size( ttriangles, 2 ) / 3 );

    for i=1:numTriangles

        tl = ttriangles(:, (i-1)*3 + 1 : (i-1)*3 + 3);
        v1 = tl(:,1)-tl(:,2);
        v2 = tl(:,3)-tl(:,2);

        sW = sW + abs( v1(1)*v2(2) - v1(2)*v2(1) );
    end
    
end