function [ P ] = CoinsU2PolysFast( x, numCoeffs )
%COINSU2POLYS Summary of this function goes here
%   Detailed explanation goes here

    numSamples = size( x, 2 );
    numSignals = size( x, 1 );
    numVariablesToGuess = numSignals * numSignals * numCoeffs;

    dmax = max( dot(x,x) );
    di = 0.5 * dmax * ones( numSignals, numSignals * numCoeffs + 1);
    
    % initialize the vectors
    mx = mean( x, 2 );
    P = zeros(numSignals, numSignals * numCoeffs + 1);
    for i=1:numSignals
        P(i, 1 + numCoeffs * i ) = dmax;
    end
    
    numPointsPerNode = 32;
    gridSize = 32;
    pq = PointQuadTree( numPointsPerNode, x );
    carre = CreateMesh( gridSize, gridSize );
    tcarre = ApplyDistortion(P, carre);
    [ nW, sW ] = CompterPointsCalculerSurface( pq, tcarre );
    
    de = di;
    maxSteps = 1000;
    step = 0;
    while( step < maxSteps )
        
        Delta = zeros(numSignals, numSignals * numCoeffs + 1);
        Pt = zeros( numSignals, numSignals * numCoeffs + 1 );
        
        for i=1:numSignals
            for k=1:numCoeffs*numSignals
                for m=0:1% a changer
                    Pt(i,k+1) = P(i,k+1) + 2.0 * (m-0.5) * de(i,k+1);
                    Pt(:,1) = InfererCoefficientsConstants( Pt, mx );% infere the constants
                    
                    tcarre = ApplyDistortion(Pt, carre);
                    [ nWt, sWt ] = CompterPointsCalculerSurface( pq, tcarre );

                    if( Meilleure( nW, sW, nWt, sWt) )
                        %found a better config
                        nW = nWt;
                        sW = sWt;
                        P = Pt;
                        Delta(i,k+1) = 2.0 * (m-0.5) * de(i,k+1);
                    end
                end
            end
        end
        
        if( any(Delta) )
            subStep = 0;
            maxSubSteps = 10;
            
            while( (subStep < maxSubSteps) )
                for i=1:numSignals
                    for k=1:numCoeffs*numSignals %continue along Delta
                        Pt(i,k+1) = P(i,k+1) + Delta(i,k+1);
                    end
                end
                Pt(:,1) = InfererCoefficientsConstants( Pt, mx );% infere the constants
                
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

function C = InfererCoefficientsConstants( P, mx )
    co = [ 1 1/3 1 1/3; 1 1/3 1 1/3 ];
    p = P(:,2:end);
    s = sum( co .* p, 2 );
    C = mx - s;
end

function [ meilleure ] = Meilleure( nWprev, sWprev, nW, sW )
    l1 = nWprev/sWprev;
    l2 = nW/sW;
    meilleure = l2 > l1;
end

function [ nW, sW ] = CompterPointsCalculerSurface( pq, ttriangles )

    nW = PointQuadTreeCompterPoints( pq, ttriangles );
    
    % compute the surface
    sW = 0;
    numTriangles = floor( size( ttriangles, 2 ) / 3 );
    for i=1:numTriangles

        tl = ttriangles(:, (i-1)*3 + 1 : (i-1)*3 + 3);
        v1 = tl(:,1)-tl(:,2);
        v2 = tl(:,3)-tl(:,2);

        sW = sW + abs( v1(1)*v2(2) - v1(2)*v2(1) );
    end
    
end