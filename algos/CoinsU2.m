function [ W ] = CoinsU2( x )
%COINSU2 Summary of this function goes here
%   Detailed explanation goes here

    numSamples = size( x, 2 );
    numSignals = size( x, 1 );
    
    c = mean(x, 2);
    di = 2.0 * c;
    
    dmax = max( dot(x,x) );
    
    % initialize the vectors
    W = zeros( numSignals, numSignals );
    for i = 1:numSignals
        W(i,i) = di(i);
    end
    
    [ nW, sW ] = CompterPointsCalculerSurface( W, x );
    
    de = c;
    maxSteps = 1000;
    step = 0;
    while( step < maxSteps )
        
        Delta = zeros(numSignals,1);
        
        for i=1:numSignals
            
            nWprev = nW;
            sWprev = sW;
            
            for l=0:1
                Wt(:,1) = W(:,1) + 2.0 * (l-0.5) * de;
                Wt(:,2) = di - Wt(:,1);

                if( (dot(Wt(:,1),Wt(:,1)) <= dmax) && (dot(Wt(:,2),Wt(:,2)) <= dmax) )
                    [ nWt, sWt ] = CompterPointsCalculerSurface( Wt, x );

                    if( Meilleure( nWprev, sWprev, nWt, sWt ) )
                        %found a better config
                        nW = nWt;
                        sW = sWt;
                        W = Wt;
                        Delta(i) = 2.0 * (l-0.5) * de(i);
                    end
                end
            end
            
        end
        
        if( any(Delta) )
            subStep = 0;
            maxSubSteps = 20;
            
            while( (subStep < maxSubSteps) )%continue along Delta
                Wt(:,1) = W(:,1) + Delta;
                Wt(:,2) = di - Wt(:,1);

                if( (dot(Wt(:,1),Wt(:,1)) <= dmax) && (dot(Wt(:,2),Wt(:,2)) <= dmax) )
                    [ nWt, sWt ] = CompterPointsCalculerSurface( Wt, x );

                    if( Meilleure( nWprev, sWprev, nWt, sWt ) )
                        %found a better config
                        nW = nWt;
                        sW = sWt;
                        W = Wt;
                        Delta = 1.1 * Delta;% increase delta
                        de = 1.1 * de;
                    else
                        break;
                    end
                end
                subStep = subStep + 1;
            end
        else
            de = 0.9 * de;%decrease de.
        end
        
        step = step + 1;
    end
    
    
end

function [ nW, sW ] = CompterPointsCalculerSurface( W, x )
    numSamples = size( x, 2 );
    
    nW = 0;
    sW = abs(det(W));%surface
    
    
    if( sW ~= 0 )
        tx = W \ x;
        uv = tx - 0.5 * ones( 2, numSamples );
        auv = abs(uv);
        
        for i = 1:numSamples
            if( auv(1,i) <= 0.5 && auv(2,i) <= 0.5 )
                nW = nW + 1;
            end
        end
    end
end

function [ meilleure ] = Meilleure( nWprev, sWprev, nW, sW )
    l1 = nWprev^3/sWprev;
    l2 = nW^3/sW;
    meilleure = l2 > l1;
    %meilleure = (nW > nWprev) || ( (nW == nWprev) && (sW < sWprev));
end
