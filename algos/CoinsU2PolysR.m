function [ P, M ] = CoinsU2PolysR( x, numCoeffs )
%COINSU2POLYS Summary of this function goes here
%   Detailed explanation goes here

    numSamples = size( x, 2 );
    numSignals = size( x, 1 );
    
    carre = GenerateUniformSources( numSignals, numSamples );

    dmax = max( dot(x,x) );
    de = 0.5 * dmax * ones( numSignals, numSignals, numCoeffs );
    dm = 0.5 * dmax * ones( 1, numSignals );
    
    % initialize the vectors
    M = mean( x, 2 );
    P = zeros( numSignals, numSignals, numCoeffs );
    for i = 1:numSignals
        for k = 1:numSignals
            for l = 1:numCoeffs
                P(i,k,l) = de(i,k,l);
            end
        end
    end
    
    kP = DistanceCompacte( ApplyDistortionR(P, M, carre), x);
    
    maxSteps = 100;
    step = 0;
    while( step < maxSteps )
        
        kPprev = kP;
        Delta = zeros(numSignals, numSignals, numCoeffs);
        DeltaM = zeros(1, numSignals);
        
        Pt = zeros( numSignals, numSignals, numCoeffs );
        Mt = zeros( 1, numSignals );
        
        for i=1:numSignals
            for k=1:numSignals
                for l=1:numCoeffs
                    for m=0:1% a changer
                        Pt(i,k,l) = P(i,k,l) + 2.0 * (m-0.5) * de(i,k,l);

                        kPt = DistanceCompacte( ApplyDistortionR(Pt, Mt, carre), x);

                        if( Meilleure( kPprev, kPt ) )
                            %found a better config
                            kP = kPt;
                            P = Pt;
                            Delta(i,k,l) = 2.0 * (m-0.5) * de(i,k,l);
                        end
                    end
                end
            end
        end
        
        for i=1:numSignals
            for m=0:1% a changer
                Mt(i) = Mt(i) + 2.0 * (m-0.5) * de(i,k,l);

                kPt = DistanceCompacte( ApplyDistortionR(Pt, Mt, carre), x);

                if( Meilleure( kPprev, kPt ) )
                    %found a better config
                    kP = kPt;
                    P = Pt;
                    DeltaM(i) = 2.0 * (m-0.5) * dm(i);
                end
            end
        end
        
        aD = any( any( any(Delta) ) );
        aDm = any(DeltaM);
        if( aD || aDm )
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
                
                for i=1:numSignals
                    for m=0:1% a changer
                        Mt(i) = M(i) + DeltaM(i);
                    end
                end
                
                kPt = DistanceCompacte( ApplyDistortionR(Pt, Mt, carre), x);

                if( Meilleure( kPprev, kPt ))
                    %found a better config
                    kP = kPt;
                    P = Pt;
                    Delta = 1.1 * Delta;% increase delta
                    DeltaM = 1.1 * DeltaM;
                    de = 1.1 * de;
                    dm = 1.1 * dm;
                else
                    break;
                end
                subStep = subStep + 1;
            end
        else
            de = 0.9 * de;%decrease de.
            dm = 0.9 * dm;
        end
        
        step = step + 1;
    end
    
    
end

function [ meilleure ] = Meilleure( kprev, k )
    meilleure = k < kprev;
end

function [k] = DistanceCompacte ( y, x )
    numSignals = size( x, 1 );
    numSamples = size( x, 2 );

    yy = y;
    k = 0;
    for i = 1:numSamples
        t = (i-1);
        
        dm = +Inf;
        ind = 1;
        
        xx = ones( numSignals, numSamples -t );
        for l=1:numSignals
            xx(l,:) = x(l,i) .* xx(l,:);
        end
        zz = xx - yy(:,i:numSamples);
        dists = dot( zz, zz );%distance between x(:,i) and all the y
        
        for l = 1:(numSamples -t)
            m = dists( l );
            if ( (m < dm) )
                dm = m;
                ind = l;
            end
        end
        
        % swap to the beginning
        temp = yy(:,ind+t);
        yy(:,ind+t) = yy(:,i);
        yy(:,i) = temp;
        
        if (dm>k)
            k=dm;
        end
    end
end