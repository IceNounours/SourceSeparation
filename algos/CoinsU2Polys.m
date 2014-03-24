function [ P ] = CoinsU2Polys( x, numCoeffs )
%COINSU2POLYS Summary of this function goes here
%   Detailed explanation goes here

    numSamples = size( x, 2 );
    numSignals = size( x, 1 );
    
    carre = GenerateUniformSources( numSignals, numSamples );

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
    
    kP = DistanceCompacte( ApplyDistortion(P, carre), x);
    
    de = di;
    maxSteps = 1000;
    step = 0;
    while( step < maxSteps )
        
        kPprev = kP;
        Delta = zeros(numSignals, numSignals, numCoeffs);
        Pt = zeros( numSignals, numSignals, numCoeffs );
        
        for i=1:numSignals
            for k=1:numSignals
                for l=1:numCoeffs
                    for m=0:1% a changer
                        Pt(i,k,l) = P(i,k,l) + 2.0 * (m-0.5) * de(i,k,l);

                        kPt = DistanceCompacte( ApplyDistortion(Pt, carre), x);

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
                
                kPt = DistanceCompacte( ApplyDistortion(Pt, carre), x);

                if( Meilleure( kPprev, kPt ))
                    %found a better config
                    kP = kPt;
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

function [ meilleure ] = Meilleure( kprev, k )
    meilleure = k < kprev;
end

function [k] = DistanceCompacte ( y, x )
    numSignals = size( x, 1 );
    numSamples = size( x, 2 );
    
%     yy = y;
%     k = 0;
%     ctrl = zeros (1, numSamples);
%     for i = 1:numSamples
%         dm = +Inf;
%         ind = 1;
%         
%         xx = ones( numSignals, numSamples );
%         for j=1:numSignals
%             xx(j,:) = x(j,i) .* xx(j,:);
%         end
%         zz = xx - yy;
%         dists = dot( zz, zz );%distance between x(:,i) and all the y
%         
%         for j = 1:numSamples
%             m = dists( j );
%             if ( (m < dm) && (ctrl(j) == 0) )
%                 dm = m;
%                 ind = j;
%             end
%         end
%         
%         ctrl(ind) = 1;
%         if (dm>k)
%             k=dm;
%         end
%     end

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