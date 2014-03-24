function [ W ] = Coins3( x, p, q, r, s )
%COINSLIKELIHOOD Summary of this function goes here
%   Detailed explanation goes here
    numSamples = size( x, 2 );
    numSignals = size( x, 1 );
    
    
    c = mean( x, 2 );% center of the rect
    
    % centered signal
    z = zeros( numSignals, numSamples );
    z(1,:) = x(1,:) - c(1,:);
    z(2,:) = x(2,:) - c(2,:);
    z(3,:) = x(3,:) - c(3,:);
    
    m = 2.0 * c;
    B0 = c;% (A1+A2+A3)/2
    n0 = normc(c);%normalized vector in one of the separating planes
    
    B = zeros( numSignals, 4);
    B(:,1) = B0;
    N= zeros( numSignals, 4 );
    N(:,1)= n0;
    
    [ U, O, V ] = SeparatePoints( z, c, 6.0 * s );
    V0 = V;
    for i=1:3
        
        figure;
        PlotSignal( V );
        title('V avant création du plan');
        
        n = FindPlaneNormal( V, V0, n0 );
        N(:,i+1) = n;
        
        % vector orienting the count
        d = cross( n, n0 );
        d = normc(d);
        
        [ Ip, Op ] = SeparatePointsFromPlane( V, n, 6.0*s );
        u = FindDiagonal( Ip, d );
        B(:,i+1) = u;
        
        figure;
        PlotSignal( U, 'b' );
        hold on;
        PlotSignal( O, 'r' );
        hold on;
        PlotSignal( Ip, 'y' );
        hold on;
        PlotSignal( Op, 'g' );
        hold on;
        PlotVectors( n, 'r' );
        hold on;
        PlotVectors( u, 'b' );
        title('Etape');
        
        
        V = Op;
    end
    
    d = Inf;
    jopt = 8;
    Copt = zeros( numSignals );
    
    for l=0:7
        
        S = zeros( numSignals, 1 );
        C = zeros(numSignals);
        
        if( bitand(l, 1) ~= 0 )
            S = S + B(:,2);
            C(:,1) = B(:,2);
        else
            S = S - B(:,2);
            C(:,1) = - B(:,2);
        end
        
        if( bitand(l, 2) ~= 0 )
            S = S + B(:,3);
            C(:,2) = B(:,3);
        else
            S = S - B(:,3);
            C(:,2) = - B(:,3);
        end
        
        if( bitand(l, 4) ~= 0 )
            S = S + B(:,4);
            C(:,3) = B(:,4);
        else
            S = S - B(:,4);
            C(:,3) = - B(:,4);
        end
        
        d0 = dot( S-B(:,1), S-B(:,1) );
        
        if( d0 < d )
            d = d0;
            jopt = l;
            Copt = C;
        end
    end
    
    
    A(:,1) = Copt(:,1) + Copt(:,2);
    A(:,2) = Copt(:,2) + Copt(:,3);
    A(:,3) = Copt(:,3) + Copt(:,1);
    
    
    % Convergence step
    scale = 4.0;
    steps = 0;
    maxSteps = 100;
    
    e = m - sum(A,2);
    
    Br = zeros(numSignals, numSamples );
    Ul = zeros(numSignals, numSamples );
    Ap = zeros(numSignals, numSamples );
    
    while( (steps < maxSteps) && (norm(e) > norm(m)*0.0001) )
        
        prevE = e;
        
        numPointsBr = size(Br,2);
        numPointsUl = size(Ul,2);
        numPointsAp = size(Ap,2);

        nUl = 0;
        nBr = 0;
        nAp = 0;
    
        if( numPointsAp > numSamples / 100 )
            for j=1:numSamples
                xx = x(1,j);
                xy = x(2,j);
                xz = x(3,j);

                if ( ((xx-A(1,3))/(scale*s(1)))^2 + ((xy-A(2,3))/(scale*s(2)))^2 + ((xz-A(3,3))/(scale*s(3)))^2 < 1 )
                    Ap(:,nAp+1) = x(:,j);
                    nAp = nAp + 1;
                end
            end
            
            Ap = Ap(:,1:nAp);
        end
        
        if( numPointsBr > numSamples / 100 )
            for j=1:numSamples
                xx = x(1,j);
                xy = x(2,j);
                xz = x(3,j);

                if ( ((xx-A(1,2))/(scale*s(1)))^2 + ((xy-A(2,2))/(scale*s(2)))^2 + ((xz-A(3,2))/(scale*s(3)))^2 < 1 )
                    Br(:,nBr+1) = x(:,j);
                    nBr = nBr + 1;
                end
            end
            
            Br = Br(:,1:nBr);
        end
              
        if( numPointsUl > numSamples / 100 )
            for j=1:numSamples
                xx = x(1,j);
                xy = x(2,j);
                xz = x(3,j);

                if ( ((xx-A(1,1))/(scale*s(1)))^2 + ((xy-A(2,1))/(scale*s(2)))^2 +((xz-A(3,1))/(scale*s(3)))^2 < 1 )
                    Ul(:,nUl+1) = x(:,j);
                    nUl = nUl + 1;
                end
            end
            
            Ul = Ul(:,1:nUl);
        end
        
        % compute the vectors
        A(:,3) = mean( Ap, 2 );
        A(:,2) = mean( Br, 2 );
        A(:,1) = mean( Ul, 2 );

        steps = steps + 1;
        e = m - sum(A,2);
        
        %if( dot(e,e) < dot(prevE,prevE) )
            scale = scale * 0.9;
        %end
    end
    
    figure;
    PlotSignal( z, 'b' );
    hold on;
    PlotVectors([ N(:,2) N(:,3) N(:,4) ], 'r' );
    hold on;
    PlotVectors( A, 'b' );
    title('regions');
    
    W = A;
end

function [U, O, V ] = SeparatePoints( z, u, s )
% Separate the points along the diagonal u
    numSamples = size( z, 2 );
    numSignals = size( z, 1 );
    
    % points near u
    U = zeros( numSignals, floor(numSamples/8));
    nU = 0;
    % points near the -u
    O = zeros( numSignals, floor(numSamples/8));
    nO = 0;
    % other points
    V = zeros( numSignals, floor(6*numSamples/8) );
    nV = 0;
    
    for i=1:numSamples
        
        xx = z(1,i);
        yy = z(2,i);
        zz = z(3,i);
        
        if( abs((xx + u(1))/s(1))^2 + abs((yy + u(2))/s(2))^2 + abs((zz+u(3))/s(3))^2 < 1 )
            nO = nO + 1;
            O(:,nO) = z(:,i);
        elseif( abs((xx - u(1))/s(1))^2 + abs((yy - u(2))/s(2))^2 + abs((zz-u(3))/s(3))^2 < 1 )
            nU = nU + 1;
            U(:,nU) = z(:,i);
        else
            nV = nV + 1;
            V(:,nV) = z(:,i);
        end
    end
    
    U = U(:,1:nU);
    O = O(:,1:nO);
    V = V(:,1:nV);
end

function [I, O ] = SeparatePointsFromPlane( z, n, s )
% Separate the points, with n the normal of the plane
    numSamples = size( z, 2 );
    numSignals = size( z, 1 );
    
    % points near the plane
    I = zeros( numSignals, floor(numSamples/8) );
    nI = 0;
    % points outside the plane
    O = zeros( numSignals, floor(numSamples/8) );
    nO = 0;
    
    nz = normc(z);
    
    for i=1:numSamples
        
        % we should take infos from the variance but, would
        % need to compute the density of the square projected on sphere...
        if( abs(dot( nz(:,i), n )) < 0.1 )
            nI = nI + 1;
            I(:, nI) = z(:,i);
        else
            nO = nO + 1;
            O(:, nO) = z(:,i);
        end
    end
    
    I = I(:,1:nI);
    O = O(:,1:nO);
end

function [ uu ] = FindPlaneNormal( z, z0, n )
% Find the normal of a plane passing by the center of the cube and
% containing the vector n
    numSignals = size( z, 1 );
    nV = size( z0, 2 );
    
    cV = z0;
    cV = normc(cV);
    
    N = ones( numSignals, nV);
    N(1,:) = N(1,:) .* n(1);
    N(2,:) = N(2,:) .* n(2);
    N(3,:) = N(3,:) .* n(3);
    we = cross(cV, N);
    
    cv0 = cross( z(:,1), n );
    uu = normc( cv0(:,1) );
    
    maxSteps = 8;
    for i = 1:maxSteps
        we1 = ones( numSignals, nV);
        we1(1,:) = we1(1,:) .* uu(1,1);
        we1(2,:) = we1(2,:) .* uu(2,1);
        we1(3,:) = we1(3,:) .* uu(3,1);

        dwe = sign( dot( we, we1 ) );
        we(1,:) = dwe .* we(1,:);
        we(2,:) = dwe .* we(2,:);
        we(3,:) = dwe .* we(3,:);

        uu = mean(we,2);% normal of the plane passing by 4 of the 8 points
        %uu = uu - dot( uu, n ) * n;
        uu = normc(uu);
    end
end

function [ u ] = FindDiagonal( z, n )
% n, vector orienting the counting.
    numSignals = size(z,1);
    numSamples = size(z,2);
    
    N = ones( numSignals, numSamples);
    N(1,:) = N(1,:) .* n(1);
    N(2,:) = N(2,:) .* n(2);
    N(3,:) = N(3,:) .* n(3);
    
    D = sign( dot( z, N ) );
    DD = ones( numSignals, numSamples);
    DD(1,:) = D(1,:) .* DD(1,:);
    DD(2,:) = D(1,:) .* DD(2,:);
    DD(3,:) = D(1,:) .* DD(3,:);
    
    P =  DD .* z;
    u = mean( P, 2 ); 
end