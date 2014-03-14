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
    B0 = m;% A1+A2+A3
    n0 = normc(c);%normalized vector in one of the separating planes
    
    B = zeros( numSignals, 4);
    B(:,1) = B0;
    N= zeros( numSignals, 4 );
    N(:,1)= n0;
    
    [ U, O, V ] = SeparatePoints( z, c, 6.0 * s );
    for i=1:3
        
        figure;
        PlotSignal( V );
        title('V avant création du plan');
        
        n = FindPlaneNormal( V, n0 );
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
    
    figure;
    PlotSignal( z, 'b' );
    hold on;
    PlotVectors([ N(:,2) N(:,3) N(:,4) ], 'r' );
    hold on;
    PlotVectors([ B(:,2) B(:,3) B(:,4) ], 'b' );
    title('regions');
    W = B(:,2:4);
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

function [ uu ] = FindPlaneNormal( z, n )
% Find the normal of a plane passing by the center of the cube and
% containing the vector n
    numSignals = size( z, 1 );
    nV = size( z, 2 );
    
    cV = z;
    cV = normc(cV);
    
    N = ones( numSignals, nV);
    N(1,:) = N(1,:) .* n(1);
    N(2,:) = N(2,:) .* n(2);
    N(3,:) = N(3,:) .* n(3);
    
    we = cross(cV, N);
    uu = normc( we(:,1) );
    
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