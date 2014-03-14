function [ W ] = CoinsLikelihood( x, p, q, s1, s2 )
%COINSLIKELIHOOD Summary of this function goes here
%   Detailed explanation goes here
    numSamples = size( x, 2 );
    numSignals = size( x, 1 );
    
    
    c = mean( x, 2 );% center of the rect
    n = perp(c);
    n = normc(n);
    
    m = 2 .* c;% diagonal of the rect.
    
    % find a starting point
    W = zeros( numSignals );
    
    % Bl: list of the points near the origin "bottom left"
    Bl = zeros( numSignals, numSamples/4);
    nBl = 0;
    
    % Ur: list of points near A1+A2 "up right"
    Ur = zeros( numSignals, numSamples/4);
    nUr = 0;
    
    % Ul: list of the points near the "up left"
    Ul = zeros( numSignals, numSamples/4);
    nUl = 0;
    % Br: list of the points near the "bottom right"
    Br = zeros( numSignals, numSamples/4);
    nBr = 0;
    
    for i=1:numSamples
        xx = x(1,i);
        xy = x(2,i);
        % test if the point is near the origin
        if ( (xx/(5*s1))^2 + (xy/(5*s2))^2 < 1 )
            Bl(:,nBl+1) = x(:,i);
            nBl = nBl + 1;
        % test if the point is near A1+A2
        elseif ( ((xx-m(1))/(5*s1))^2 + ((xy-m(2))/(5*s2))^2 < 1 )
            Ur(:,nUr+1) = x(:,i);
            nUr = nUr + 1;
        % the point is one of the others.
        else
            xc = x(:,i) - c;
            
            % the point is near "bottom right"
            if( dot(xc,n) < 0 )
                Br(:,nBr+1) = x(:,i);
                nBr = nBr + 1;
            else
                Ul(:,nUl+1) = x(:,i);
                nUl = nUl + 1;
            end
        end
    end
    
    % shrink the lists to the actual size.
    Bl = Bl(:,1:nBl);
    Br = Br(:,1:nBr);
    Ul = Ul(:,1:nUl);
    Ur = Ur(:,1:nUr);
    
    % compute the vectors
    A2 = mean( Br, 2 );
    A1 = mean( Ul, 2 );
    
%     % compute the error and distribute part of it to the vectors
%     e = m - (A1+A2);
%     A1 = A1 + 0.33 * e;
%     A2 = A2 + 0.33 * e;
    
    scale = 4.0;
    steps = 0;
    maxSteps = 100;
    
    e = m - (A1+A2);
    
    while( (steps < maxSteps) && (norm(e) > norm(m)*0.0001) )
        
        prevE = e;
        
        numPointsBr = size(Br,2);
        numPointsUl = size(Ul,2);

        nUl = 0;
        nBr = 0;
    
        if( numPointsBr > numSamples / 100 )
            for j=1:numSamples
                xx = x(1,j);
                xy = x(2,j);

                if ( ((xx-A2(1))/(scale*s1))^2 + ((xy-A2(2))/(scale*s2))^2 < 2 )
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

                if ( ((xx-A1(1))/(scale*s1))^2 + ((xy-A1(2))/(scale*s2))^2 < 2 )
                    Ul(:,nUl+1) = x(:,j);
                    nUl = nUl + 1;
                end
            end
            
            Ul = Ul(:,1:nUl);
        end
        
        % compute the vectors
        A2 = mean( Br, 2 );
        A1 = mean( Ul, 2 );

%         % compute the error and distribute part of it to the vectors
%         e = m - (A1+A2);
%         A1 = A1 + 0.33 * e;
%         A2 = A2 + 0.33 * e;

        steps = steps + 1;
        e = m - (A1+A2);
        
        %if( dot(e,e) < dot(prevE,prevE) )
            scale = scale * 0.9;
        %end
    end
    
    W = [ A1, A2 ];
    
    figure;
    plot( Bl(1,:), Bl(2, :), '.r' );
    hold on;
    plot( Ur(1,:), Ur(2, :), '.g' );
    hold on
    plot( Br(1,:), Br(2, :), '.b' );
    hold on
    plot( Ul(1,:), Ul(2, :), '.y' );
    hold on;
    PlotVectors( W, 'r' );
    hold on;
    PlotVectors( m, 'b' );
    title('Regions');
    
    
%     prevW = zeros( numSignals );
%     
%     U = zeros( numSignals, 1 );
%     nU = 0;
%     D = zeros( numSignals, 1 );
%     nD = 0;
%     for i=1:numSamples
%         
%         xc = x(:,i);
%         xc = normc(xc);
%         d = det( [xc c]);
%         
%         if( d > 0 )
%             nU = nU + 1;
%             U = U + xc;
%         elseif ( d < 0 )
%             nD = nD + 1;
%             D = D + xc;
%         end
%     end
%     
%     U = U / nU;
%     D = D / nD;
%     
%     W = [ U D ];
%     
%     figure;
%     PlotSignal( x );
%     hold on;
%     Plotvectors( W );
%     title('Vecteurs originaux');
%     
%     steps = 0;
%     maxSteps = 100;
%     while( (min( abs(dot(W,prevW)) ) < 0.9999) && (steps < maxSteps) )
%         prevW = W;
%         W(:,1) = W(:,1) + 0.5 * dl( x, p, q, s1, s2, W(:,1), m );
%         W(:,2) = m - W(:,1);
%         steps = steps + 1 ;
%     end
end

function [ v ] = perp(x)
    v(1) = x(2);
    v(2) = -x(1);
end

function [ r ] = pn( x, s )
    r = gaussmf( x./s, [1 0 ] ) ./ s;
end
function [ r ] = px1( x, p, q, s1, s2, a)

    numSamples = size(x, 2 );
    
    x1 = x(1,:);
    x2 = x(2,:);
    xb1 = x1-a(1);
    xb2 = x2-a(2);
    p1 = pn(xb1, s1);
    p2 = pn(xb2, s2);
    r = q*(1-p).*p2.*p1;
    
    for i= 1:numSamples
        if( r(i) < 0.1 )
            r(i) = 0.1;
        end
    end
end

function [ r ] = px2( x, p, q, s1, s2, a, m )

    numSamples = size(x, 2 );
    
    p1 = pn(x(1,:)+a(1)-m(1), s1);
    p2 = pn(x(2,:)-m(2)+a(2), s2);
    r = p*(1-q).*p2.*p1;
    
    for i= 1:numSamples
        if( r(i) < 0.1 )
            r(i) = 0.1;
        end
    end
end

function [ r ] = px3( x, p, q, s1, s2 )

    numSamples = size(x, 2 );
    
    p1 = pn(x(1,:), s1);
    p2 = pn(x(2,:), s2);
    r = p*q.*p1.*p2;
    
    for i= 1:numSamples
        if( r(i) < 0.1 )
            r(i) = 0.1;
        end
    end
end

function [ r ] = px4( x, p, q, s1, s2, m )

    numSamples = size(x, 2 );
    
    p1 = pn(x(1,:)-m(1), s1);
    p2 = pn(x(2,:)-m(2), s2);
    r = (1-p)*(1-q).*p1.*p2;
    
    for i= 1:numSamples
        if( r(i) < 0.1 )
            r(i) = 0.1;
        end
    end
end

function [ r ] = px( x, p, q, s1, s2, a, m )

    p1 = px1( x, p, q, s1, s2, a );
    p2 = px2( x, p, q, s1, s2, a, m );
    p3 = px3( x, p, q, s1, s2 );
    p4 = px4( x, p, q, s1, s2, m );
    
    r = p1 + p2 + p3 + p4;
end

function [ r ] = dpxda1( x, p, q, s1, s2, a, m )

    r = px1( x, p, q, s1, s2, a ) .* (x(1,:) - a(1)) ./ (s1*s1);
    r = r - px2( x, p, q, s1, s2, a, m ) .* (x(1,:) + a(1) - m(1) ) ./ (s1*s1);
end

function [ r ] = dpxda2( x, p, q, s1, s2, a, m )

    r = px1( x, p, q, s1, s2, a ) .* (x(2,:) - a(2)) ./ (s2*s2);
    r = r - px2( x, p, q, s1, s2, a, m ) .* (x(2,:) + a(2) - m(2) ) ./ (s2*s2);
end

function [ G ]= dl( x, p, q, s1, s2, a, m )

    numSamples = size( x, 2 );
    numSignals = size( x, 1 );
    
    R = px( x, p, q, s1, s2, a, m );
    T = zeros( numSignals, numSamples );
    T(1,:) = dpxda1( x, p, q, s1, s2, a, m ) ./ R;
    T(2,:) = dpxda2( x, p, q, s1, s2, a, m ) ./ R;
    
    G = sum( T, 2 );
end