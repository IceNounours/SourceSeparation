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

        steps = steps + 1;
        e = m - (A1+A2);
        
        scale = scale * 0.9;
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
end

function [ v ] = perp(x)
    v(1) = x(2);
    v(2) = -x(1);
end