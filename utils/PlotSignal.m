function [ ] = PlotSignal( s, params )
%PLOTSIGNAL Plot the signal s.
%   We suppose that s is a 2D or 3D signal

    D = size( s, 1 );
    
    if nargin == 1
        params = '.';
    else
        params = strcat( params, '.' ); 
    end
    
    if( D == 2 )
        plot( s(1,:), s(2, :), params );
    elseif( D == 3)
        plot3( s(1,:),s(2, :),s(3,:), params );   
    end
    
    grid on;
    axis square;
end

