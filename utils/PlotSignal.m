function [ ] = PlotSignal( s )
%PLOTSIGNAL Plot the signal s.
%   We suppose that s is a 2D or 3D signal

    D = size( s, 1 );
    
    if( D == 2 )
        plot( s(1,:), s(2, :), '.' );
    elseif( D == 3)
        plot3( s(1,:),s(2, :),s(3,:), '.' );   
    end
    
    grid on;
    axis square;
end

