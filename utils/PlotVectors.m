function [ ] = PlotVectors( W, params, ps )
%PLOTVECTORS Plot the vectors of s.
%   We suppose that s is a 2D or 3D signal

    O = zeros(1, size(W,2));
    
    if(nargin == 1)
        params = 'r';
        xs = O;
        ys = O;
        zs = O;
    elseif(nargin == 2 )
        xs = O;
        ys = O;
        zs = O;
    else
        xs = ps(1,:);
        ys = ps(2,:);
        zs = ps(3,:);
    end
    
    D = size( W, 1 );
    
    if( D == 2 )
        quiver( xs, ys, W(1,:), W(2, :), params, 'LineWidth', 1 );
    elseif( D == 3)
        quiver3( xs, ys, zs, W(1,:),W(2, :),W(3,:), params, 'LineWidth', 1 );   
    end
    
    grid on;
end

