classdef TransformedGrid
    %TRANSFORMEDGRID Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        quadTree
        surface;
    end
    
    methods
        
        function grid = TransformedGrid( triangles, trianglesPerNode )
            
            if( nargin == 1 )
                trianglesPerNode = 32;% by default, maximum 32 triangles per node
            end
            
            grid.surface = ComputeSurface( triangles );
            grid.quadTree = QuadTree( trianglesPerNode, triangles, 5 );
        end
        
        function isIn = TransformedGridIsPointIn( grid, point )
            canBeIn = QuadTreePointCanBeIn( grid.quadTree, point );
            if( canBeIn )
                isIn = QuadTreeIsPointIn( grid.quadTree, point );
            else
                isIn = false;
            end
        end
        
    end
    
end

function s = ComputeSurface( triangles )

    numTriangles = floor( size( triangles, 2 ) / 3 );

    alls = zeros(1, numTriangles);
    for i=1:numTriangles

        tl = triangles(:, (i-1)*3 + 1 : (i-1)*3 + 3);
        v1 = tl(:,1)-tl(:,2);
        v2 = tl(:,3)-tl(:,2);

        alls(i) = v1(1)*v2(2) - v1(2)*v2(1);
    end
    
    s = sum(abs(alls));
end

