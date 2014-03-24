classdef AABB
    
    properties
        Center       % Center of the AABB
        HalfDiagonal % half diagonal of the AABB
    end
    
    methods
        
        function obj = AABB( pointList )
            numPoints = size(pointList, 2 );
            
            if( numPoints ~= 0 )
                maxv = max( pointList, [], 2 );
                minv = min( pointList, [], 2 );
                
                obj.Center = 0.5 * (maxv + minv);
                obj.HalfDiagonal = 0.5 * (maxv - minv);
            else
                obj.Center = zeros( 2, 1);
                obj.HalfDiagonal = zeros( 2, 1 );
            end
            
        end
        
        function DisplayAABB(aabb, color)
            
            line( [ aabb.Center(1)-aabb.HalfDiagonal(1), aabb.Center(1)-aabb.HalfDiagonal(1) ], [ aabb.Center(2)-aabb.HalfDiagonal(2), aabb.Center(2)+aabb.HalfDiagonal(2) ], 'Color',color);
            line( [ aabb.Center(1)+aabb.HalfDiagonal(1), aabb.Center(1)+aabb.HalfDiagonal(1) ], [ aabb.Center(2)-aabb.HalfDiagonal(2), aabb.Center(2)+aabb.HalfDiagonal(2) ], 'Color',color);
            line( [ aabb.Center(1)-aabb.HalfDiagonal(1), aabb.Center(1)+aabb.HalfDiagonal(1) ], [ aabb.Center(2)+aabb.HalfDiagonal(2), aabb.Center(2)+aabb.HalfDiagonal(2) ], 'Color',color);
            line( [ aabb.Center(1)-aabb.HalfDiagonal(1), aabb.Center(1)+aabb.HalfDiagonal(1) ], [ aabb.Center(2)-aabb.HalfDiagonal(2), aabb.Center(2)-aabb.HalfDiagonal(2) ], 'Color',color);
        end
        
    end
end

