classdef PointQuadTree
    %POINTQUADTREE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        aabb;
        Triangles;
        Points;
        SubNodes;
    end
    
    methods
        
        % indices are triangle indices, not point indices
        function q  = PointQuadTree( pointsPerNode, pointList, maxDepth, curDepth )
            
            if( nargin ~= 0 )
                if( nargin <= 3 )
                    curDepth = 1;
                end
                if( nargin == 2 )
                    maxDepth = 7;
                end


                numPoints = size( pointList, 2 );

                q.aabb = AABB( pointList );
                q.Triangles = zeros(2,0);
                q.SubNodes = PointQuadTree.empty(4,0);

                if( (numPoints <= pointsPerNode) || (maxDepth == curDepth) )
                    q.Points = pointList;
                else

                    Ur = zeros(2, numPoints);
                    nUr = 0;
                    Ul = zeros(2, numPoints);
                    nUl = 0;
                    Br = zeros(2, numPoints);
                    nBr = 0;
                    Bl = zeros(2, numPoints);
                    nBl = 0;

                    for i=1:numPoints

                        if( ((pointList(1, i ) >= q.aabb.Center(1)) && (pointList(2, i ) >= q.aabb.Center(2))) )
                             nUr = nUr + 1;
                             Ur(:, nUr ) = pointList(:,i);
                        elseif( ((pointList(1, i ) >= q.aabb.Center(1)) && (pointList(2, i ) <= q.aabb.Center(2))) )
                             nBr = nBr + 1;
                             Br(:, nBr ) = pointList(:,i);
                        elseif( ((pointList(1, i ) <= q.aabb.Center(1)) && (pointList(2, i ) <= q.aabb.Center(2))) )
                             nBl = nBl + 1;
                             Bl(:, nBl ) = pointList(:,i);
                        elseif( ((pointList(1, i ) <= q.aabb.Center(1)) && (pointList(2, i ) >= q.aabb.Center(2))) )
                             nUl = nUl + 1;
                             Ul(:, nUl ) = pointList(:,i);
                        end
                    end

                    Ur = Ur(:, 1:nUr);
                    q.SubNodes(1) = PointQuadTree( pointsPerNode, Ur, maxDepth, curDepth + 1 );

                    Ul = Ul(:, 1:nUl);
                    q.SubNodes(2) = PointQuadTree( pointsPerNode, Ul, maxDepth, curDepth + 1 );

                    Br = Br(:, 1:nBr);
                    q.SubNodes(4) = PointQuadTree( pointsPerNode, Br, maxDepth, curDepth + 1 );

                    Bl = Bl(:, 1:nBl);
                    q.SubNodes(3) = PointQuadTree( pointsPerNode, Bl, maxDepth, curDepth + 1 );

                    q.Points = zeros(0,0);
                end
            end
        end
        
        function N = PointQuadTreeCompterPoints( q, triangleList )
            
            N = 0;
            numTriangles = floor( size( triangleList, 2 ) / 3);
            if( size(q.Points,1) ~= 0 )
                 % leaf node: must direclty count the points in the triangles
                 numPoints = size(q.Points, 2 );
                 
                 for l =1:numPoints
                     point = q.Points(:,l);
                     
                     for i=1:numTriangles

                        A = triangleList(:, (i-1)*3 + 1);
                        B = triangleList(:, (i-1)*3 + 2);
                        C = triangleList(:, (i-1)*3 + 3);

                        % we test if the point is in the triangle
                        v1 = A - B;
                        v2 = C - B;
                        tp = point - B;

                        if( det( [ v1 v2 ] ) ~= 0 )
                            uv = [ v1 v2 ] \ tp;
                            uv = uv - [ 0.5; 0.5 ];

                            if( (abs(uv(1)) <= 0.5) && (abs(uv(2)) <= 0.5) && (uv(1)+uv(2) <= 0) )
                                N = N + 1;%point is in triangle
                                break;
                            end
                        end
                     end
                 end
                 
                 
            else
                
                % must make lists of triangles for each subcell
                Ur = zeros(2, numTriangles*3);
                nUr = 0;
                Ul = zeros(2, numTriangles*3);
                nUl = 0;
                Br = zeros(2, numTriangles*3);
                nBr = 0;
                Bl = zeros(2, numTriangles*3);
                nBl = 0;

                for i=1:numTriangles

                    if( ((triangleList(1, (i-1)*3 + 1 ) >= q.aabb.Center(1)) && (triangleList(2, (i-1)*3 + 1 ) >= q.aabb.Center(2))) || ((triangleList(1, (i-1)*3 + 2 ) >= q.aabb.Center(1)) && (triangleList(2, (i-1)*3 + 2 ) >= q.aabb.Center(2))) || ((triangleList(1, (i-1)*3 + 3 ) >= q.aabb.Center(1)) && (triangleList(2, (i-1)*3 + 3 ) >= q.aabb.Center(2))) )
                         nUr = nUr + 1;
                         Ur(:, (nUr-1)*3 + 1 ) = triangleList(:, (i-1)*3 + 1 );
                         Ur(:, (nUr-1)*3 + 2 ) = triangleList(:, (i-1)*3 + 2 );
                         Ur(:, (nUr-1)*3 + 3 ) = triangleList(:, (i-1)*3 + 3 );
                    end
                    if( ((triangleList(1, (i-1)*3 + 1 ) >= q.aabb.Center(1)) && (triangleList(2, (i-1)*3 + 1 ) <= q.aabb.Center(2))) || ((triangleList(1, (i-1)*3 + 2 ) >= q.aabb.Center(1)) && (triangleList(2, (i-1)*3 + 2 ) <= q.aabb.Center(2))) || ((triangleList(1, (i-1)*3 + 3 ) >= q.aabb.Center(1)) && (triangleList(2, (i-1)*3 + 3 ) <= q.aabb.Center(2))) )
                         nBr = nBr + 1;
                         Br(:, (nBr-1)*3 + 1 ) = triangleList(:, (i-1)*3 + 1 );
                         Br(:, (nBr-1)*3 + 2 ) = triangleList(:, (i-1)*3 + 2 );
                         Br(:, (nBr-1)*3 + 3 ) = triangleList(:, (i-1)*3 + 3 );
                    end
                    if( ((triangleList(1, (i-1)*3 + 1 ) <= q.aabb.Center(1)) && (triangleList(2, (i-1)*3 + 1 ) <= q.aabb.Center(2))) || ((triangleList(1, (i-1)*3 + 2 ) <= q.aabb.Center(1)) && (triangleList(2, (i-1)*3 + 2 ) <= q.aabb.Center(2))) || ((triangleList(1, (i-1)*3 + 3 ) <= q.aabb.Center(1)) && (triangleList(2, (i-1)*3 + 3 ) <= q.aabb.Center(2))) )
                         nBl = nBl + 1;
                         Bl(:, (nBl-1)*3 + 1 ) = triangleList(:, (i-1)*3 + 1 );
                         Bl(:, (nBl-1)*3 + 2 ) = triangleList(:, (i-1)*3 + 2 );
                         Bl(:, (nBl-1)*3 + 3 ) = triangleList(:, (i-1)*3 + 3 );
                    end
                    if( ((triangleList(1, (i-1)*3 + 1 ) <= q.aabb.Center(1)) && (triangleList(2, (i-1)*3 + 1 ) >= q.aabb.Center(2))) || ((triangleList(1, (i-1)*3 + 2 ) <= q.aabb.Center(1)) && (triangleList(2, (i-1)*3 + 2 ) >= q.aabb.Center(2))) || ((triangleList(1, (i-1)*3 + 3 ) <= q.aabb.Center(1)) && (triangleList(2, (i-1)*3 + 3 ) >= q.aabb.Center(2))) )
                         nUl = nUl + 1;
                         Ul(:, (nUl-1)*3 + 1 ) = triangleList(:, (i-1)*3 + 1 );
                         Ul(:, (nUl-1)*3 + 2 ) = triangleList(:, (i-1)*3 + 2 );
                         Ul(:, (nUl-1)*3 + 3 ) = triangleList(:, (i-1)*3 + 3 );
                    end
                end

                Ur = Ur(:, 1:nUr*3);
                Ul = Ul(:, 1:nUl*3);
                Br = Br(:, 1:nBr*3);
                Bl = Bl(:, 1:nBl*3);
                
                if( nUr ~= 0 )
                    N = PointQuadTreeCompterPoints( q.SubNodes(1), Ur );
                end
                
                if( nUl ~= 0 )
                    N = N + PointQuadTreeCompterPoints( q.SubNodes(2), Ul );
                end
                
                if( nBr ~= 0 )
                    N = N + PointQuadTreeCompterPoints( q.SubNodes(4), Br );
                end
                
                if( nBl ~= 0 )
                    N = N + PointQuadTreeCompterPoints( q.SubNodes(3), Bl );
                end
                
            end
        end
        
        function DisplayPointQuadTree(q, color)
            
            if( nargin  == 1 )
                color = [ 1 0 0 ];
            end
            
            a = q.aabb;
            di = a.HalfDiagonal;
            if( any(di) )
                hold on;
                DisplayAABB(q.aabb, color);

                if( size(q.SubNodes) ~= 0 )
                    DisplayPointQuadTree(q.SubNodes(1), 0.9 * color);
                    DisplayPointQuadTree(q.SubNodes(2), 0.9 * color);
                    DisplayPointQuadTree(q.SubNodes(3), 0.9 * color);
                    DisplayPointQuadTree(q.SubNodes(4), 0.9 * color);
                end
            else
                title('Points');
            end
        end
        
        function depth = PointQuadTreeDepth( q )
            subDepths = zeros(1,4);
            if( size(q.SubNodes) ~= 0 )
                subDepths(1) = PointQuadTreeDepth(q.SubNodes(1));
                subDepths(2) = PointQuadTreeDepth(q.SubNodes(2));
                subDepths(3) = PointQuadTreeDepth(q.SubNodes(3));
                subDepths(4) = PointQuadTreeDepth(q.SubNodes(4));
            end
            
            depth = 1 + max( subDepths );
        end
        
    end
    
end

