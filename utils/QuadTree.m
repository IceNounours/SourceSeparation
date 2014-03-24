classdef QuadTree
    %QUADTREE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        aabb;
        Data;
        SubNodes;
    end
    
    methods
        
        % indices are triangle indices, not point indices
        function q  = QuadTree( trianglesPerNode, triangleList, maxDepth, curDepth )
            
            if( nargin ~= 0 )
                if( nargin <= 3 )
                    curDepth = 1;
                end
                if( nargin == 2 )
                    maxDepth = 7;
                end


                numTriangles = floor( size( triangleList, 2 ) / 3 );

                q.aabb = AABB( triangleList );

                q.SubNodes = QuadTree.empty(4,0);

                if( (numTriangles <= trianglesPerNode) || (maxDepth == curDepth) )
                    q.Data = triangleList;
                else

                    Ur = zeros(2, numTriangles*3);
                    nUr = 0;
                    Ul = zeros(2, numTriangles*3);
                    nUl = 0;
                    Br = zeros(2, numTriangles*3);
                    nBr = 0;
                    Bl = zeros(2, numTriangles*3);
                    nBl = 0;

                    O = zeros(2, numTriangles*3);
                    nO = 0;

                    for i=1:numTriangles

                        if( ((triangleList(1, (i-1)*3 + 1 ) >= q.aabb.Center(1)) && (triangleList(2, (i-1)*3 + 1 ) >= q.aabb.Center(2))) && ((triangleList(1, (i-1)*3 + 2 ) >= q.aabb.Center(1)) && (triangleList(2, (i-1)*3 + 2 ) >= q.aabb.Center(2))) && ((triangleList(1, (i-1)*3 + 3 ) >= q.aabb.Center(1)) && (triangleList(2, (i-1)*3 + 3 ) >= q.aabb.Center(2))) )
                             nUr = nUr + 1;
                             Ur(:, (nUr-1)*3 + 1 ) = triangleList(:, (i-1)*3 + 2 );
                             Ur(:, (nUr-1)*3 + 2 ) = triangleList(:, (i-1)*3 + 2 );
                             Ur(:, (nUr-1)*3 + 3 ) = triangleList(:, (i-1)*3 + 3 );
                        elseif( ((triangleList(1, (i-1)*3 + 1 ) >= q.aabb.Center(1)) && (triangleList(2, (i-1)*3 + 1 ) <= q.aabb.Center(2))) && ((triangleList(1, (i-1)*3 + 2 ) >= q.aabb.Center(1)) && (triangleList(2, (i-1)*3 + 2 ) <= q.aabb.Center(2))) && ((triangleList(1, (i-1)*3 + 3 ) >= q.aabb.Center(1)) && (triangleList(2, (i-1)*3 + 3 ) <= q.aabb.Center(2))) )
                             nBr = nBr + 1;
                             Br(:, (nBr-1)*3 + 1 ) = triangleList(:, (i-1)*3 + 2 );
                             Br(:, (nBr-1)*3 + 2 ) = triangleList(:, (i-1)*3 + 2 );
                             Br(:, (nBr-1)*3 + 3 ) = triangleList(:, (i-1)*3 + 3 );
                        elseif( ((triangleList(1, (i-1)*3 + 1 ) <= q.aabb.Center(1)) && (triangleList(2, (i-1)*3 + 1 ) <= q.aabb.Center(2))) && ((triangleList(1, (i-1)*3 + 2 ) <= q.aabb.Center(1)) && (triangleList(2, (i-1)*3 + 2 ) <= q.aabb.Center(2))) && ((triangleList(1, (i-1)*3 + 3 ) <= q.aabb.Center(1)) && (triangleList(2, (i-1)*3 + 3 ) <= q.aabb.Center(2))) )
                             nBl = nBl + 1;
                             Bl(:, (nBl-1)*3 + 1 ) = triangleList(:, (i-1)*3 + 2 );
                             Bl(:, (nBl-1)*3 + 2 ) = triangleList(:, (i-1)*3 + 2 );
                             Bl(:, (nBl-1)*3 + 3 ) = triangleList(:, (i-1)*3 + 3 );
                        elseif( ((triangleList(1, (i-1)*3 + 1 ) <= q.aabb.Center(1)) && (triangleList(2, (i-1)*3 + 1 ) >= q.aabb.Center(2))) && ((triangleList(1, (i-1)*3 + 2 ) <= q.aabb.Center(1)) && (triangleList(2, (i-1)*3 + 2 ) >= q.aabb.Center(2))) && ((triangleList(1, (i-1)*3 + 3 ) <= q.aabb.Center(1)) && (triangleList(2, (i-1)*3 + 3 ) >= q.aabb.Center(2))) )
                             nUl = nUl + 1;
                             Ul(:, (nUl-1)*3 + 1 ) = triangleList(:, (i-1)*3 + 2 );
                             Ul(:, (nUl-1)*3 + 2 ) = triangleList(:, (i-1)*3 + 2 );
                             Ul(:, (nUl-1)*3 + 3 ) = triangleList(:, (i-1)*3 + 3 );
                        else
                             nO = nO + 1;
                             O(:, (nO-1)*3 + 1 ) = triangleList(:, (i-1)*3 + 1 );
                             O(:, (nO-1)*3 + 2 ) = triangleList(:, (i-1)*3 + 2 );
                             O(:, (nO-1)*3 + 3 ) = triangleList(:, (i-1)*3 + 3 );
                        end
                    end

                    Ur = Ur(:, 1:nUr*3);
                    q.SubNodes(1) = QuadTree( trianglesPerNode, Ur, maxDepth, curDepth + 1 );

                    Ul = Ul(:, 1:nUl*3);
                    q.SubNodes(2) = QuadTree( trianglesPerNode, Ul, maxDepth, curDepth + 1 );

                    Br = Br(:, 1:nBr*3);
                    q.SubNodes(4) = QuadTree( trianglesPerNode, Br, maxDepth, curDepth + 1 );

                    Bl = Bl(:, 1:nBl*3);
                    q.SubNodes(3) = QuadTree( trianglesPerNode, Bl, maxDepth, curDepth + 1 );

                    O = O(:,1:nO*3);
                    q.Data = O;
                end
            end
        end
        
        function DisplayQuadTree(q, color)
            
            if( nargin  == 1 )
                color = [ 1 0 0 ];
            end
            
            hold on;
            DisplayAABB(q.aabb, color);
            
            if( size(q.SubNodes) ~= 0 )
                DisplayQuadTree(q.SubNodes(1), 0.9 * color);
                DisplayQuadTree(q.SubNodes(2), 0.9 * color);
                DisplayQuadTree(q.SubNodes(3), 0.9 * color);
                DisplayQuadTree(q.SubNodes(4), 0.9 * color);
            end
        end
        
        function depth = QuadTreeDepth( q )
            subDepths = zeros(1,4);
            if( size(q.SubNodes) ~= 0 )
                subDepths(1) = QuadTreeDepth(q.SubNodes(1));
                subDepths(2) = QuadTreeDepth(q.SubNodes(2));
                subDepths(3) = QuadTreeDepth(q.SubNodes(3));
                subDepths(4) = QuadTreeDepth(q.SubNodes(4));
            end
            
            depth = 1 + max( subDepths );
        end
        
        function canBeIn = QuadTreePointCanBeIn( q, point )
            
            delta = point - q.aabb.Center;
            canBeIn =  abs(delta(1)) <= q.aabb.HalfDiagonal(1) && abs(delta(2)) <= q.aabb.HalfDiagonal(2);
        end
        
        % we assume that the point can be in.
        function isIn = QuadTreeIsPointIn( q, point )
            
            numTriangles = size(q.Data, 2 );
            for i=1:numTriangles

                A = q.Data(:, (i-1)*3 + 1);
                B = q.Data(:, (i-1)*3 + 2);
                C = q.Data(:, (i-1)*3 + 3);

                % we test if the point is in the triangle
                v1 = A - B;
                v2 = C - B;
                tp = point - B;

                uv = [ v1 v2 ] \ tp;
                uv = uv - [ 0.5; 0.5 ];

                isIn = (abs(uv(1)) <= 0.5) && (abs(uv(2)) <= 0.5) && (uv(1)+uv(2) <= 0);

                if( isIn == true )
                    break;
                end
            end

            if( size(q.SubNodes) ~= 0 )
                for i=1:4
                	canBeIn = QuadTreePointCanBeIn(q.SubNodes(i), point);
                    if( canBeIn )
                        isIn = QuadTreeIsPointIn( q.SubNodes(i), point );
                        break;% can't be in another node
                    end
                end
            end
            
        end
        
    end
    
end

