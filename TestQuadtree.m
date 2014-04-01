close all;

% generate the triangles

numVariables = 2;
numSamples = 10000;
numCoeffs = 2;
A = randn( numVariables, numVariables, numCoeffs );
B = randn( 1, numVariables );

ppoints = GenerateUniformSources( numVariables, numSamples );
tppoints = ApplyDistortionR( A, B, ppoints );
figure;
PlotSignal( tppoints );
title( 'Points' );

pq = PointQuadTree( 32, tppoints );
DisplayPointQuadTree( pq );

triangles = CreateMesh( 32, 32 );
ttriangles = ApplyDistortionR( A, B, triangles );

N = PointQuadTreeCompterPoints( pq, ttriangles );

M = 0;
numTriangles = floor( size( triangles, 2 ) / 3 );

for i=1:numSamples
    point = tppoints(:,i);
    
    for j=1:numTriangles
        
        A = ttriangles(:, (j-1)*3 + 1);
        B = ttriangles(:, (j-1)*3 + 2);
        C = ttriangles(:, (j-1)*3 + 3);

        % we test if the point is in the triangle
        v1 = A - B;
        v2 = C - B;
        tp = point - B;

        if( det( [ v1 v2 ] ) ~= 0 )
            uv = [ v1 v2 ] \ tp;
            uv = uv - [ 0.5; 0.5 ];

            if( (abs(uv(1)) <= 0.5) && (abs(uv(2)) <= 0.5) && (uv(1)+uv(2) <= 0) )
                M = M + 1;%point is in triangle
                break;
            end
        end
    end
end