close all;

% generate the triangles

profile on;
for i=1:10
    points = CreateMesh( 32*i, 32*i );
    grid = TransformedGrid( points );
end
profile viewer;


points = CreateMesh( 32, 32 );
figure;
PlotSignal( points );
title( 'Triangles' );

grid = TransformedGrid( points );
q = grid.quadTree;
depth = QuadTreeDepth( q );
DisplayQuadTree( q );

isIn1 = TransformedGridIsPointIn( grid, [2; 0] );
isIn2 = TransformedGridIsPointIn( grid, [0.5; 0.5 ]);