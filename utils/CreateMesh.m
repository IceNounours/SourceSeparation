function [ points ] = CreateMesh( numRows, numColumns )
%CREATEMESH Create a rectangle subdivised with triangles
%   Detailed explanation goes here
    numTriangles = (numRows-1) * (numColumns-1);
    numPoints = numTriangles * 6;

    points = zeros( 2, numPoints );

    dx = 1.0/ (numColumns-1);
    dy = 1.0/ (numRows-1);

    for i=0:(numColumns-2)
        for l=0:(numRows-2)

            points(:, 6 * (l * (numColumns-2) + i) + 1 ) = [ (i+1) * dx, (l+0) * dy ];
            points(:, 6 * (l * (numColumns-2) + i) + 2 ) = [ (i+0) * dx, (l+0) * dy ];
            points(:, 6 * (l * (numColumns-2) + i) + 3 ) = [ (i+0) * dx, (l+1) * dy ];

            points(:, 6 * (l * (numColumns-2) + i) + 4 ) = [ (i+1) * dx, (l+0) * dy ];
            points(:, 6 * (l * (numColumns-2) + i) + 5 ) = [ (i+1) * dx, (l+1) * dy ];
            points(:, 6 * (l * (numColumns-2) + i) + 6 ) = [ (i+0) * dx, (l+1) * dy ];
        end
    end
end

