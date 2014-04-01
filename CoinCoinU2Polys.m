close all;
clear;

% demonstration settings
numVariables = 2;
numSamples = 1000;
numCoeffs = 3;

% generate the signals
s = GenerateUniformSources( numVariables, numSamples );

A = randn( numVariables, numVariables, numCoeffs );
A( 1, 1, 3 ) = 0;
A( 1, 2, 3 ) = 0;
A( 2, 1, 3 ) = 0;
A( 2, 2, 3 ) = 0;

B = randn( 1, numVariables );

% mixed signal
x = ApplyDistortion(A, s);

% plot signals
figure;
PlotSignal( s );
title('Original signal');

figure;
PlotSignal( x );
title('Transformed signal');

% find matrix back
[ P]  = CoinsU2PolysFast( x, numCoeffs );

figure;
PlotSignal( ApplyDistortion(P, s) );
title('Original mixed with the guesssed coefficients');