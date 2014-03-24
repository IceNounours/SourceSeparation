close all;
clear;

% demonstration settings
numVariables = 2;
numSamples = 1000;
numCoeffs = 3;

% generate the signals
s = GenerateUniformSources( numVariables, numSamples );

A = randn( numVariables, numVariables, numCoeffs );

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
W  = CoinsU2PolysFast( x, numCoeffs );

figure;
PlotSignal( ApplyDistortion(W, s) );
title('Original mixed with the gusssed coefficients');