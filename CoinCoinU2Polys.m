close all;
clear;

% demonstration settings
numVariables = 2;
numSamples = 1000;
numCoeffs = 2;

% generate the signals
s = GenerateUniformSources( numVariables, numSamples );

A = randn( numVariables, ( numVariables * numCoeffs + 1) );

% mixed signal
x = ApplyDistortion(A, s);

figure;
PlotSignal(x);
title('Distorted signal');

P = CoinsU2PolysFast( x, numCoeffs );
gx = ApplyDistortion(P, x );
figure;
PlotSignal(gx);
title('Distorted signal with the guessed coefficients');