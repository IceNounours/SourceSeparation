close all;
clear;

% demonstration settings
numVariables = 2;
numSamples = 10000;

% generate the signals
s = GenerateUniformSources( numVariables, numSamples );

A = randn( numVariables );

% mixed signal
x = A*s;

% plot signals
figure;
PlotSignal( s );
title('Original signal');

figure;
PlotSignal( A*s );
title('Transformed signal');

% find matrix back
W  = CoinsU2( x );
B = A \ W;

z = W \ x;
figure;
axis equal;
PlotSignal(x);
hold on;
axis equal;
PlotVectors(W);
title('Vectors');
