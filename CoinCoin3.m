close all;
clear;

% demonstration settings
numVariables = 3;
numSamples = 10000;

% generate the signals
sigma = 0.01;
s1 = sigma;
s2 = sigma;
s3 = sigma;
p = 0.5;
q = 0.5;
r = 0.5;
s = GenerateBernouilliSources( numVariables, numSamples, [ p q r ] );
n = sigma .* randn( numVariables, numSamples );

A = randn( numVariables );
A = [ 1 0 0; 0 1 0; 0 0 1 ];

% mixed signal
x = A*s + n;

% plot signals
figure;
PlotSignal( s );
title('Original signal');

figure;
PlotSignal( A*s );
title('Transformed signal without noise');

figure;
PlotSignal( x );
title('Mixed signal with noise');

% find matrix back
W  = Coin3( x, p, q, r, [ sigma sigma sigma ] );