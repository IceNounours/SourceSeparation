close all;
clear;

% demonstration settings
numVariables = 2;
numSamples = 10000;

% generate the signals
sigma = 0.1;
s1 = sigma;
s2 = sigma;
p = 0.5;
q = 0.5;
s = GenerateBernouilliSources( numVariables, numSamples, [ p q ] );
n = sigma .* randn( numVariables, numSamples );

A = randn( numVariables );

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

figure;
PlotSignal( A\x );
title('Original transformed back signal with noise');

% find matrix back
W  = CoinsLikelihood( x, p, q, sigma, sigma );
B = A \ W;

z = W \ x;
figure;
PlotSignal(z);
title('Transformed back signal');
