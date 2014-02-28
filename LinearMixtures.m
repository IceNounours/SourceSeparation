close all;
clear;

% Demonstration settings
numVariables = 2;
numSamples = 10000;
fNumSamples = single(numSamples);
fInvNumSamples = 1/fNumSamples;

% Generate the signals
s = GenerateUniformSources( numVariables, numSamples );
Cs = s*s'*fInvNumSamples;

% Wether or not we want to show the original data
showOriginalSignals = 1;
if( showOriginalSignals )
    figure;
    PlotSignal( s );
    title('Original signal');
end

% Mixture
A = [ 1, 0.5; 1, 1 ];

x = GenerateLinearMixture( A, s );
Cx = x*x'*fInvNumSamples;

% Wether or not we want to show the mixture
showMixture= 1;
if( showMixture )
    figure;
    PlotSignal( x );
    title('Mixture');
end

% Whitening
z = Whiten( x );
Cz = z*z'*fInvNumSamples;

showWhitened = 1;
if( showWhitened )
    figure;
    PlotSignal(z);
    title('Whitened mixture');
end
