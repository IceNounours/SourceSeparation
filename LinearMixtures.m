close all;
clear;

% Demonstration settings
numVariables = 3;
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
A = [ 1, 0.5, 0.5; 0, 1, 0; 0, 0, 1 ];

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
    axis equal;
    title('Whitened mixture');
end

% R�cup�ration des vecteurs de la nouvelle base
W = Kurt( z );
figure;
PlotSignal(z);
hold on;
PlotVectors(W);
title('Basis vector');

% R�cup�ration des vecteurs de la nouvelle base
Wb = KurtSymmetricOrtho( z );
figure;
PlotSignal(z);
hold on;
PlotVectors(Wb);
title('Basis vector (Symmetric orthogonalization)');

% R�cup�ration des vecteurs de la nouvelle base
Wc = Negen( z );
figure;
PlotSignal(z);
hold on;
PlotVectors(Wc);
title('Basis vector Negentropy');

% R�cup�ration des vecteurs de la nouvelle base
Wd = NegenSymmetricOrtho( z );
figure;
PlotSignal(z);
hold on;
PlotVectors(Wd);
title('Basis vector Negentropy (Symmetric orthogonalization)');

% Donn�es d�mix�es
x = W' *  z;
figure;
PlotSignal(x);
title('Donnees demixees');

% Donn�es d�mix�es
xb = Wb' *  z;
figure;
PlotSignal(xb);
title('Donnees demixees (Symmetric orthogonalization)');

% Donn�es d�mix�es
xc = Wc' *  z;
figure;
PlotSignal(xc);
title('Donnees demixees Negentropy');

% Donn�es d�mix�es
xd = Wd *  z;
figure;
PlotSignal(xd);
title('Donnees demixees Negentropy (Symmetric orthogonalization)');
