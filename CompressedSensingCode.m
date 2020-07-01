close all;
clear all;
clc;

N = 1000;
k = 150;  % As k < N
Fs = 200;
t = (0:N-1)*1/Fs;
freq = 25;
x = sin(2*pi*freq*t);
figure;
subplot(3,2,1);
plot(t, x); title('Sine wave'), xlim([0, 1.3]);
xlabel("Time"), ylabel("Amplitude")

% For calculating the fourier transform, Orthogonal basis matrix
discreteFourierTransform = dftmtx(N);
inverseFourierTransformMatrix = conj(discreteFourierTransform)/N;

% fourier transform
X = discreteFourierTransform*x'; % sparse signal
subplot(3,2,2);
freq1 = -Fs/2:Fs/length(X):Fs/2-(Fs/length(X));
plot(freq1, abs(fftshift(X))/N); title('X(t)');
xlabel("Frequency"), ylabel("Amplitude");

% Generating measurement matrix to show how it looks on 
% plot after cs acquisition

% set RNG seed,  initializes the Mersenne Twister generator 
% using a seed of 0.
rng(0);

% initialising the array
measurmentValues = zeros(N,1);

% taking random values till N, this is important
randomMeasures = randperm(N); 
randomMeasures = randomMeasures(1:k);
measurmentValues(randomMeasures) = x(randomMeasures)';
subplot(3,2,3);
plot(measurmentValues); title('Randomly chosen values'); xlim([1 , k]);
xlabel("Number of measurements"), ylabel("Amplitude");

% Generating sensing matrix
% M x N sensing matrix
A = inverseFourierTransformMatrix(randomMeasures, :);

% Generating measurement values/ observations
y = A*X; % measurement values.
subplot(3,2,4);
plot(abs(y)); title('Compressed value/Measurement values (A * X(t))'); xlim([1 , k]);
xlabel("Number of measurements"), ylabel("Amplitude");

% Recovery of the signal begins from here
% Using the formula y = Ax0, calculating x0
x0 = A'*y; % initial guess

% Using l1 norm(uses primal dual algorithm) to get the sparse signal
% x0: the initial guess, 
% A: our measurement matrix,
% [] : third argument which is used to specify the transpose 
% of the measurement matrix is left empty as we are providing A explicitly,
% y:  is our measured values, 
% 1e-3: is the precission to which we want to solve the problem
XSparseRecovered = l1eq_pd(x0, A, [], y, 1e-3);
subplot(3,2,5);
plot(freq1, abs(fftshift(XSparseRecovered))/N);
title('Recovered sparse signal using L1 magic box');
xlabel("Frequency"), ylabel("Amplitude");

% Using transform matrix to get the original signal out of the sparsed
% signal
originalSIgnal = inverseFourierTransformMatrix*XSparseRecovered;
subplot(3,2,6);
plot(t, abs(originalSIgnal)');
title('Recovered signal'),  xlim([0, 1.3]);
xlabel("Time"), ylabel("Amplitude");

