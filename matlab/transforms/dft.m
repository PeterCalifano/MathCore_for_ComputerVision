function [G, f] = dft(g, t, Nf)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DISCREET FOURIER TRANSFORM
% g = input signal
% t = input time axis
% Nf = number of frequency samples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The execution time for fft depends on the length of the transform. 
% The time is fastest for powers of two and almost as fast for lengths 
% that have only small prime factors. The time is typically several 
% times slower for lengths that are prime, or which have large prime factors.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if not(exist('Nf', 'var'))
    Nf = length(t);
end

if size(g,1)==1
    g = g.';
end

G = fft(g,Nf,1); % transform along columns by default
G = fftshift(G,1);

% frequency axis
dt = t(2)-t(1);
if mod(Nf,2) % odd number of samples
    f = (-(Nf-1)/2:(Nf-1)/2)/Nf/dt;
else   % even number of samples
    f = (-Nf/2:Nf/2-1)/Nf/dt;
end

% fft delay compensation
phase = -2*pi*f(:)*t(1);
w = exp(1i*phase(:))*ones(1,size(g,2));
G = G.*w;

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% example
clear, clc
t = linspace(-1,1,101);
g = randn(length(t),2) + 1i*randn(length(t),2);

Nf = 2^ceil(log2(length(t)));
[G,f] = dft(g,t,Nf);

W = exp(-1i*2*pi*f(:)*t);
G2 = W*g;

figure, 
subplot(2,1,1), plot(f,abs(G2-G)), grid

g2 = idft(G,f,t);
subplot(2,1,2), plot(t,abs(g2-g)), grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
