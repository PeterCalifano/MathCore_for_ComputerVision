function g = idft(G,f,t)

% inverse Fourier Transform using fft
fase = -2*pi*f(:)*t(1);
w = exp(1i*fase(:))*ones(1,size(G,2));
g = ifft(ifftshift(G.*conj(w),1));
if length(f)>length(t)
    g = g(1:length(t),:);
end