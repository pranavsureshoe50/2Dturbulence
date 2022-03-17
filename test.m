L =2*pi;
nx = 64;
x = (0:nx-1)*L/nx ;

kx = [0:nx/2-1,  -nx/2:-1]*(2*pi)/L;  %wavenumbers corresponding to x

u = cos(2*x);
%plot(x,u)
dudx = 1i*kx.*fft(u) ;
du = ifft(dudx);
plot(x,du)