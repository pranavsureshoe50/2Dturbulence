function turbulence2D
%simulates 2d turbulence in a periodic box
clear; clc;

Lx = 2*pi; %domain length in x-direction
Ly = 2*pi; %domain length in y-direction

Nx = 128; %number of grid points in x-direction
Ny = 128; %number of grid points in y-direction

Re = 5e4;  %Reynods number of the flow

x = (0:Nx-1)*Lx/Nx;
y = (0:Ny-1)*Ly/Ny;
[X, Y] = meshgrid(x, y);

k = [0:Nx/2-1, -Nx/2:-1]*(2*pi)/Lx;  %wavenumbers corresponding to x
l = [0:Ny/2-1, -Ny/2:-1]*(2*pi)/Ly;  %wavenumbers corresponding to y
[K, L] = meshgrid(k, l);

mu = 0.0;    %mean
sigma = 1.0; %standard deviation

%stream function in x-direction
psi = normrnd(mu, sigma, [Ny, Nx]);
psi_hat = fft2(psi);

u = ifft2(1i*L.*psi_hat);
sigma = real(sqrt(sum(sum(u.^2)))/(Nx*Ny));

psi_hat = psi_hat/sigma;
delt = .025; %initial time step

solve_rk3w_theta(psi_hat, Re, delt, K, L, X, Y, Lx, Ly, Nx, Ny);
%surf(X, Y, psi);
end

function solve_rk3w_theta(psi_hat, Re, delt, K, L, X, Y, Lx, Ly,Nx, Ny)
                                                                 
%RK3W-\theta parameters
alpha1 = 29/96; alpha2 = -3/40; alpha3 = 1/6;
beta1 = 37/160; beta2 = 5/24;   beta3 = 1/6;
gamma1 = 8/15;  gamma2 = 5/12;  gamma3 = 3/4;
zeta1 = -17/60; zeta2 = -5/12;

%defining del2 = K^2+L^2
del2 = K.*K+L.*L; 

Px = int32(Nx/2);   Py = int32(Ny/2);   %zero padding in x-andy-directions
%Mx = int32(3*Nx/2); My = int32(3*Ny/2);

%some indices used for zero padding
sx1 = 1; ex1 = Nx/2; sx2 = Nx/2+Px+1; ex2 = Nx+Px;
sy1 = 1; ey1 = Ny/2; sy2 = Ny/2+Py+1; ey2 = Nx+Py;

%some mzero-matrices used for zero paddin
C1 = zeros(Ny/2, Px)+1i*zeros(Ny/2, Px); 
C2 = zeros(Py, Nx+Px)+1i*zeros(Py, Nx+Px); 

omg = normrnd(0, 1, [Ny, Nx]);
omg_hat = fft2(omg);
u_hat = 1i*L.*omg_hat./del2; u_hat(1, 1) = 0.0+1i*0.0;

v_hat = -1i*K.*omg_hat./del2; v_hat(1, 1) = 0.0+1i*0.0;
enr0 = 0.5*sum(sum(u_hat.*conj(u_hat)+v_hat.*conj(v_hat)))/(Nx*Ny);

% sigma = sqrt(2*enr0);
% omg_hat = omg_hat/sigma;

enr = enr0;
fprintf('total energy at t=0 is %f\n', enr0);

fp = fopen('transient.dat', 'w');  %writing enry and enstrophy at each time
time = 0.0; n = 0;
while(time < 500)
    time = time+delt;
    fprintf('time: %f, delt: %12.8f, enr frac= %f\n', time, delt, enr/enr0);
    %**********************************************************************
    %************************FIRST SUBSTEP*********************************
    %**********************************************************************
    %vorticity field
    
    temp = (1i*K).*omg_hat;  
    temp_extd = [temp(1:Ny/2, 1:Nx/2), C1, temp(1:Ny/2, Nx/2+1:Nx);... 
                                          C2;...
                 temp(Ny/2+1:Ny, 1:Nx/2), C1, temp(Ny/2+1:Ny, Nx/2+1:Nx)];
    omg_conv_extd = ifft2(temp_extd);
    
    temp = (1i*L).*omg_hat./del2; temp(1, 1) = 0.0+1i*0.0;
    temp_extd = [temp(1:Ny/2, 1:Nx/2), C1, temp(1:Ny/2, Nx/2+1:Nx);... 
                                          C2;...
                 temp(Ny/2+1:Ny, 1:Nx/2), C1, temp(Ny/2+1:Ny, Nx/2+1:Nx)];
    psi_conv_extd = ifft2(temp_extd);
    conv_extd = -omg_conv_extd.*psi_conv_extd;
    
    temp = (1i*L).*omg_hat;
    temp_extd = [temp(1:Ny/2, 1:Nx/2), C1, temp(1:Ny/2, Nx/2+1:Nx);... 
                                          C2;...
                 temp(Ny/2+1:Ny, 1:Nx/2), C1, temp(Ny/2+1:Ny, Nx/2+1:Nx)];
    omg_conv_extd = ifft2(temp_extd);
    
    temp = (1i*K).*omg_hat./del2; temp(1, 1) = 0.0+1i*0.0;
    temp_extd = [temp(1:Ny/2, 1:Nx/2), C1, temp(1:Ny/2, Nx/2+1:Nx);... 
                                          C2;...
                 temp(Ny/2+1:Ny, 1:Nx/2), C1, temp(Ny/2+1:Ny, Nx/2+1:Nx)];
    psi_conv_extd = ifft2(temp_extd);
    conv_extd = conv_extd + omg_conv_extd.*psi_conv_extd;
    
    conv_hat_extd = fft2(conv_extd);
    NL1 = [conv_hat_extd(sy1:ey1, sx1:ex1), conv_hat_extd(sy1:ey1, sx2:ex2); ...
           conv_hat_extd(sy2:ey2, sx1:ex1), conv_hat_extd(sy2:ey2, sx2:ex2)];
    
    omg_hat = (omg_hat+delt*(-alpha1*(del2.*omg_hat)/Re+...
              (gamma1*NL1)))./(1.0+delt*beta1*del2/Re);
    %**********************************************************************
    %************************SECOND SUBSTEP********************************
    %**********************************************************************
    temp = (1i*K).*omg_hat;
    temp_extd = [temp(1:Ny/2, 1:Nx/2), C1, temp(1:Ny/2, Nx/2+1:Nx);... 
                                          C2;...
                 temp(Ny/2+1:Ny, 1:Nx/2), C1, temp(Ny/2+1:Ny, Nx/2+1:Nx)];
    omg_conv_extd = ifft2(temp_extd);
    
    temp = (1i*L).*omg_hat./del2; temp(1, 1) = 0.0;
    temp_extd = [temp(1:Ny/2, 1:Nx/2), C1, temp(1:Ny/2, Nx/2+1:Nx);... 
                                          C2;...
                 temp(Ny/2+1:Ny, 1:Nx/2), C1, temp(Ny/2+1:Ny, Nx/2+1:Nx)];
    psi_conv_extd = ifft2(temp_extd);
    conv_extd = -omg_conv_extd.*psi_conv_extd;
    
    temp = (1i*L).*omg_hat;
    temp_extd = [temp(1:Ny/2, 1:Nx/2), C1, temp(1:Ny/2, Nx/2+1:Nx);... 
                                          C2;...
                 temp(Ny/2+1:Ny, 1:Nx/2), C1, temp(Ny/2+1:Ny, Nx/2+1:Nx)];
    omg_conv_extd = ifft2(temp_extd);
    
    temp = (1i*K).*omg_hat./del2; temp(1, 1) = 0.0;
    temp_extd = [temp(1:Ny/2, 1:Nx/2), C1, temp(1:Ny/2, Nx/2+1:Nx);... 
                                          C2;...
                 temp(Ny/2+1:Ny, 1:Nx/2), C1, temp(Ny/2+1:Ny, Nx/2+1:Nx)];
    psi_conv_extd = ifft2(temp_extd);
    conv_extd = conv_extd + omg_conv_extd.*psi_conv_extd;
    
    conv_hat_extd = fft2(conv_extd);
    NL2 = [conv_hat_extd(sy1:ey1, sx1:ex1), conv_hat_extd(sy1:ey1, sx2:ex2); ...
           conv_hat_extd(sy2:ey2, sx1:ex1), conv_hat_extd(sy2:ey2, sx2:ex2)];
    
    omg_hat = (omg_hat+delt*(-alpha2*(del2.*omg_hat)/Re+...
              (gamma2*NL2+zeta1*NL1)))./(1.0+delt*beta2*del2/Re);
    
    NL1 = NL2;  %saving for the next step
    
    %**********************************************************************
    %**************************THIRD STEP**********************************
    %**********************************************************************
    
    temp = (1i*K).*omg_hat;
    temp_extd = [temp(1:Ny/2, 1:Nx/2), C1, temp(1:Ny/2, Nx/2+1:Nx);... 
                                          C2;...
                 temp(Ny/2+1:Ny, 1:Nx/2), C1, temp(Ny/2+1:Ny, Nx/2+1:Nx)];
    omg_conv_extd = ifft2(temp_extd);
    
    temp = (1i*L).*omg_hat./del2; temp(1, 1) = 0.0;
    temp_extd = [temp(1:Ny/2, 1:Nx/2), C1, temp(1:Ny/2, Nx/2+1:Nx);... 
                                          C2;...
                 temp(Ny/2+1:Ny, 1:Nx/2), C1, temp(Ny/2+1:Ny, Nx/2+1:Nx)];
    psi_conv_extd = ifft2(temp_extd);
    conv_extd = -omg_conv_extd.*psi_conv_extd;
    
    temp = (1i*L).*omg_hat;
    temp_extd = [temp(1:Ny/2, 1:Nx/2), C1, temp(1:Ny/2, Nx/2+1:Nx);... 
                                          C2;...
                 temp(Ny/2+1:Ny, 1:Nx/2), C1, temp(Ny/2+1:Ny, Nx/2+1:Nx)];
    omg_conv_extd = ifft2(temp_extd);
    
    temp = (1i*K).*omg_hat./del2; temp(1, 1) = 0.0;
    temp_extd = [temp(1:Ny/2, 1:Nx/2), C1, temp(1:Ny/2, Nx/2+1:Nx);... 
                                          C2;...
                 temp(Ny/2+1:Ny, 1:Nx/2), C1, temp(Ny/2+1:Ny, Nx/2+1:Nx)];
    psi_conv_extd = ifft2(temp_extd);
    conv_extd = conv_extd + omg_conv_extd.*psi_conv_extd;
    
    conv_hat_extd = fft2(conv_extd);
    NL2 = [conv_hat_extd(sy1:ey1, sx1:ex1), conv_hat_extd(sy1:ey1, sx2:ex2); ...
           conv_hat_extd(sy2:ey2, sx1:ex1), conv_hat_extd(sy2:ey2, sx2:ex2)];
    
    omg_hat = (omg_hat+delt*(-alpha3*(del2.*omg_hat)/Re+...
              (gamma3*NL2+zeta2*NL1)))./(1.0+delt*beta3*del2/Re);
    
    %**********************************************************************
    %**************************PLOTTING VORTICITY**************************
    %**********************************************************************
    n = n+1;  %number of time cycle
    
    %calculating total energy in the system, \sum(0.5*(u^2+v^2))
    u_hat = 1i*L.*omg_hat./del2; u_hat(1, 1) = 0.0+1i*0.0;

    v_hat = -1i*K.*omg_hat./del2; v_hat(1, 1) = 0.0+1i*0.0;
    enr = 0.5*sum(sum(u_hat.*conj(u_hat)+v_hat.*conj(v_hat)))/(Nx*Ny);
    
    %calculating total enstrophy
    ens = 0.5*(sum(sum(omg_hat.*conj(omg_hat))))/(Nx*Ny);
    
    %data output
    fprintf(fp, '%d    %14.8e    %14.8e    %14.8e\n', n, time, enr, ens);
    
    if(mod(n, 100)==0)
        omg = ifft2(omg_hat);
%         myFile = sprintf('omega_%d.dat', n);
%         dlmwrite(myFile, real(omg), 'delimiter', '\t', 'precision', '%14.8e');
        
        [hC, hC] = contourf(X, Y, real(omg), 50); axis equal;
        set(hC,'LineStyle','none'); colorbar; drawnow
    end
end
fclose(fp);
end