clc;clear; close all
% Code for 2D turbulence in a periodic box
format long

% Physical domain details
nx =1024;
ny =1024;
cfl = sqrt(2);
Re = 5.0e4; 
r = 1.0e-10;
lx = 8*pi; ly = 8*pi; dt = .01 ; end_t = 1500; cc=3;
plot_initial = 0;
% Transient data
t = zeros(end_t/dt , 1); dt_list = t ; energy_frac = t; ens = t;
dx = lx/nx; dy = ly/ny;
x = [0: nx-1]*lx/nx ;
y =  [0: ny-1]*ly/ny ;
[X,Y] = meshgrid(x,y);

kx = [0:nx/2-1, -nx/2:-1]*(2*pi)/lx;  %wavenumbers corresponding to x
ky = [0:ny/2-1, -ny/2:-1]*(2*pi)/ly;  %wavenumbers corresponding to y

[Kx,Ky] = meshgrid(kx,ky);
clear kx ky
% Matrices intialisation
omega = zeros(nx,ny);

% Intialiasing coefficients
g1 = 8/15; g2 = 5/12;g3 = 3/4; % gamma s
a1 = 29/96 ; a2 = -3/40 ; a3 = 1/6 ; % alpha s
b1 = 37/160 ; b2= 5/24; b3 = 1/6; % beta s
eta1 = -17/60 ;eta2 = -5/12; % eta s\

% Assign Initial condition
omega_new = intial_condition(omega , lx , ly , nx , ny, cc);
if plot_initial ==1
    m = mean(mean(abs(omega_new))); s = std(std(abs(omega_new))) ;
    contourf(y,x,real(omega_new),10,'EdgeColor','none'), colormap jet , title('Intial distribution'),    pause
end
omega_hat = fft2(omega_new);

% Wave numbers
K = Kx.^2 + Ky.^2 ;
linear = -K./Re ;
var=0.0+ 1i*0.0;
[E0,~] = energy(omega_hat,K,Kx,Ky) ;
% Stream function
psi_hat = omega_hat./K ;
psi_hat(1,1) =var ; psi_hat(nx/2 +1,nx/2 +1) =var ;psi_hat(nx/2 +1,1) = var;psi_hat(1,nx/2 +1) = var;
time = 0.0; n = 0; m = 0;

while time < end_t
    n=n+1 ; 
   
    % ----------------------------------------------------------------------------------------------
    % Post processing and time march
    [E_n,Ens_n] = energy(omega_hat,K,Kx,Ky) ;
    t(n) = time ; dt_list(n) = dt; energy_frac(n) = E_n/E0 ; ens(n) = Ens_n;
    % Plot the contour
    if time>2*m 
         omega = real(ifft2(omega_hat));
         
%         m1 = max(max((omega)))/1.2; %m2 = min(min((omega)))/1.2;
%         contourf(y,x,(omega),30,'EdgeColor','none'),title(['time:',num2str(time)]), colormap jet  ;
%         colorbar, caxis([-m1 m1])
%         
%         set(gcf,'Visible', 'off'); 
%         print(['image',num2str(floor(n/printevery)),'.png'],'-dpng','-r0')
%         hold on
        save(['data_',num2str(m)],'omega','time','E_n','Ens_n')
        m=m+1;
        
    end
    %-----------------------------------------------------------------------------------------------
    time = time+dt; % Next time pls
    
    %-----------------------------------------------------------------------------------------------
    % Substep1 : IMEX - RK3-CN\theta
    %-----------------------------------------------------------------------------------------------
    [Nu_hat,um] = nonlinearv5(psi_hat,omega_hat,nx,ny,Kx,Ky); % fft( (d\psi /dx * dw/dy) - (d\psi /dy * dw/dx))
    % dt update
    dt = 3*cfl*dx/(4*pi*abs(um)) ; 
    
    Lu_hat = linear .* omega_hat ;          % fft(-nu*k^2 \omega )
    
    omega_1 = omega_hat + dt*( (a1*Lu_hat) + (g1*Nu_hat) ) ;
    temp1 = 1.0 + r - ((b1*dt)*linear) ;
    omega_1 = omega_1./(temp1 ) ; 
    
    psi_hat = omega_1./(K +r) ; 
    psi_hat(1,1) =var ; psi_hat(nx/2 +1,nx/2 +1) =var ;psi_hat(nx/2 +1,1) = var;psi_hat(1,nx/2 +1) = var;

    %-----------------------------------------------------------------------------------------------
    % Substep2 : IMEX - RK3-CN\theta
    %-----------------------------------------------------------------------------------------------
    Nu_hat_prev = Nu_hat ;
    [Nu_hat,~] = nonlinearv5(psi_hat,omega_1,nx,ny,Kx,Ky); % fft( (d\psi /dx * dw/dy) - (d\psi /dy * dw/dx))
    Lu_hat = linear .* omega_1 ;                       % fft(-nu*k^2 \omega )
     
    omega_2 = omega_1 + dt*( (a2*Lu_hat) + (g2*Nu_hat) + (eta1*Nu_hat_prev) ) ;
    temp1 = 1.0 +r - ((b2*dt)*linear);
    omega_2 = omega_2./(temp1 );
 
    psi_hat = omega_2./(K+r) ; 
    psi_hat(1,1) =var ; psi_hat(nx/2 +1,nx/2 +1) =var ;psi_hat(nx/2 +1,1) = var;psi_hat(1,nx/2 +1) = var;
    
    %-----------------------------------------------------------------------------------------------
    % Substep3 : IMEX - RK3-CN\theta
    %-----------------------------------------------------------------------------------------------
    Nu_hat_prev = Nu_hat ;
    [Nu_hat,~] = nonlinearv5(psi_hat,omega_2,nx,ny,Kx,Ky); % fft( (d\psi /dx * dw/dy) - (d\psi /dy * dw/dx))
    Lu_hat = linear .* omega_2 ;                      % fft(-nu*k^2 \omega )
     
    omega_hat = omega_2 + dt*( (a3*Lu_hat) + (g3*Nu_hat) + (eta2*Nu_hat_prev) ) ;
    temp1 = 1.0 +r - ((b3*dt)*linear);
    omega_hat = omega_hat./(temp1 );
   
    psi_hat = omega_hat./(K+r) ; 
    psi_hat(1,1) =var ; psi_hat(nx/2 +1,nx/2 +1) =var ;psi_hat(nx/2 +1,1) = var;psi_hat(1,nx/2 +1) = var;
    if E_n>10*E0 || E_n ~= E_n, error('code blows') ; end
end
% Temporal data writing
t = t(1:n); dt_list = dt_list(1:n) ; energy_frac = energy_frac(1:n) ; ens = ens(1:n) ;
save('temporal','t','dt_list','energy_frac','ens')
