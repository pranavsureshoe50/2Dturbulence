function omega_new = intial_condition(omega , lx , ly , nx , ny,c)
%         %Use all variables in fourier domain

    x =  [0: nx-1]*lx/nx ;
    y =  [0: ny-1]*ly/ny ;
    [X,Y] = meshgrid(x,y);
    
    kx = [0:nx/2-1,  -nx/2:-1]*(2*pi)/lx;  %wavenumbers corresponding to x
    ky = [0:ny/2-1,  -ny/2:-1]*(2*pi)/ly;  %wavenumbers corresponding to y

    [Kx,Ky] = meshgrid(kx,ky);
    K = sqrt(Kx.^2 + Ky.^2) ;
    clear kx ky
% Initial condition 1
switch c
    case 1
    for i=1:nx
        for j=1:ny
            if omega(i,j) > 1,  omega(i,j) = 1.0/omega(i,j) ; end
            omega(i,j) = omega(i,j) * exp(1i*pi*rand(1,1)) ;
        end
    end
    omega_new = omega;

% Initial condition 2
    case 2
     e = 0.3; 
    for i=1:nx
         for j=1:ny
            omega(i,j) = exp (-((x(i)-pi+pi/5).^2+(y(j)-pi+pi/5).^2)/(0.3))-exp(-((x(i)-pi-pi/5).^2+(y(j)-pi+pi/5).^2)/(0.2))+ ...
            exp( -(( x(i)-pi-pi/5 ).^2+( y(j)-pi- pi/5 ).^2 ) / ( 0.4 ) ) ;
         end
    end
    Noise = random ('unif',-1,1,nx,ny) ;
    omega_new = omega + e.*Noise ; 

% Initial condition 3
    case 3
    mu = 0.0;    %mean
    sigma = 1; %standard deviation
    
    omega_new = 1.25*normrnd(mu, sigma, [ny, nx]);
    
% Initial condition 4
    case 4
        norm1 = mvnpdf([X(:) Y(:)],[0.3,0.5],[0.01,0.01]);
        norm1 = reshape(norm1,nx,ny);

        norm2 = mvnpdf([X(:) Y(:)],[0.7,0.5],[0.01,0.01]);
        norm2 = reshape(norm2,nx,ny);

        omega_new    = norm1 - norm2;
        
% Initial condition 5 - Random condition
    case 5
        a = -3;
        b = 3;
        omega_new = (b-a).*rand(nx,ny) +a; 
    case 6
        a = 0; r=1.0e-6;
        b = 2*pi;
        q = (b-a).*rand(nx,ny) + a;
        phi = (b-a).*rand(nx,ny) + a;
        E=K;
        for i=1:nx
           for j=1:ny 
               if K(i,j)<30
                   E(i,j) =0;
               else
                   E(i,j) =  K(i,j)^4 /((1+K(i,j)^2)^3);%(K(i,j)^4)*exp(-(K(i,j)^2));% ;
               end
           end
        end
        alpha = sqrt(E./(4*pi*(K.^2 +r))).*exp(1i*q).*cos(phi) ;%sqrt(E./(2*pi*(K+r))).*exp(1i*q).*cos(phi) ;
        omega = -1i*(alpha.*K);
        omega_new = 3*(10^4)*ifft2(omega);
       % omega = -1i*sqrt(E./(4*pi)).*exp(1i*q) .* cos(phi);
       % omega_new = (10^6)*ifft2(omega);
    case 7
        a = 0;
        b = 2*pi;
        q = (b-a).*rand(nx,ny) + a;
        sigma =20;
        %omega_new =omega ;
        for i=1:nx
            for j=1:ny
            omega(i,j) = sigma^2 * (2*lx/pi)* (1/ (1+ (1.339*lx*K(i,j))^2)^.5) ;
            omega(i,j) = omega(i,j) * exp(1i*q(i,j)) ; 
            end
        end
        omega_new = 100*ifft2(omega) ;
    
end   
end
