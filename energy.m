function [E,enstrpy] = energy(omega_hat,K,Kx,Ky)
    
    nx = length(omega_hat(:,1)) ;     ny = length(omega_hat(1,:)) ;
    
    u_hat = 1i*Ky.*omega_hat./K; u_hat(1, 1) = 0.0+1i*0.0;
    v_hat = -1i*Kx.*omega_hat./K; v_hat(1, 1) = 0.0+1i*0.0;

    E = 0; enstrpy = 0;
    for i=1:nx
        for j=1:ny
            E = E + 0.5*(u_hat(i,j).*conj(u_hat(i,j))+v_hat(i,j).*conj(v_hat(i,j)));
            enstrpy = enstrpy + 0.5*(omega_hat(i,j) .*conj(omega_hat(i,j)));
        end
    end
    E = E/(nx*ny) ; enstrpy = enstrpy/(nx*ny) ;
end