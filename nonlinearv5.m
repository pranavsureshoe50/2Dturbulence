function [Nu,u_max] = nonlinearv5(psi_hat,omega_hat,nx,ny,Kx,Ky)
    bx1 = nx/2 ; bx2 = nx/2 + (1.5*nx - nx) +1;
    by1 = ny/2 ; by2 = ny/2 + (1.5*ny - ny) +1;
    
    Nu = zeros(nx,ny);
    padded_mat(1:1.5*nx,1:1.5*ny) = 0.0 + 1i*0.0 ;
        
    % dw/dx calculation
    temp = (1i*Kx).*omega_hat;  
    padded_mat(1:bx1,1:by1) = temp(1:nx/2 , 1:ny/2) ;
    padded_mat(bx2:end,1:by1) = temp(nx/2+1:end , 1:ny/2) ;
    padded_mat(1:bx1,by2:end) = temp(1:nx/2,ny/2+1:end) ;
    padded_mat(bx2:end,by2:end) = temp(nx/2+1:end, ny/2+1:end);
          
    domegadx= ifft2(padded_mat);
    
    % dw/dy calculation
    temp = (1i*Ky).*omega_hat;  
    padded_mat(1:bx1,1:by1) = temp(1:nx/2 , 1:ny/2) ;
    padded_mat(bx2:end,1:by1) = temp(nx/2+1:end , 1:ny/2) ;
    padded_mat(1:bx1,by2:end) = temp(1:nx/2,ny/2+1:end) ;
    padded_mat(bx2:end,by2:end) = temp(nx/2+1:end, ny/2+1:end);         
    
    domegady= ifft2(padded_mat);
    
    % dpsi/dx calculation
    temp = (1i*Kx).*psi_hat;  
    padded_mat(1:bx1,1:by1) = temp(1:nx/2 , 1:ny/2) ;
    padded_mat(bx2:end,1:by1) = temp(nx/2+1:end , 1:ny/2) ;
    padded_mat(1:bx1,by2:end) = temp(1:nx/2,ny/2+1:end) ;
    padded_mat(bx2:end,by2:end) = temp(nx/2+1:end, ny/2+1:end); 
             
    dpsidx= ifft2(padded_mat);
    
    % dpsi/dx calculation
    temp = (1i*Ky).*psi_hat;  
    padded_mat(1:bx1,1:by1) = temp(1:nx/2 , 1:ny/2) ;
    padded_mat(bx2:end,1:by1) = temp(nx/2+1:end , 1:ny/2) ;
    padded_mat(1:bx1,by2:end) = temp(1:nx/2,ny/2+1:end) ;
    padded_mat(bx2:end,by2:end) = temp(nx/2+1:end, ny/2+1:end); 
             
    dpsidy= ifft2(padded_mat);
    u = abs(dpsidx) + abs(dpsidy) ;
    u_max = (max(max(u)));
    
    temp = -(dpsidy.*domegadx) + (dpsidx.*domegady) ;
    Nu_temp = fft2(temp);
    Nu(1:nx/2 , 1:ny/2) = Nu_temp(1:bx1,1:by1);
    Nu(nx/2+1:end , 1:ny/2) = Nu_temp(bx2:end,1:by1);
    Nu(1:nx/2,ny/2+1:end) = Nu_temp(1:bx1,by2:end) ;
    Nu(nx/2+1:end, ny/2+1:end) = Nu_temp(bx2:end,by2:end);

end