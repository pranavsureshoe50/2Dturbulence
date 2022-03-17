% Post_process data
clc;clear; close all

% Code for 2D turbulence in a periodic box
format long
%plot location
data_dir = 'data_plots_512x512/data/';
run_dir = 'image_';
%cd(data_dir)
% Physical domain details
nx =1024;
ny =1024;

lx = 8*pi; ly =8*pi; 
% Transient data

dx = lx/nx; dy = ly/ny;
x = [0: nx-1]*lx/nx ;   

y =  [0: ny-1]*ly/ny ;
[X,Y] = meshgrid(x,y);
load('data_0.mat');
E0 = E_n ;
t_scale = (lx/sqrt(2*E0));
 
st =269; n =436; % num of files
E_vec = zeros(n,1); Ens_vec = E_vec;time1 =E_vec;

%system('mkdir plots')

for (i=st:n)
        
        [E_n,Ens_n,time] = mapper(i,t_scale,lx,ly,x,y);
        
         set(gcf,'Visible', 'off'); 
         print([run_dir,num2str(i),'.png'],'-dpng','-r0')
         hold on
        time1(i+1) = time;
        E_vec(i+1) = E_n/E0 ;
        Ens_vec(i+1) = Ens_n ;
end

function [e1,e2,time] = mapper(i,t_scale,lx,ly,x,y)
        load(['data_',num2str(i),'.mat'])
        e1=E_n;e2=Ens_n;
         time = time/ t_scale;  
        m1 = max(max((omega*t_scale)))/1.2; %m2 = min(min((omega)))/1.2;
        contourf(y/ly,x/lx,(omega*t_scale),30,'EdgeColor','none'),title(['time (tu_{rms}/L):',num2str(time)]), colormap jet  ;
        colorbar, caxis([-m1 m1])
        
end
