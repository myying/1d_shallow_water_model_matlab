clear all; close all;

show_animation=1;      %if 1, show animated results after model run
animation_stride=1;    %if =n, skip every n frames
show_analytic=0;       %show analytical solution (linear advection eqn)
show_wind=1;           %show wind vector

%forcing term switches
nifcor=0;   %if 1, include Coriolis
nifwind=1;  %if 1, include other terms besides Coriolis in u,v eqns.
nifdif=0;   %if 1, include diffusion
nifad=1;    %if 1, include advection

%%%%%%%PARAMETERS
%numerical scheme switches
advection_scheme=1; %1:CTCS, 2:FTCS(unstable), 3: FTBS numerical schemes
diffusion_scheme=2; %1:CTCS (unstable), 2:FTCS
x_order=2; %order of accuracy for CS (2,4)
bc_type=1; %0:cyclic; 1:reflective; 2:radiative 
bc_order=2; %order of accuracy (1,2)

%computation grid
dx=100000;       %grid spacing (m)
nx=51;           %number of grid points
dt=100;          %time step (s)
nt=48*3600/dt;   %number of time steps
t(1:nt+1)=0:dt:nt*dt;
x(1:nx)=(0:nx-1)*dx/1000;
L=(nx-1)*dx/1000; %domain length

%initial condition
pi=4*atan(1.0);
C=0;       %advection velocity
hbar=8000;  %base level of height
amp=100;    %amplitude of hill
hc=15;      %center location of hill
hw=10;      %width of hill

%physics parameters
g=9.8;      %gravity
f=0.0001;   %Coriolis
kdif=0.1*dx*dx/(2*dt);  %diffuction coef

%%%%%%%%INITIALIZE
u(1:nt+1,1:nx)=NaN; v(1:nt+1,1:nx)=NaN; h(1:nt+1,1:nx)=NaN;
u(1,1:nx)=C;
v(1,1:nx)=0.0;
h(1,1:nx)=hbar;
ind=ceil(hc-hw/2):floor(hc+hw/2);
h(1,ind)=hbar+(amp/2)*(1+cos(2*pi*(ind-hc)/hw));

offset=@(j,d) mod(j-1+d,nx)+1;
ddx_CS2=@(a,j) (a(offset(j,1))-a(offset(j,-1)))/(2*dx);
ddx_FS1=@(a,j) (a(offset(j,1))-a(j))/dx;
ddx_BS1=@(a,j) (a(j)-a(offset(j,-1)))/dx;
ddx2_CS2=@(a,j) (a(offset(j,1))+a(offset(j,-1))-2*a(j))/(dx*dx);
d_h=@(u,v,h,j,ddx,dt) dt*(-nifad*u(j).*ddx(h,j)-h(j).*ddx(u,j));
d_u=@(u,v,h,j,ddx,dt) dt*nifwind*(-nifad*u(j).*ddx(u,j)-g*ddx(h,j)+nifcor*f*v(j));
d_v=@(u,v,h,j,ddx,dt) dt*nifwind*(-nifad*u(j).*ddx(v,j)-nifcor*f*u(j));
d_dif=@(a,j,ddx2,dt) dt*nifdif*kdif*ddx2(a,j);

%%%%%%%%%RUN MODEL
t_start=cputime;
for n=1:nt
  h_inc(1:nx)=0; u_inc(1:nx)=0; v_inc(1:nx)=0;
  
  if(advection_scheme==1) %CTCS
    nb=max(1,n-1); na=n+1;
    dtcalc=dt*(na-nb);
    ddx=ddx_CS2;
  elseif(advection_scheme==2) %FTCS
    nb=n; na=n+1;
    dtcalc=dt;
    ddx=ddx_CS2;
  elseif(advection_scheme==3) %FTBS
    nb=n; na=n+1;
    dtcalc=dt;
    ddx=ddx_BS1;
  end
  
  if(bc_type==0) %cyclic
    j=1:nx;
  else
    j=2:nx-1;
  end

  h(na,j)=h(nb,j)+d_h(u(n,:),v(n,:),h(n,:),j,ddx,dtcalc)+d_dif(h(nb,:),j,ddx2_CS2,dtcalc);
  u(na,j)=u(nb,j)+d_u(u(n,:),v(n,:),h(n,:),j,ddx,dtcalc)+d_dif(u(nb,:),j,ddx2_CS2,dtcalc);
  v(na,j)=v(nb,j)+d_v(u(n,:),v(n,:),h(n,:),j,ddx,dtcalc)+d_dif(v(nb,:),j,ddx2_CS2,dtcalc);
  
  if(bc_type~=0)
    h(na,1)=h(nb,1)+d_h(u(n,:),v(n,:),h(n,:),1,ddx_FS1,dtcalc);
    u(na,1)=u(nb,1)+d_u(u(n,:),v(n,:),h(n,:),1,ddx_FS1,dtcalc);
    v(na,1)=v(nb,1)+d_v(u(n,:),v(n,:),h(n,:),1,ddx_FS1,dtcalc);
    h(na,nx)=h(nb,nx)+d_h(u(n,:),v(n,:),h(n,:),nx,ddx_BS1,dtcalc);
    u(na,nx)=u(nb,nx)+d_u(u(n,:),v(n,:),h(n,:),nx,ddx_BS1,dtcalc);
    v(na,nx)=v(nb,nx)+d_v(u(n,:),v(n,:),h(n,:),nx,ddx_BS1,dtcalc);
  end
  
  if(bc_type==1)  %reflective bc
    u(n+1,1)=0; u(n+1,nx)=0;
    if(bc_order==1)
      h(n+1,1)=h(n+1,2); h(n+1,nx)=h(n+1,nx-1);
    elseif(bc_order==2)
      h(n+1,1)=(4*h(n+1,2)-h(n+1,3))/3;
      h(n+1,nx)=(4*h(n+1,nx-1)-h(n+1,nx-2))/3;
    end
  elseif(bc_type==2)  %radiative bc
    
  end
  
  if(any(h(n+1,:)>1e10)>0) 
    disp(['model exploded at t=' num2str(t(n)/3600) 'h (' num2str(n) ' steps)']);
    break; 
  end
end
t_end=cputime;
disp(['model integration takes ' num2str(t_end-t_start) ' s']);

%%%%%%%%%DIAGNOSTICS
% m=20; %trace max(h) across m grid points
% loc1=find(h(1,:)==max(h(1,:)),1); 
% loc2=find(h(floor(m*dx/C/dt),:)==max(h(floor(m*dx/C/dt),:)),1);
% Cd=(loc2-loc1)*dx/(floor(m*dx/C/dt)*dt); disp(['Cd=' num2str(Cd)]);
% maxh=max(h');  %time series of peak height amplitude

%%%%%%%%%PLOT RESULTS
figure('Position',[100 100 560 560]);
x0(1:nx)=NaN; h0(1:nx)=NaN;
if(show_animation==1)
  for n=1:animation_stride:nt+1
    if(any(h(n,:)>1e10) || any(isnan(h(n,:)))) break; end
    subplot(2,1,1)
    plot(x,h(n,:)); hold on
    if(show_analytic==1) %analytic solution for advection equation
      for k=1:nx
        xloc=mod(x(k)+C*t(n)/1000,L);
        jloc=xloc*1000/dx;
        x0(floor(jloc)+1)=xloc;
        h0(floor(jloc)+1)=h(1,k);
      end
      plot(x0,h0,':','color',[.5 .5 .5],'linewidth',2);
    end
    hold off;
    axis([min(x) max(x) hbar-100 hbar+120]);
    ylabel('height'); title(['t = ' num2str(t(n)/3600) ' h'])
    xlabel('x (km)')
    if(show_wind==1)
      subplot(2,1,2)
      quiver(x,zeros(1,nx),u(n,:)*100,v(n,:)*100,'autoscale','off');
      axis([min(x) max(x) -max(x)/5 max(x)/5]);
      xlabel('x (km)'); ylabel('y (km)');
    end
    pause(0.001);
  end
end

%%%%%Hovemoller diagram of height
% figure('Position',[660 100 560 560])
% contourf(x,t,h,'linestyle','none'); colorbar
% xlabel('x (km)'); ylabel('t (h)'); title('height')

