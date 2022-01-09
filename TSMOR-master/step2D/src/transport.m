function [error] = step_map(cc,i,q_exp,Nx,Ny,N,Lx,Ly,x,y,dx,dy,x_grid,y_grid,m_bct,m,final,Mu,ng,dd)

u1=reshape(q_exp(1:N,i),Nx,Ny);
u1t=u1';

for j=ng%1:size(q_exp,2)
    dMu = (Mu(j)-Mu(i))/dd;
    
    cx=(cc(1,1)+cc(1,2)*sin(pi*x/Lx)+cc(1,3)*sin(2*pi*x/Lx)+...
        +cc(1,4)*sin(3*pi*x/Lx)+cc(1,5)*sin(4*pi*x/Lx)...
        +cc(1,6)*sin(5*pi*x/Lx)+cc(1,7)*sin(6*pi*x/Lx)...
        +cc(1,8)*sin(7*pi*x/Lx)+cc(1,9)*sin(8*pi*x/Lx))*dMu+...
        (cc(1,10)+cc(1,11)*sin(pi*x/Lx)+cc(1,12)*sin(2*pi*x/Lx)+...
        +cc(1,13)*sin(3*pi*x/Lx)+cc(1,14)*sin(4*pi*x/Lx)...
        +cc(1,15)*sin(5*pi*x/Lx)+cc(1,16)*sin(6*pi*x/Lx)...
        +cc(1,17)*sin(7*pi*x/Lx)+cc(1,18)*sin(8*pi*x/Lx))*dMu^2;

    x1=x-cx; %transport in x
    
    cy=(cc(2,1)+cc(2,2)*sin(pi*y/Ly)+cc(2,3)*sin(2*pi*y/Ly)+...
        +cc(2,4)*sin(3*pi*y/Ly)+cc(2,5)*sin(4*pi*y/Ly)...
        +cc(2,6)*sin(5*pi*y/Ly)+cc(2,7)*sin(6*pi*y/Ly)...
        +cc(2,8)*sin(7*pi*y/Ly)+cc(2,9)*sin(8*pi*y/Ly))*dMu+...
        (cc(2,10)+cc(2,11)*sin(pi*y/Ly)+cc(2,12)*sin(2*pi*y/Ly)+...
        +cc(2,13)*sin(3*pi*y/Ly)+cc(2,14)*sin(4*pi*y/Ly)...
        +cc(2,15)*sin(5*pi*y/Ly)+cc(2,16)*sin(6*pi*y/Ly)...
        +cc(2,17)*sin(7*pi*y/Ly)+cc(2,18)*sin(8*pi*y/Ly))*dMu^2;
    
    y1=y-cy; %transport in x
    
    %Interpolation over whole domain using meshgrid-----------------------
    if final==0  
    [x1_grid,y1_grid]=meshgrid(x1,y1);
    U{j}=(interp2(x1_grid,y1_grid,u1t,x_grid,y_grid,'spline'))'; %Nx*Ny
    U{j}(m)=0;% Most probably makes all the boundary point velocities 0 
    end
    
    %Interpolation with special treatment at boundary---------------------
    if final==1
        u1tt=u1t;
        % Levelset for 1st order outflow----------------------
        levelset=NaN*ones(Nx,Ny);
        for kk=1:Nx
            for jj=1:Ny
                if x1(kk)<=0.6
                    levelset(kk,jj)=dy-y1(jj);
                end
                if x1(kk)==0.6
                    if y1(jj)<=0.2 && y1(jj)>=dy
                        levelset(kk,jj)=y1(jj)-y1(jj);
                    end
                end
                if x1(kk)>0.6
                    levelset(kk,jj)=0.2-y1(jj);
                end
                if y1(jj)>=1-dy
                    levelset(kk,jj)=y1(jj)-(1-dy);        %For upper reflecting boundary
                end
            end
        end
        %------------------------------------------
        m_out=(levelset'<=0);
        m_intersect=m'+m_out;
        m_ghosttt=find(m_intersect==2);
        [x1_grid,y1_grid]=meshgrid(x1,y1);
        for kk=1:length(m_ghosttt);
            dist=sqrt((x1_grid(m_ghosttt(kk))-x_grid(m_bct(:))).^2+(y1_grid(m_ghosttt(kk))-y_grid(m_bct(:))).^2)...
              +sqrt((x_grid(m_ghosttt(kk))-x_grid(m_bct(:))).^2+(y_grid(m_ghosttt(kk))-y_grid(m_bct(:))).^2);
            [~,jmin]=min(dist);%jmin is the index of the min value
            u1tt(m_ghosttt(kk))=u1tt(m_bct(jmin));
        end
        U{j}=(interp2(x1_grid,y1_grid,u1tt,x_grid,y_grid,'spline'))';
        U{j}(m)=0;
    end
    
    
    q_exp(m,j)=0;
    error(j)=norm(U{j}(:)-q_exp(1:N,j),'fro'); %error between the transported snapshot and neighboring snapshot
    
end

error=norm(error,'fro');
end
