function [phi] = step_basis(C,q_exp,Nx,Ny,N,Lx,Ly,x_grid,y_grid,final,x,y,dy,m_bct,m,dMu,jerror)

iter=1;
for j=jerror%1:size(q_exp,2)
    
    cx=(C(1,1,j)+C(1,2,j)*sin(pi*x/Lx)+C(1,3,j)*sin(2*pi*x/Lx)+...
        +C(1,4,j)*sin(3*pi*x/Lx)+C(1,5,j)*sin(4*pi*x/Lx)...
        +C(1,6,j)*sin(5*pi*x/Lx)+C(1,7,j)*sin(6*pi*x/Lx)...
        +C(1,8,j)*sin(7*pi*x/Lx)+C(1,9,j)*sin(8*pi*x/Lx))*dMu(j)+...
        (C(1,10,j)+C(1,11,j)*sin(pi*x/Lx)+C(1,12,j)*sin(2*pi*x/Lx)+...
        +C(1,13,j)*sin(3*pi*x/Lx)+C(1,14,j)*sin(4*pi*x/Lx)...
        +C(1,15,j)*sin(5*pi*x/Lx)+C(1,16,j)*sin(6*pi*x/Lx)...
        +C(1,17,j)*sin(7*pi*x/Lx)+C(1,18,j)*sin(8*pi*x/Lx))*dMu(j)^2;

    x1=x-cx;
    
    cy=(C(2,1,j)+C(2,2,j)*sin(pi*y/Ly)+C(2,3,j)*sin(2*pi*y/Ly)+...
        +C(2,4,j)*sin(3*pi*y/Ly)+C(2,5,j)*sin(4*pi*y/Ly)...
        +C(2,6,j)*sin(5*pi*y/Ly)+C(2,7,j)*sin(6*pi*y/Ly)...
        +C(2,8,j)*sin(7*pi*y/Ly)+C(2,9,j)*sin(8*pi*y/Ly))*dMu(j)+...
        (C(2,10,j)+C(2,11,j)*sin(pi*y/Ly)+C(2,12,j)*sin(2*pi*y/Ly)+...
        +C(2,13,j)*sin(3*pi*y/Ly)+C(2,14,j)*sin(4*pi*y/Ly)...
        +C(2,15,j)*sin(5*pi*y/Ly)+C(2,16,j)*sin(6*pi*y/Ly)...
        +C(2,17,j)*sin(7*pi*y/Ly)+C(2,18,j)*sin(8*pi*y/Ly))*dMu(j)^2;
    
    y1=y-cy;
    
    for k=1:4
    u1=reshape(q_exp(N*(k-1)+1:N*(k-1)+N,j),Nx,Ny);
    u1t=u1';
    
    %Interpolation over whole domain using meshgrid%-----------------------
    if final==0  
    [x1_grid,y1_grid]=meshgrid(x1,y1);
    U{j}=(interp2(x1_grid,y1_grid,u1t,x_grid,y_grid,'spline'))';
    end
    
    % Interpolation with special treatment for boundary--------------------
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
            [~,jmin]=min(dist);
            u1tt(m_ghosttt(kk))=u1tt(m_bct(jmin));
        end
        U{j}=(interp2(x1_grid,y1_grid,u1tt,x_grid,y_grid,'spline'))';
    end
    
    phi(N*(k-1)+1:N*(k-1)+N,iter)=U{j}(:);
    end
    iter=iter+1;
end

end
