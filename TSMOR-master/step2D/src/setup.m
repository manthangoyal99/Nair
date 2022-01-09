function [Nx,Ny,N,Lx,Ly,x,y,dx,dy,x_grid,y_grid,m_bct,m,m_inv,Dx,Dy,Dxx,Dyy,gamma,nu,options_lspg,options_transport] = setup()

%Grid setup (specific to forward facing step geometry)---------------------
Lx=3;
Ly=1;
dx=0.0025;
dy=0.0025;
Nx=Lx/dx+1;
Ny=Ly/dy+1;
N=Nx*Ny;
gamma=1.4;
nu=0.003;

x=0:dx:Lx;
y=Ly:-dy:0;
[x_grid,y_grid]=meshgrid(x,y);
xs=find(x==0.6);    ys=find(y==0.2);    %Identifying grid location of corner of the step

%Derivative operators
Dx = sparse(fd_normal(Nx,3,x,1));
Dy = sparse(fd_normal(Ny,3,y,1));
%Making the first derivative first order at the outflow boundary
Dx(Nx,Nx-2)=0;
Dx(Nx,Nx-1)=-1/dx;
Dx(Nx,Nx)=1/dx;
Dxx = sparse(fd_normal(Nx,3,x,2));
Dyy = sparse(fd_normal(Ny,3,y,2));

%Levelset for 1st order outflow
for i=1:Nx
    for j=1:Ny
    if x(i)<=0.6
        levelset(i,j)=dy-abs(y(j));
    end
    if x(i)==0.6
        if y(j)<=0.2 && y(j)>0
            levelset(i,j)=y(j)-abs(y(j));
        end
    end
    if x(i)>0.6
        levelset(i,j)=0.2-abs(y(j));
    end
    if y(j)>=1-dy
    levelset(i,j)=abs(y(j))-(1-dy);        %For upper reflecting boundary
    end
    end
end  

m=(levelset>=0);        %m=1: boundary locations; m=0: wherever fluid exists including inlet
m_inv=(levelset<0);
m_bct=find(levelset'==0);
%m_interp=find(levelset<=0);
%m_ghost=find(levelset>0);
%x_gridt=x_grid';
%y_gridt=y_grid';
%--------------------------------------------------------------------------


%Optimization options:-----------------------------------------------------
%options1: Online stage -- performing LSPG
%options2: Offline stage -- computing transport field
    
options_lspg = optimoptions(...
        'fminunc',...
        'Algorithm',    'quasi-newton',...
        'HessUpdate',   'dfp',...
        'Display',      'iter',...
        'MaxFunEvals',  100000,...
        'MaxIter',      20,...
        'FinDiffType',  'central',...
        'TolFun',       1e-12,...
        'FiniteDifferenceStepSize', 1e-12);
    
options_transport = optimoptions(...
        'fmincon',...
        'Algorithm',    'sqp',...
        'Display',      'iter',...
        'MaxFunEvals',  100000,...
        'MaxIter',      50,...
        'FinDiffType',  'central',...
        'TolFun',       1e-12,...
        'ConstraintTolerance',  1e-12,...
        'FiniteDifferenceStepSize', 1e-12);
%-------------------------------------------------------------------------- 