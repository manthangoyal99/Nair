% Transported Snapshot Model Order Reduction (TSMOR) 
% for forward facing step
% Nirmal Nair and Maciej Balajewicz
% Based on "Transported snapshot model order reduction approach for parametric, 
% steady‐state fluid flows containing parameter‐dependent shocks." 
% International Journal for Numerical Methods in Engineering 117.12 (2019): 1234-1262.

clc; clear all;

%Setup grid, optimization options etc.
[Nx,Ny,N,Lx,Ly,x,y,dx,dy,x_grid,y_grid,m_bct,m,m_inv,Dx,Dy,Dxx,Dyy,gamma,nu,options_lspg,options_transport] = setup;

%Offline stage of TSMOR: load snapshots and compute transports=============

%Load snapshots corresponding to parameters Mu (training) and Mupredict (testing)
%The parameters are inlet Mach numbers
fprintf('Offline stage: Loading snapshots ...\n');
load('../snapshots/training_snapshots.mat') %This loads 'q_exp' which is used to cinstruct transport maps
load('../snapshots/testing_snapshots.mat') %This loads 'q_exp_predict' which is used to test the accuracy
Mu = 3.3:0.15:3.9;                                  %5 sampled snapshots for training
Mupredict = 3.3:0.05:3.9; Mupredict(1:3:end) = [];  %8 snapshots for prediction

%Construct transport field ------------------------------------------------
load_saved_transport = 1;   %Set it to 1 to load saved transport fields; 
                            %Set it to 0 to evaluate the tranports
nfc = 9;                    %Number of Fourier coefficients in each spatial dimension.                
                            %Choose from 1 to 9 (9 is used in the paper)
boundary_interpolation = 1; %Set it to 1 for special treatment at boundary while transporting the snapshot (recommended)
                            %Set it to 0 for regular interpolation over the entire domain
dd = Mu(2)-Mu(1);           %normalization constant
                          
if load_saved_transport==1
    fprintf('Offline stage: Loading transport field ...\n');
    load(['../transport_fields/fc' num2str(nfc) '.mat']) %Loads Fourier coeffs in 'C'
                    
elseif load_saved_transport==0
    C=zeros(2,18,size(q_exp,2)); %array of Fourier coefficients for each training snapshot and size(q_exp,2) gives the number of snapshots and 2 = (x,y) and 18 = (Np*Nq) 
    
    %Indices of 2 nearest neighboring snapshpots for each training snapshot
    %Manually coded for simplicity
    neighbor{1}=[2,3]; neighbor{2}=[1,3]; neighbor{3}=[2,4]; neighbor{4}=[2,3]; neighbor{5}=[3,4]; 
    
    %Loop over each sampled snapshot for evaluating transport fields
    for i=1:size(q_exp,2)
        fprintf('Offline stage: Constructing transport field for sampled snapshot #: %i\n',i);
        cc = C(:,:,i);
        ng = neighbor{i};
        [Error] = @(cc)transport(cc,i,q_exp,Nx,Ny,N,Lx,Ly,x,y,dx,dy,x_grid,y_grid,m_bct,m,boundary_interpolation,Mu,ng,dd);
        constraint = @(cc)transport_constraint(cc,i,Lx,Ly,x,y,Mu,ng,dd);
        cc = fmincon(Error,cc,[],[],[],[],[],[],constraint,options_transport);
        C(:,:,i)=cc; 
    end
end
%End of offline stage======================================================

%Online stage==============================================================
for ii=1:length(Mupredict)       %Loop over predictions
    fprintf('Online stage: Prediction at unsampled snapshot #: %i out of 8\n',ii);
    
    %Identifying indices of 2 nearest neighboring snapshpots
    amu = find(Mu<Mupredict(ii),1,'last'); bmu=amu+1;
    ng = [amu,bmu];
    
    %Initial guess of generzlized coordinates for LSPG based on the
    %paramettrical distance of sampled snapshot from the unsampled one
    dMu = zeros(1,size(q_exp,2)); w = zeros(1,size(q_exp,2));
    dMu(ng) = (Mupredict(ii)-Mu(ng))/dd; DMu = abs(dMu);
    w(ng) = (1./DMu(ng))/sum(1./DMu(ng));
    
    %Basis calculation by transporting the snapshot------------------------
    phi = basis(C,q_exp,Nx,Ny,N,Lx,Ly,x_grid,y_grid,boundary_interpolation,x,y,dy,m_bct,m,dMu,ng);

    %LSPG------------------------------------------------------------------
    a = w(ng)'; %setting initial guess
    [R] = @(a)residual(a,phi,Nx,Ny,N,Dx,Dy,Dxx,Dyy,gamma,nu,m,ii,Mupredict); %residual evaluation
    a = fminunc(R,a,options_lspg); %minimize residual
    [~,q_new] = residual(a,phi,Nx,Ny,N,Dx,Dy,Dxx,Dyy,gamma,nu,m,ii,Mupredict); %compute final solution in 'q_new'
    
    %Computing error---------------------------------------------------
    MM=[m_inv(:);m_inv(:);m_inv(:);m_inv(:)];
    for k=1:4 %loop over state variables (\rho, \rho*u, \rho*v, \rho*E)
        q_predict_TSMOR(N*(k-1)+1:N*(k-1)+N,ii)=q_new{k}(:);
        q_exp_predict_samp=q_exp_predict(N*(k-1)+1:N*(k-1)+N,ii); q_exp_predict_samp(m)=0;
        Error_TSMOR{k,ii}=abs(q_exp_predict(N*(k-1)+1:N*(k-1)+N,ii)-q_predict_TSMOR(N*(k-1)+1:N*(k-1)+N,ii));
        Error_TSMOR{k,ii}(m)=0;
        error_TSMOR(k,ii)=norm(Error_TSMOR{k,ii},'fro')./norm(q_exp_predict_samp,'fro')*100;
    end
    error_TSMOR_weighted(ii)=norm(q_exp_predict(MM,ii)-q_predict_TSMOR(MM,ii),'fro')/norm(q_exp_predict(MM,ii),'fro')*100;
    
    %Plotting--------------------------------------------------------------
    k=1; %1: \rho; 2: u; 3: v, 4: \rhoE
    figure
    subplot(2,1,1)
    q_interpp3=reshape(q_exp_predict(N*(k-1)+1:N*(k-1)+N,ii),[Nx,Ny]);
    if k==2 || k==3
    q_interpp3=reshape(q_exp_predict(N*(k-1)+1:N*(k-1)+N,ii)./q_exp_predict(N*(1-1)+1:N*(1-1)+N,ii),[Nx,Ny]);
    end
    q_interpp3(m)=NaN;
    pcolor(x,y,(q_interpp3)');
    axis equal;
    xlim([0 3]); ylim([0 1]);
    shading interp;
    colormap jet
    xlabel('x'); ylabel('y'); title('FOM')
    
    subplot(2,1,2)
    q_plot = reshape(q_predict_TSMOR(N*(k-1)+1:N*(k-1)+N,ii),[Nx,Ny]);
    if k==2 || k==3
    q_plot = reshape(q_predict_TSMOR(N*(k-1)+1:N*(k-1)+N,ii)./q_predict_TSMOR(N*(1-1)+1:N*(1-1)+N,ii),[Nx,Ny]);
    end
    q_plot(m)=NaN;
    pcolor(x,y,(q_plot)');
    axis equal;
    xlim([0 3]); ylim([0 1]);
    shading interp;
    colormap jet
    xlabel('x'); ylabel('y'); title('TSMOR')
    %----------------------------------------------------------------------
    
end

fprintf('Computed all predictions \n');
error_mean=mean(error_TSMOR_weighted);
error_max=max(error_TSMOR_weighted);
fprintf('Mean error: %f %% \n',error_mean);
fprintf('Max error : %f %% \n',error_max);




