function [R,q_new] = step_a(a,phi,Nx,Ny,N,Dx,Dy,Dxx,Dyy,gamma,nu,m,ii,Mupredict)

Q_new=phi*a;
for k=1:4
q_new{k}=reshape(Q_new(N*(k-1)+1:N*(k-1)+N,1),Nx,Ny);
end

% % %Enforcing inlet BC-------------------
% q_new{1}(1,:)=gamma;
% q_new{2}(1,:)=gamma*Mupredict(ii);
% q_new{3}(1,:)=0;
% q_new{4}(1,:)=1/(gamma-1)+0.5./q_new{1}(1,:).*(q_new{2}(1,:).^2+q_new{3}(1,:).^2);
% % %--------------------------------------   


u_new=q_new{2}./q_new{1};
v_new=q_new{3}./q_new{1};
E_new=q_new{4}./q_new{1};
p_new=(gamma-1)*q_new{1}.*(E_new-1/2*(u_new.^2 + v_new.^2));

sigma=100000; scale=1;
R1=Dx*q_new{2}+q_new{3}*Dy'-nu*Dxx*q_new{1}-nu*q_new{1}*Dyy'; R1(1,:)=scale*R1(1,:);
R2=Dx*(q_new{2}.*u_new+p_new)+(q_new{2}.*v_new)*Dy'-nu*Dxx*q_new{2}-nu*q_new{2}*Dyy'; R2(1,:)=scale*R2(1,:);
R3=Dx*(q_new{3}.*u_new)+(q_new{3}.*v_new+p_new)*Dy'-nu*Dxx*q_new{3}-nu*q_new{3}*Dyy'; R3(1,:)=scale*R3(1,:);
R4=Dx*((q_new{4}+p_new).*u_new)+((q_new{4}+p_new).*v_new)*Dy'-nu*Dxx*q_new{4}-nu*q_new{4}*Dyy'; R4(1,:)=scale*R4(1,:);
R1(m)=0;R2(m)=0;R3(m)=0;R4(m)=0;

R1(1,:)=sigma*(q_new{1}(1,:)-gamma);
R2(1,:)=sigma*(q_new{2}(1,:)-gamma*Mupredict(ii));
R3(1,:)=sigma*(q_new{3}(1,:)-0);
R4(1,:)=sigma*(q_new{4}(1,:)-(1/(gamma-1)+0.5*gamma*Mupredict(ii)^2));

% maxR=max([max(abs(R1(:))),max(abs(R2(:))),max(abs(R3(:))),max(abs(R4(:)))]);
% R=[R1(:)/max(abs(R1(:)));R2(:)/max(abs(R2(:)));R3(:)/max(abs(R3(:)));R4(:)/max(abs(R4(:)))];
% R=[R1(:)/q_new{1}(1,1);R2(:)/q_new{2}(1,1);R3(:)/q_new{2}(1,1);R4(:)/q_new{4}(1,1)];
R=[R1(:);R2(:);R3(:);R4(:)];
% R=[R1(:)*maxR/max(abs(R1(:)));R2(:)*maxR/max(abs(R2(:)));R3(:)*maxR/max(abs(R3(:)));R4(:)*maxR/max(abs(R4(:)))];

% R=sqrt(abs(R));%L1 for lsqnonlin
R=norm(R,1);%L1 for fminunc

end