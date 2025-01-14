function [F,Phi,epsil]=matrixFandPhiandEpsil(A,B,C,Np,phi,n)
%Creating matrix F
h(1,:)=C;
F(1,:)=C*A;
for kk=2:Np
h(kk,:)=h(kk-1,:)*A;
F(kk,:)= F(kk-1,:)*A;
end

%Creating matrix Phi
v=h*B;
Phi=zeros(Np,Np); %declare the dimension of Phi
Phi(:,1)=v; % first column of Phi
for i=2:Np
Phi(:,i)=[zeros(i-1,1);v(1:Np-i+1,1)]; %Toeplitz matrix
end

%Creating matrix E
e=h*phi(:,:,n-1);
epsil=zeros(Np,Np); %declare the dimension of Phi
epsil(:,1)=e; % first column
for i=2:Np
epsil(:,i)=[zeros(i-1,1);e(1:Np-i+1,1)]; %Toeplitz matrix
end