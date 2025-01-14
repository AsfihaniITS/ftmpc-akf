function [S1,T1]=constraintInput(T,S,Np)

%Define matriX S1
S1=zeros(2*Np,Np);
S1(1:2,1)=S;
for i=2:(Np)
S1(2*i-1:2*i,i)=S;
end

%Define matrix T1
T12 = T;
T1 = zeros(2*Np,1);
T1(1:2,:) = [T12];
for i = 2:Np
T1(2*i-1:2*i,1)= [T12]; 
end