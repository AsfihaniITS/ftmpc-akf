function [H,f]=quadraticfunction(F,Phi,xt,yref,Q,R,n,epsil2)
cul1=F*xt;
cul=cul1-yref(:,n-1)+epsil2(:,:,n-1);
H=2*(Phi'*Q*Phi+R);
f=2*Phi'*Q*cul;