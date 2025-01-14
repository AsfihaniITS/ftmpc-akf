function [S2,J12]=constraintDeltaInput(S,J1,usb,Np,n)

%define the matrix E1
S2 = zeros(2*Np,Np);
S2(1:4,1) = [S; -S];
dev = 1;
for i = 2: Np-1
    baw = i+dev;
    bak = i+dev+3;
    S2(baw:bak,i) = [S; -S];
    dev = dev + 1;
end
baw = 2*Np-1;
S2(baw:2*Np,Np) = [S];

%define the matrix F1+Ub in FTMPC-AKF
up=[usb(n-1);-usb(n-1)];
J12=[J1+up];
for i=2:Np
J12=[J12;J1];
end