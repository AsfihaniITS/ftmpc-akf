clc;
clear all;
close all;

Delta_t=0.25;
[A,B,C]=modelMPCDavidson(Delta_t);%A is the state transition matrix, B is the input matrix and C represents the measurement matrix

%==========Define prediction time(t), sampling time(ts), Np=Nc, Total time(TotTime)
Ts=1;         %sampling time
t=1000;
TotTime=t/Ts; 
time=0:TotTime-1;
Np=15;        %prediction horizon

Skf(:,:,1) = 0.1*eye(1);%initialization value for the matrix S(k) in AKF

teta1(1:1000)=0;
teta1(200:1000)=0.2;
teta(:,:)=[teta1];%initial value of actuator fault

%define the initial condition of the state
x1(1)=0; %initial value of sway rate
x2(1)=0; %initial value of yaw rate
x3(1)=deg2rad(0) ;%initial value of yaw angle
x131(1)=rad2deg(x3(1));
xlt=[x1(1);x2(1);x3(1)];
xt=xlt;
xa(:,1)=xt;
xr(:,1)=xt;


%Define weighting matrices Q and R for the output and input 
R=0.01 *eye(Np,Np);
Q=1*eye(Np,Np);

S=[1;-1]; %matrix S in FTMPC-AKF
T=[deg2rad(35);deg2rad(35)];%The constraint on the rudder angle, matrix T in FTMPC-AKF
J1=[deg2rad(5);deg2rad(5)];%The constraint on the rate of change of rudder angle, matrix V in FTMPC-AKF

state=3;
q1 = 10^(-10); %model error of sway rate
q2 = 10^(-10); %model error of yaw rate
q3 = 10^(-10); %model error of yaw angle
r = 10^(-10); %measurement error

pa(:,:,1) = [1 0 0;0 1 0;0 0 1]; %initial value of the error covariance matrix
z(:,1) = C*xr(:,1) + normrnd(0,sqrt(r),1,1);%output vector
%define references for the yaw angle
reff=linspace(1,360,TotTime);
for i=1:TotTime
yref(:,i)=deg2rad(reff(:,i))*ones(Np,1);
end
 
%===initial optimal solution  
usb(1)=deg2rad(0);
usbk(1)=rad2deg(usb(1));

lambda=0.9;%initial value of forgetting factor
phi(:,:,1)=-B*diag(usb(1));%matrix corresponding to fault
gamma(:,:,1)=[0;0;0];

tetatopi(:,1)=teta(:,1);%initial value of fault estimation
tetatopi2=zeros(Np,1);
tetatopi2(:,1)=teta(:,1);

epsil2(:,:,1)=zeros(Np,1);

%===================================================================
%================== simulation is started from here ================
%===================================================================

for n=2:TotTime

%===Creating Matrix F, Phi and E in FTMPC-AKF
[F,Phi,epsil]=matrixFandPhiandEpsil(A,B,C,Np,phi,n);
epsil2(:,:,n)=epsil*tetatopi2(:,n-1);

%creating objective function
[H,f]=quadraticfunction(F,Phi,xt,yref,Q,R,n,epsil2);

%define the matrix S1 and T1 in FTMPC-AKF (The constraint on the rudder angle)
[S1,T1]=constraintInput(T,S,Np);

%define the matrix E1 and F1+Ub in FTMPC-AKF (The constraint on the rate of change of rudder angle)
[S2,J12]=constraintDeltaInput(S,J1,usb,Np,n);

%=====Defines a Matrix to use in a quadratic programing
Aken=[S1;S2];
Bken=[T1;J12];

%=====initial optimal solution
m=zeros(Np,1);
options=optimset('largescale','off');
ut=quadprog(H,f,Aken,Bken,[],[],[],[],m,options);

usb(n)=ut(1,1);%we only take the first input of the sequences

%---Adaptive Kalman Filter---
    phi(:,:,n)=-B*diag(usb(n-1));%matrix corresponding to fault
    xr(:,n) = A*xr(:,n-1) + B*usb(n) +phi(:,:,n-1)*teta(:,n-1)+ eye(state)*[normrnd(0,sqrt(q1),1,1); normrnd(0,sqrt(q2),1,1); normrnd(0,sqrt(q3),1,1)];%system model
    z(:,n) = C*xr(:,n) + normrnd(0,sqrt(r),1,1);      %measurement model
    
    pf(:,:,n) = A*pa(:,:,n-1)*A' + q1*eye(state);
    sigma(:,:,n)=C*pf(:,:,n)*C'+ eye(1)*r;
    k(:,:,n) = pf(:,:,n)*C'*inv(sigma(:,:,n));%Kalman Gain
    pa(:,:,n) = (eye(3) - k(:,:,n)*C)*pf(:,:,n);%error covariance matrix
    
    gamma(:,:,n)=(eye(state) - k(:,:,n)*C)*A*gamma(:,:,n-1)+(eye(state) - k(:,:,n)*C)*phi(:,:,n-1);
    omega(:,:,n)=C*A*gamma(:,:,n-1)+C*phi(:,:,n-1);
    LAMBDA(:,:,n)=inv(lambda*sigma(:,:,n)+omega(:,:,n)*Skf(:,:,n-1)*omega(:,:,n)');
    tho(:,:,n)=Skf(:,:,n-1)*omega(:,:,n)'*LAMBDA(:,:,n);%the fault estimation gain matrix
    Skf(:,:,n)=(1/lambda)*Skf(:,:,n-1)-(1/lambda)*Skf(:,:,n-1)*omega(:,:,n)'*LAMBDA(:,:,n)*omega(:,:,n)*Skf(:,:,n-1);
    
    ytopi(:,n)=z(:,n)-C*(A*xa(:,n-1)+B*usb(n)+phi(:,:,n-1)*tetatopi(:,n-1));%measurement error
    
    tetatopi(:,n)=tetatopi(:,n-1)+tho(:,:,n)*ytopi(:,n);%fault estimation
    xa(:,n)=A*xa(:,n-1)+B*usb(n)+phi(:,:,n-1)*tetatopi(:,n-1)+k(:,:,n)*ytopi(:,n)+gamma(:,:,n)*(tetatopi(:,n)-tetatopi(:,n-1));%state estimation

%=====save the state estimation
x1(n)=xa(1,n);
x2(n)=xa(2,n);
x3(n)=xa(3,n);
xt=[x1(n);x2(n);x3(n)];
 
usbk(n)=rad2deg(usb(n)); %save the optimal every iteration
x131(n)=rad2deg(x3(n));

squareError(1,n)=[(x3(n)-yref(1,n))^2];%calculating error of yaw angle and reference angle
squareError1(1,n)=[(tetatopi(:,n)-teta(:,n))^2];%calculating error of fault estimation and actuator fault
tetatopi2(:,n)=tetatopi(:,n);
end

for n=2:TotTime
    beda_U(n)=abs(usbk(n)-usbk(n-1));%calculating change of rudder angle
end

%Calculating Root Mean Square Error
meanSquareError=mean(squareError,2);
rootMeanSquareError = sqrt(meanSquareError);

meanSquareError1=mean(squareError1,2);
rootMeanSquareError1 = sqrt(meanSquareError1);

%=====only for plotting trajectory
L=101.07;%LPP-Length (DWL)- meters
Lx=L*cosd(x131);
Ly=L*sind(x131);
R=101.07;%LPP-Length (DWL)- meters
Rx=R*cosd(reff);
Ry=R*sind(reff);

clf;
figure(1)
plot(time,reff,'--r',time,x131,'k','linewidth',2),xlabel('Time (s)'),ylabel('\psi (deg)')
hleg1 = legend('reference','FTMPC-AKF');
set(gcf,'color','w')
 
figure(2)
plot(time,usbk,'k','linewidth',2),xlabel('Time(s)'),ylabel('u (deg)')
set(gcf,'color','w')

figure(3)
plot(time,x1,'k','linewidth',2),xlabel('Time (s)'),ylabel('Sway rate (m/s)')
set(gcf,'color','w')

figure(4)
plot(time,x2,'k','linewidth',2),xlabel('Time (s)'),ylabel('r (rad/s)')
set(gcf,'color','w')

figure(5)
plot(time,tetatopi(1,:),'k',time,teta(1,:),'--r','linewidth',2)
grid on;
xlabel('Time (s)')
ylabel('\theta')
legend('fault estimation','actual')

figure(6)
plot(Rx,Ry,'--r',Lx,Ly,'k','linewidth',2)
grid on;
xlabel('x-axis')
ylabel('y-axis')
legend('reference','FTMPC-AKF')
