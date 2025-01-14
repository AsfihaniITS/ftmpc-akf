function [A,B,C]=modelMPCDavidson(Delta_t)
%Using the Strip Theory approach
L=101.07;       %LPP-Length (DWL)- meters
B=14;           %Width - meters
T=3.7;          %Design Draught (DWL) - meters
m=2423000;      %Displacement - kg
U=15.4;         %Ship Service Speed - meters/second
CB=0.65;        %Block Coefficient
xG=5.25;        %Center of gravity
rho=1024;       %Sea water mass density - kg/m^3
Adelta=5.7224;  %rudder width
r=0.156*L;      %radius of gyration

%Hydrodynamic Coefficient
Yvdot=-1*((1+0.16*CB*(B/T)-5.1*(B/L)^2)*pi*(T/L)^2);
Yrdot=-1*((0.67*(B/L)-0.0033*(B/T)^2)*pi*(T/L)^2);
Nvdot=-1*((1.1*(B/L)-0.041*(B/T))*pi*(T/L)^2);
Nrdot=-1*(((1/12)+0.017*(CB*B/T)-0.33*(B/L))*pi*(T/L)^2);
Yv=-1*((1+0.4*(CB*B/T))*pi*(T/L)^2);
Yr=-1*((-0.5+2.2*(B/L)-0.08*(B/T))*pi*(T/L)^2);
Nv=-1*((0.5+2.4*(T/L))*pi*(T/L)^2);
Nr=-1*((0.25+0.039*(B/T)-0.56*(B/L))*pi*(T/L)^2);
Ydelta=rho*pi*Adelta/(4*L*T);
Ndelta=-0.5*Ydelta;
Ir=(m*r^2)/(0.5*rho*L^5);
Iz=(m*(xG^2))/(0.5*rho*L^5)+Ir;

%====DAVIDSON AND SCHIFF'S LINEAR MODEL
uo=1;
m_non=m/(0.5*rho*L^3);
xG_non=xG/L;

M=[m_non-Yvdot (m_non*xG_non)-Yrdot;(m_non*xG_non)-Nvdot Iz-Nrdot];
N=[-Yv (m_non*uo)-Yr;-Nv (m_non*xG_non*uo)-Nr];
b1=[0.01;1];

A1=-inv(M)*N;

%=====Davidson and Schiff's Motion Model
Ad=[A1 zeros(2,1);0 1 0];
Bd=[b1;0];
Cd=[0 0 1];
Dd=0;
%=====Check Controllability
M1=[Bd Ad*Bd Ad^2*Bd];
M2=[Bd Ad*Bd Ad^2*Bd Ad^3*Bd Ad^4*Bd];
cek_kontrol1=rank(M1);
cek_kontrol=rank(M2);

%=====Check Observability
O=[Cd;Cd*Ad;Cd*(Ad)^2];
cek_teramat=rank(O);

[A,B,C,D]=c2dm(Ad,Bd,Cd,Dd,Delta_t);
