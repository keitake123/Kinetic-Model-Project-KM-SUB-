% X = gas molecules (e.g., ozone), Y = matrix (e.g., organics)
clear all;
%%%%%%%%%%%%%%%%%%%constants%%%%%%%%%%%%%%%%%%
global pi R W Xgs k11 sozone T H
global k2 Td alphas0
global rp L delta_O3 delta_ole r A0 delta
global Dx Dy kbs kbss kssb kbbx kbby
global A V
pi=3.14159;
R=82.0578;%gas constant
%Radius and Number of Layers in the bulk
rp=2.5e-5;%[cm] particle radius (effective thickness 8.4e-6) (real thickness 0.01045)
L=100;% number of bulk layers
T=296;%[K] Temperature
Xgs=5.658e15;%[cm-3] gas phase concentration of X (230ppmv:5.658e15) (24.5ppmv:6.02e14 )
W=3.6e4;%[cm s-1] thermal velocity of X
sozone=1.724e-15;%[cm2] effective molecular cross section of X
delta_O3=0.4e-7;%[cm] molecular diameter of X (thickness of surface layer)
delta_ole=1.2e-7;%cm molecular diameter of Y (thickness of quasi-static surface layer)
%Input kinetic parameters
alphas0=1e-3;%surface accomodation coefficient of X on free substrate
Td=1e-8;%[s] desorption lifetime of X
H=4.8e-4*R*T;%Henry's law coefficient of X in Y
k11=1e-12;%[cm2 s-1] surface reaction rate coefficient between X and Y
k2=1e-17;%[cm3 s-1] bulk reaction rate coefficient between X and Y
Dx=1e-8;%[cm2 s-1] diffusion coefficient of X
Dy=1e-17;%[cm2s-1] diffusion doefficient of Y
%--------------------------------------------------------
r=rp-delta_ole;%calculate bulk radius
A0=4*pi*r^2;%calculate surface ares
delta=r/L;%[cm] thickness of one bulk layer
kbs=2*Dx/((delta_O3+delta)/2+delta_ole);%(8/pi)*Dx/((delta_O3+delta)/2+delta_ole);%[cm s-1] transport velocity of X from bulk to surface
kbss=2*Dy/(delta_ole+delta);%(8/pi)*Dy/(delta_ole+delta);%[cm s-1]transport velocity of Y from bulk to surafce
kssb=kbss/delta_ole;%[cm s-1] transport velocity of Y from surface to bulk
kbbx=Dx/delta;%(4/pi)*Dx/delta;%[cm s-1] bulk transport velocity of X
kbby=Dy/delta;%(4/pi)*Dy/delta;%[cm s-1] bulk transport velocity of Y
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n=1:L
    %spherical KM-SUB(particle, aerosol flow tube experiments)
    V(n)=4/3*pi*((r-(n-1)*delta)^3-(r-n*delta)^3);%calculate volume of bulk layer n
    A(n)=4*pi*(r-(n-1)*delta)^2;%calculate surface area of bulk layer n
    %Plane KM-SUB (thin film, coated wall flow tube experiments)
    %A(n)=1;
    %V(n)=A(n)*delta;
    Ratio(n)=A(n)./V(n);
end
n=999;
% ***** Initial Conditions **********************
%t=linspace(0.001,180,n);%input the simulation time
t=logspace(-3,2,n);
%y(1):[O3]s,y(2):[oleic]ss,y(3):[O3]bs,y(4):[oleic]bs
%y(5):[O3]b1,y(6):[oleic]b1,y(7):[O3]b2,y(8):[oleic]b2
%...
%y(21):[O3]b9,y(22):[oleic]b9,y(23):[O3]b10,y(24):[oleic]b10
%Initial concentration of X on surface(y0(1)) and in the each bulk layer
for i=1:L+1
    y0(2*i-1)=0;
end
%Initial surface concentration of Y
y0(2)=9.82e13; %9.82e13
%Initial concentration of Y in each bulk layer
for i=2:L+1
    y0(2*i)=8.1872e20;
end
%Call function to solve differential equations
[X,Y] = ode23tb(@KMSUB_F,t,y0);
%calculate transport velocity of X from surface to bulk (this is time-dependent)
ksb=H.*kbs./Td./(W.*alphas0.*(1-sozone.*Y(1:n,1))./4);
%********
time=transpose(t);
tetas=Y(1:n,1)*sozone;%surface coverage of X
tetass=Y(1:n,2)*1.8e-15;%surface coverage of Y
alphas=alphas0*(1-tetas);%surface accomodation coefficient of X
ka=alphas*W/4;%adsorption coefficient of X
ka0=alphas0*W/4;%initial adsorption coefficient of X
ks=k11*Y(1:n,2);%pseudo-first order surface reaction rate coefficient
kd=1/Td;%desorption coefficient of X
Jssb=kssb*Y(1:n,2);%Flux of Y from surface to bulk
Jbss=kbss*Y(1:n,4);%Flux of Y from bulk to surface
Jssbnet=Jssb-Jbss;%Net flux of Y from surface to bulk
Jsb=ksb.*Y(1:n,1);%Flux of X from surface to bulk
Jbs=kbs.*Y(1:n,3);%Flux of X from bulk to surface
Jdes=Y(1:n,1)/Td;%Desorption flux of X
Jcoll=W/4*Xgs;%Collision flux of X
Jads=alphas*Jcoll;%Adsorption flux of X
Ls=k11*Y(1:n,1).*Y(1:n,2);%Loss rate of surface reaction
alphab=alphas.*Jsb./(Jsb+Jdes+Ls);%bulk accomodation coefficient of X
gamma=(Jads-Jdes)/Jcoll;%Uptake coefficient of X
G1e9=gamma*1e9;%Uptake coefficient*1e9
%loss rate (turnover) in layer n =1-6 [cm-3 s-1]
for i=1:6
    if i==1
        Lb6(1:n,i)=k2*Y(1:n,3).*Y(1:n,4);
    else
        %Lb6(1:n,i)=k2*Y(1:n,1+2*L*((i-1)/5)).*Y(1:n,2+2*L*((i-1)/5));
        Lb6(1:n,i)=k2*Y(1:n,2*i+1).*Y(1:n,2*i+2);
    end
end
%Absolute loss rate in layer n =1-6 [s-1]
for i=1:6
    if i==1
        Lb_S(1:n,i)=k2*Y(1:n,3).*Y(1:n,4)*V(1);
    else
        Lb_S(1:n,i)=k2*Y(1:n,2*i+1).*Y(1:n,2*i+2)*V(i);
    end
end
%Absolute loss rate in each bulk layer [s-1]
for i=1:L
    LbV(1:n,i)=k2*Y(1:n,1+2*i).*Y(1:n,2+2*i)*V(i);
end
%Loss rate in each bulk layer [cm-3 s-1]
for i=1:L
    Lb(1:n,i)=k2*Y(1:n,1+2*i).*Y(1:n,2+2*i);
end
for i=1:L
    Lb_t(:,i)=Lb(:,L+1-i);
end
for i=1:n
    LbV_T(i,1)=LbV(i,1);
    for j=2:L
        LbV_T(i,j)=LbV_T(i,j-1)+LbV(i,j);
    end
end
%create r/rp
for i=1:L
    R_Rp(i)=(r-delta*(i-1))/rp;
    R_Rp_t(i)=delta*(i-1)/rp;
end
for i=1:n
    R_LbV(i)=interp1q(LbV_T(i,1:L),R_Rp',LbV_T(i,L)*(1-1/2.718));
end
for i=1:n
    R_LbV2(i)=interp1q(LbV_T(i,1:L),R_Rp',LbV_T(i,L)/2);
end
for i=1:n
    R_Lb(i)=interp1q(Lb_t(i,1:L),R_Rp_t',(Lb_t(i,L)+Lb_t(i,1))/2);
end
R_LbV1=R_LbV';
R_LbV22=R_LbV2';
R_Lb1=R_Lb';
%Ny=Ass*Yss+V(1)*Y(4)+V(2)*Y(6)+V(3)*Y(8)...+V(n)*Ybn
Ny=A0.*Y(1:n,2);
for i=1:L
    Ny=Ny+V(i).*Y(1:n,2*i+2);
end
r_cat=Ny./Ny(1);
%Creat separate vector for ozone and oleic acid
Y_O3=Y(1:n,3);
for i=2:L
    Y_O3=[Y_O3 Y(1:n,2*i+1)];
end
Y_oleic=Y(1:n,4);
for i=2:L
    Y_oleic=[Y_oleic Y(1:n,2*i+2)];
end
%Change the row of metrix of X to the other way around
for i=1:L
    Y_O3_t(:,i)=Y_O3(:,L+1-i);
end
%r/rp value when bulk concentration of X is average of surface and core
%concentration of X
for i=1:n
    O3_half(i)=interp1q(Y_O3_t(i,1:L),R_Rp_t',(Y_O3_t(i,1)+Y_O3_t(i,L))/2);
end
%r/rp value when bulk concentration of X is 1/e of surface and core
%concentration of X
for i=1:n
    O3_halfe(i)=interp1q(Y_O3_t(i,1:L),R_Rp_t',(Y_O3_t(i,1)+Y_O3_t(i,L))*(1-1/2.718));
end
O3_half2=O3_halfe';
%Figures
figure(1);
loglog(t,Y(1:n,1),t,Y(1:n,2));
title('Surface depletion of ozone and oleic acid over time');
xlabel('time (sec)');ylabel ('Xs, Yss (molecule cm^-2)');
%axis([1e-2,1e4,1e9,1e16]);
legend('[O3]s','[oleic acid]ss');
figure(2);
loglog(t,Y(1:n,3),t,Y(1:n,2*L/5+1),t,Y(1:n,2*2*L/5+1),t,Y(1:n,2*3/5*L+1)...
    ,t,Y(1:n,2*4/5*L+1),t,Y(1:n,2*L+1));
title('Bulk surface depletion of ozone');
xlabel('time (sec)');ylabel ('[O3]bs (molecule cm-3)');
legend('[O3]bs','[O3]b1','[O3]b2','[O3]b3','[O3]b4','[O3]b5(core)');
%axis([1e-2,1e2,1e14,1e15]);
figure(3);
plot(t,Y(1:n,4),t,Y(1:n,2*L/5+2),t,Y(1:n,2*2*L/5+2),t,Y(1:n,2*3*L/5+2)...
    ,t,Y(1:n,2*4*L/5+2),t,Y(1:n,2*L+2));
title('Bulk depletion of oleic acid over time');
xlabel('time (sec)');ylabel ('[oleic]b (molecule cm-3)');
legend('[Y]bs','[Y]b1','[Y]b2','[Y]b3','[Y]b4','[Y]b5(core)');
figure(4);
loglog(t,alphas,t,gamma,t,alphab);
title('alphas,gamma');
legend('alphas','gamma','alphab');
%axis([1e-2,1e2,1e-5,0.1]);
figure(5);
loglog(t,Jsb,t,Jbs,t,Jdes,t,Ls);
title('Jsb,Jbs,Jdes,Ls');
legend('Jsb','Jbs','Jdes','Ls');
%Experimental data Ziemann 2005
t_data=[0 1.3 2.5 4.1 5.2 6.38 7.5 9 10.3 11.5 13.1 14.3 15.4];
Ny_data=[1 0.97516 0.93848 0.8836 0.85562 0.81 0.76563 0.7921 0.68063 0.65004 0.60062 0.54391 0.52563];
load catecholdecay.mat
xaxis=catecholdecay.xaxis;
yaxis=catecholdecay.yaxis;
load cat70RH2.mat
xaxis2=cat70RH2.xaxis2;
yaxis2=cat70RH2.yaxis2;
figure(6);
plot(t_data,Ny_data,'o',t,Ny);
title('Ny');
legend('Ny data','Ny model');
%axis([0,40,0,4.5e7]);
figure(7);
plot(xaxis,yaxis,'o',t,r_cat)
title('Ny');
legend('Ny data','Ny model');
%figure(8);
%plot(xaxis2,yaxis2,'o',t,r_cat)
%title('Ny');
%legend('Ny data','Ny model');