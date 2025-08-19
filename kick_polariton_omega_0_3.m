clc;
clearvars;
close all;
%------------------Time

[x,tt, psiC_x0_V_0, psiX_x0_V_0, psiC_t_0, psiX_t_0] = potential(0); 
% [x,tt, psiC_x0_V_1, psiX_x0_V_1, psiC_t_1, psiX_t_1] = potential(1); 
[x,tt, psiC_x0_V_2, psiX_x0_V_2, psiC_t_2, psiX_t_2] = potential(2); 

% for i=1:length(tt)
%     psiC_sum_0(i) = sum(abs(psiC_t_0(i,:)).^2);
%     psiX_sum_0(i) = sum(abs(psiX_t_0(i,:)).^2);
% 
%     psiC_sum_1(i) = sum(abs(psiC_t_1(i,:)).^2);
%     psiX_sum_1(i) = sum(abs(psiX_t_1(i,:)).^2);
% end
%%
% figure()
% plot(tt,psiC_sum_0)
% hold on 
% plot(tt,psiX_sum_0)
% plot(tt,psiC_sum_0+psiX_sum_0,LineWidth=2,Color='k')
% xlabel('t',FontSize=20)
% ylabel('Population',FontSize=20)
% title('Vx=0',FontSize=20)
% legend(["Photon","Exciton","Photon+Exciton"],FontSize=20)

% figure()
% plot(tt,psiC_sum_1)
% hold on 
% plot(tt,psiX_sum_1)
% plot(tt,psiC_sum_1+psiX_sum_1,LineWidth=2,Color='k')
% xlabel('t',FontSize=20)
% ylabel('Population',FontSize=20)
% title('Vx=1',FontSize=20)
% legend(["Photon","Exciton","Photon+Exciton"],FontSize=20)

%%
% figure()
% plot(tt,abs(psiC_x0_V_0).^2)
% hold on
% plot(tt,abs(psiC_x0_V_1).^2)
% plot(tt,abs(psiC_x0_V_2).^2)
% legend(["V=0","V=1","V=2"])
% title("|\Psi_C(x=0,t)|^2",FontSize=20)


% figure()
% plot(tt,abs(psiX_x0_V_0).^2)
% hold on
% plot(tt,abs(psiX_x0_V_1).^2)
% plot(tt,abs(psiX_x0_V_2).^2)
% legend(["V=0","V=1","V=2"])
% title("|\Psi_X(x=0,t)|^2",FontSize=20)

%% 
figure()
pcolor(x,tt,abs(psiC_t_2).^2)
shading interp
xlabel('x',FontSize=20)
ylabel('t',FontSize=20)
title("|\Psi_C(x,t)|^2",FontSize=20)

figure()
pcolor(x,tt,abs(psiX_t_2).^2)
shading interp
xlabel('x',FontSize=20)
ylabel('t',FontSize=20)
title("|\Psi_X(x,t)|^2",FontSize=20)


function [x,tt, psiC_x0, psiX_x0, psiC_t, psiX_t] = potential(V_x) 
h=0.01;
T=80;
t=0:h:T;
Nt=length(t);
%------------------grid in space and momentum
xmax=200;%200
N=400;%400
dx=xmax/N;
x=[-xmax/2:dx:xmax/2-dx];
kx=[0:N/2-1 -N/2:-1]*2*pi/xmax;
k2=kx.^2;
dk=kx(2)-kx(1);
Nx=length(x);
%------------------pumping amplitude in time
f_c=zeros(1,Nt);
f_0=10;%The amplitude is f_0*n
omega_0=0.3;%The period is 2*pi/omega_0
n=5;
for i=1:n
    f_c=f_c+ 1 +f_0*cos(i*omega_0*t);
end
%f_c and f_cc have the same output.
%in the main part of the code I use f_cc
% f_cc = @(tt) 1 + f_0 * sum(arrayfun(@(i) cos(i * omega_0 * tt), 1:5));
t0=0;
f_cc=@(tt) exp(-(tt-t0)^2/(0.5^2));
% figure;
% plot(t,f_c)

%%
%-----------------------cosntants
E_c=0;
E_x=0;
gamma_c=0;
gamma_x=0;
% k=0.1;
k=1;
Omega=10;%3
g_x=0;
% V_x=1;
alpha_c=0.5;%hbar/2m
%-----------------------Potential
for i=1:length(x)
        % V(i)=(V_x*cos(k*(x(i)-pi/k)) + V_x)/2;
        V(i)=V_x*cos(k*(x(i)-pi/k));
end

%------------------------Main
%---initial condition---
% psitC=exp(-(x+6).^2/4^2).*exp(1i*2*x);%zeros(1,N);
% psitC=psitC/sqrt(trapz(x,abs(psitC).^2));
psitC=zeros(1,Nx);
psitX=zeros(1,Nx);
p_mean_c(1)=-1i*trapz(x,conj(psitC).*gradient(psitC,dx));
p_mean_x(1)=-1i*trapz(x,conj(psitX).*gradient(psitX,dx));
x_mean_c(1)=trapz(x,conj(psitC).*x.*psitC);
x_mean_x(1)=trapz(x,conj(psitX).*x.*psitX);

%psitC=psitC/sqrt(trapz(x,trapz(x,abs(psitC).^2,2)));
%-----pulse
AA=1;
Wpuls=8;
kx_pulse=0;
ff=AA*exp(-(x.^2)/(Wpuls^2)).*exp(1i*kx_pulse*x);
%-----main------
kk=1;


figure()
plot(x,V)
hold on
plot(x,ff)


psiC_t = zeros(Nt-1,N);
psix_t = zeros(Nt-1,N);
for j=2:Nt
    psiC=psitC;
    psiX=psitX;
    k1C=(1i*alpha_c)*ifft(-k2.*fft(psiC))...
            -1i*0.5*Omega*psiX...
            -0.5*gamma_c*psiC...
            -1i*E_c*psiC...
            -1i*f_cc(t(j)).*ff;
    k1X=(-1i)*V.*psiX...
            -1i*0.5*Omega*psiC...
            -0.5*gamma_x*psiX...
            -1i*psiX*E_x...
            +(-1i*.5*g_x*(abs(psiX).^2).*psiX);

   %----
    psiC=.5*k1C*h+psitC;
    psiX=.5*k1X*h+psitX;

    k2C=(1i*alpha_c)*ifft(-k2.*fft(psiC))...
            -1i*0.5*Omega*psiX...
            -0.5*gamma_c*psiC...
            -1i*E_c*psiC...
            -1i*f_cc(t(j)+0.5*h).*ff;
    k2X=((-1i)*V.*psiX)...
            +0.5*-1i*Omega*psiC...
            -0.5*gamma_x*psiX...
            -1i*psiX*E_x...
            +(-1i*.5*g_x*(abs(psiX).^2).*psiX);
    %----  
    psiC=0.5*k2C*h+psitC;
    psiX=0.5*k2X*h+psitX;

    k3C=((1i*alpha_c)*ifft(-k2.*fft(psiC)))...
            +0.5*-1i*Omega*psiX...
            -0.5*gamma_c*psiC...
            -1i*E_c*psiC...
            -1i*f_cc(t(j)+0.5*h).*ff;
        
    k3X=((-1i)*V.*psiX)...
            +0.5*-1i*Omega*psiC...
            -0.5*gamma_x*psiX...
            -1i*psiX*E_x...
            +(-.5*g_x*1i*(abs(psiX).^2).*psiX);

%----
    psiC=psitC+k3C*h;
    psiX=psitX+k3X*h;

    k4C=((1i*alpha_c)*ifft(-k2.*fft(psiC)))...
            -1i*0.5*Omega*psiX...
            -0.5*gamma_c*psiC...
            -1i*E_c*psiC...
            -1i*f_cc(t(j)).*ff;
        
    k4X=((-1i)*V.*psiX)...
            -1i*0.5*Omega*psiC...
            -0.5*gamma_x*psiX...
            -1i*psiX*E_x...
            +(-.5*g_x*1i*(abs(psiX).^2).*psiX);

%----
    psitC=psitC+(1/6)*h*(k1C+2*k2C+2*k3C+k4C);
    psitX=psitX+(1/6)*h*(k1X+2*k2X+2*k3X+k4X);

    psiC_t(j,:) = psitC;
    psiX_t(j,:) = psitX;
    
    zero_index = find(x==0);
    psiC_x0(j) = psitC(zero_index);
    psiX_x0(j) = psitX(zero_index);

    kk=kk+1;
    tt(kk)=t(j);
end


end