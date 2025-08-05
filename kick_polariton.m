clc;
clear;
close all;
%------------------Time
h=0.001;
T=40;
t=0:h:T;
Nt=length(t);
%------------------grid in space and momentum
xmax=200;
N=400;
dx=xmax/N;
x=[-xmax/2:dx:xmax/2-dx];
kx=[0:N/2-1 -N/2:-1]*2*pi/xmax;
k2=kx.^2;
dk=kx(2)-kx(1);
Nx=length(x);
%------------------pumping amplitude in time
f_c=zeros(1,Nt);
f_0=10;%The amplitude is f_0*n
omega_0=1;%The period is 2*pi/omega_0
n=5;
for i=1:n
    f_c=f_c+ 1 +f_0*cos(i*omega_0*t);
end
%f_c and f_cc have the same output.
%in the main part of the code I use f_cc
f_cc = @(tt) 1 + f_0 * sum(arrayfun(@(i) cos(i * omega_0 * tt), 1:5));
%t0=3;
%f_cc=@(tt) exp(-(tt-t0)^2/(0.5^2));
figure;
plot(t,f_c)

%%
%-----------------------cosntants
E_c=0;
E_x=0;
gamma_c=0.1;
gamma_x=0.1;
k=0.3;
Omega=0.03;%3
g_x=0;
V_x=1;
alpha_c=0.5;%hbar/2m
%-----------------------Potential
for i=1:length(x)
        V(i)=V_x*cos(k*x(i));
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
Wpuls=4;
kx=1;
ff=AA*exp(-(x.^2)/(Wpuls^2)).*exp(1i*kx*x);
%-----main------
kk=1;
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
        %-----
        %disp(j);
        if mod(j,100)==0%this part give frames for a movie
            plot(x,abs(psitX).^2);
            pbaspect([1 1 1]);
            %view(0, 90);
            % print(gcf,'-dpng',sprintf('argzc%02d.png',kk))
            % close all;
            kk=kk+1;
            
             [p_c]=gradient(psitC,dx);
             [p_x]=gradient(psitX,dx);
             
        x_mean_c(kk)=trapz(x,conj(psiC).*x.*psiC);
        x_mean_x(kk)=trapz(x,conj(psiX).*x.*psiX);
        
         p_mean_c(kk)=-1i*trapz(x,conj(psiC).*p_c);
         p_mean_x(kk)=-1i*trapz(x,conj(psiX).*p_x);
         
         pop_c(kk)=trapz(x,abs(psitC).^2);
         pop_x(kk)=trapz(x,abs(psitX).^2);
         
         tt(kk)=t(j);
         disp(kk);
       end
end

%%
%mean value of x wavepacket center in time
figure;
plot(tt,abs(x_mean_c));
hold on;
plot(tt,abs(x_mean_x))

%mean value of momentum in time
figure;
plot(tt,abs(p_mean_c));
hold on;
plot(tt,abs(p_mean_x))

% plot(tt,pop_c);
% hold on;
% plot(tt,pop_x)

figure;
scatter(abs(x_mean_c), abs(p_mean_c))
title("phase map for photon")

figure;
scatter(abs(x_mean_x), abs(p_mean_x))
title("phase map for exciton")