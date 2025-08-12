clearvars;
clc;

h=0.02;
T=80;
t=0:h:T;
Nt=length(t);
%------------------grid in space and momentum
xmax=200;%200
N=400;%400
dx=xmax/N;
x=[-xmax/2:dx:xmax/2-dx];
% kx=[0:N/2-1 -N/2:-1]*2*pi/xmax;
kx=[0:N/2-1]*2*pi/xmax;
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
k=0.3;
Omega=5;%3
g_x=0;
V_x=0;
alpha_c=0.5;%hbar/2m
%-----------------------Potential
% for i=1:length(x)
%         V(i)=V_x*cos(k*x(i));
% end
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


for i = 1:length(kx)
    % kc = k(i);
    
    % Photon dispersion
    % omegaC = omegaC0 + hbar * kx(i)^2 / (2 * mC);  % ω_C(k) in rad/s
    omegaC = alpha_c*kx(i)^2;
    
    % Hamiltonian matrix (complex due to decay)
    H = [omegaC + E_c - 1i*gamma_c,  Omega;
         Omega,             E_x - 1i*gamma_x];
     
    % Eigenvalue decomposition
    [V, D] = eig(H);
    
    % Extract eigenvalues
    omega_vals = diag(D);
    
    % Sort by real part (lower and upper polaritons)
    [~, idx] = sort(real(omega_vals));
    omega_minus(i) = omega_vals(idx(1));
    omega_plus(i) = omega_vals(idx(2));
    
    % Corresponding eigenvectors
    v_minus(:,i) = V(:, idx(1)) / norm(V(:, idx(1)));
    v_plus(:,i) = V(:, idx(2)) / norm(V(:, idx(2)));
end


figure;
plot(kx, real(omega_minus), 'b', 'LineWidth', 2); 
hold on;
plot(kx, real(omega_plus), 'r', 'LineWidth', 2);
legend('Lower Polariton', 'Upper Polariton');
title('Exciton-Polariton Dispersion');
grid on;
axis tight


