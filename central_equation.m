clc;
clearvars;
close all;

E_c = 1;
E_x = 1;
gamma_c = 0;
gamma_x = 0;
dk = 0.01;

xmax=200;%200
N=400;%400
dx=xmax/N;
x = -xmax/2:dx:xmax/2-dx;
kx = (0:N/2-1)*2*pi/xmax;
k2=kx.^2;
a = 1;% lattice constant

NG = 10;                           % number of G-vectors on each side
Glist = (-NG:NG) * (2*pi/a);     % reciprocal lattice vectors

Omega=10;
g_x=0;
alpha_c=0.5;

V_x = 1;

% for i=1:length(x)
%         % V(i)=(V_x*cos(k*(x(i)-pi/k)) + V_x)/2;
%         V(i)=V_x*cos(k*(x(i)-pi/k));
% end


%%
bands = zeros(2*NG, length(kx));
for i=1:length(kx)
    Pblock = (alpha_c*k2(i) + E_c)*diag(ones(1,NG));

    Oblock = Omega*diag(ones(1,NG));

    Eblock = V_x/2*(diag(ones(1,NG-1),1) + diag(ones(1,NG-1),-1));
    Eblock = Eblock + E_x*diag(ones(1,NG));

    %Boundary condition
    Eblock(1,end) = V_x/2;
    Eblock(end,1) = V_x/2;

    Mat = [Pblock, Oblock;Oblock, Eblock];

    E = eig(Mat);
    E = sort(real(E));
    bands(:, i) = E;

end

%%
figure;
hold on;
for n = 1:2*NG
    plot(kx, bands(n,:), 'b');
end