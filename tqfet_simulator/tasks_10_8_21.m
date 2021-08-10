clear;
clc;
close all;

tic;
count = 0;
for Vg_norm = 0:0.05:2
count = count+1;

t = -1.6;
yso = 0.05*t*3*sqrt(3);
l = 0.4;
Ez_factor = Vg_norm;
yr_factor = 0.0;
params = [t yso l Ez_factor yr_factor];
NL = 1;
NW = 8;
NU = 4;
kx_arr = 0*pi:0.01*pi:2*pi;
tol = [1e-3 1e-3];
Earr = linspace(-0.7,0.7,71);
Elen = size(Earr);
Elen = Elen(2);
q = 1.6e-19;
hbar = 1.06e-34;
h = 2*pi*hbar;
quantum = (q*q)/h;
Temp = 300;
mu1 = 0.1*t;
mu2 = -0.1*t;


[alphaCH,betaCH,H] = channel(NL,NW,params);
[alphaLD,betaLD] = lead(NW,params);
% Ek = bandstructure(alphaCH,betaCH,kx_arr);

Tarr = [];
if 1
M1 = [];
M2 = [];
g1 = inv(Earr(1)*eye(NW*NU*2)-alphaCH);
g2 = inv(Earr(1)*eye(NW*NU*2)-alphaCH);
for E = Earr
    [sigma1,sigma2,g1new,g2new] = self_RGF(E,alphaCH,betaCH,tol,NW,g1,g2);
    g1 = g1new;
    g2 = g2new;
    M1 = [M1 sigma1];
    M2 = [M2 sigma2];
    [Gr,T] = NEGF(E,sigma1,sigma2,H,NW,NL);
    Tarr = [Tarr T];
%     disp(T)
end
end

% figure
% plot(Earr, Tarr, 'LineWidth',2)
% if count ~= 1
% for i = 1:Elen
%     sigma1 = M1(1:NW*NU*2,1+(i-1)*NW*NU*2:i*NW*NU*2);
%     sigma2 = M2(1:NW*NU*2,1+(i-1)*NW*NU*2:i*NW*NU*2);
%     [Gr,T] = NEGF(Earr(i),sigma1,sigma2,H,NW,NL);
%     Tarr = [Tarr T];
%     disp(T)
% end
% end

I(count) = current(Earr,Tarr,mu1,mu2,Temp)*quantum;
disp(I(count))
V(count) = abs(Vg_norm*yso);


end
toc;
semilogy(V,I,'LineWidth',2);