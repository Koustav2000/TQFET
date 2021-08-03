clear;
clc;
%% ZNR dimensions
NL = 1; 
NW = 8;

Ez_vector = linspace(0,2,20);
for count = 1: length(Ez_vector)
%% parameters
NU = 4;
errmax = 1e-3;
zplus=1i*1e-3;

q = 1.6e-19;
hbar = 1.06e-34;
h = 2*pi*hbar;
quantum = (q*q)/h;

sin_30 = 0.5;
cos_30 = sqrt(3)/2;
t = -1.6;
yso = 0.05*t*3*sqrt(3);
soi = yso/(3*sqrt(3));               
l = 0.4;
Ez = Ez_vector(count)*(yso/l);
deltaz = l*Ez;
%yr = 0.1*t;
yr = 0.0*deltaz;

%% pauli matrices
px = [0 1;1 0];
py = [0 -1i;1i 0];
pz = [1 0;0 -1];
po = [1 0;0 1]; 

%% NN vectors
hop_N = [0,-1];
hop_SE = [-cos_30,sin_30];
hop_SW = [cos_30,sin_30];
hop_S = -hop_N;
hop_NW = -hop_SE;
hop_NE = -hop_SW;
cross = @(hop) px*hop(2)-py*hop(1);

%% alphaU
NN_U = [0 1 0 0;1 0 1 0;0 1 0 1;0 0 1 0];
sub = [1 0 0 0;0 -1 0 0;0 0 1 0;0 0 0 -1];
NNN_U = zeros(4,4);
NNN_U(3,1) = 1i; NNN_U(1,3) = -1i;
NNN_U(4,2) = 1i; NNN_U(2,4) = -1i;
rashbaU = zeros(8,8);
rashbaU(1:2,3:4) = cross(hop_SW);rashbaU(3:4,1:2) = cross(hop_NE);
rashbaU(3:4,5:6) = cross(hop_S);rashbaU(5:6,3:4) = cross(hop_N);
rashbaU(5:6,7:8) = cross(hop_SE);rashbaU(7:8,5:6) = cross(hop_NW);
alphaU = kron(t*NN_U,po);
alphaU = alphaU + kron(deltaz*sub,po);
alphaU = alphaU +  kron(soi*NNN_U,pz);
alphaU = alphaU + 1i*yr*rashbaU;

%% betaN
NN_N = zeros(4,4); NNN_N = zeros(4,4);
NN_N(1,4) = 1;
NNN_N(1,3) = -1i; NNN_N(2,4) = -1i;
rashbaN = zeros(8,8);
rashbaN(1:2,7:8)= cross(hop_N);
betaN = kron(t*NN_N,po);
betaN = betaN + kron(soi*NNN_N,pz);
betaN = betaN + 1i*yr*rashbaN;

%% betaNE
NNN_NE = zeros(4,4);
NNN_NE(1,3) = 1i;
betaNE = kron(soi*NNN_NE,pz);

%% betaE
NN_E = zeros(4,4); NNN_E = zeros(4,4);
NN_E(1,2) = 1; NN_E(4,3) = 1;
NNN_E(1,1) = -1i;NNN_E(2,2) = 1i;NNN_E(3,3) = -1i;NNN_E(4,4) = 1i;
NNN_E(1,3) = 1i;NNN_E(4,2)= -1i;
rashbaE = zeros(8,8);
rashbaE(1:2,3:4)= cross(hop_SE); rashbaE(7:8,5:6)= cross(hop_NE);
betaE = kron(t*NN_E,po);
betaE = betaE + kron(soi*NNN_E,pz);
betaE = betaE + 1i*yr*rashbaE;

%% betaSE
NNN_SE = zeros(4,4);
NNN_SE(4,2) = -1i;
betaSE = kron(soi*NNN_SE,pz);

%% alpha
alpha = kron(eye(NW),alphaU);
alpha = alpha + kron(diag(ones(1,NW-1),1),betaN) + kron(diag(ones(1,NW-1),-1),betaN');

%% beta
beta = kron(eye(NW),betaE);
beta = beta + kron(diag(ones(1,NW-1),1),betaNE) + kron(diag(ones(1,NW-1),-1),betaSE);

%% Hamiltonian

H = kron(diag(ones(1,NL)),alpha) + kron(diag(ones(1,NL-1),1),beta) + kron(diag(ones(1,NL-1),-1),beta');

%% NEGF CALCULATIONS

TLR = 300;
mu1 = -0.1*t;
mu2 = 0.1*t;

low = -0.5;
high = 0.5;

%Define Energy grid
E = linspace(-1,1,81);
fprintf(1,'Energy range = %d eV to %d eV\n',min(E),max(E));
fprintf(1,'%s\n',['[',blanks(length(E)),']']);
fprintf(1,' ');
%Define the matrices for NEGF
T = zeros(1,length(E));
g1 = inv(E(1)*eye(NW*NU*2)-alpha);
g2 = inv(E(1)*eye(NW*NU*2)-alpha);

I(count) = 0;
%Run for each energy
for k = 1:length(E)
    
f1 = fermi(E(k), mu1, TLR);
f2 = fermi(E(k), mu2, TLR);
%Calculate for beta'
 err = 100;
 while(err>errmax)
 g1new = inv((E(k)+zplus)*eye(NW*NU*2)-alpha-beta'*g1*beta);
 err = (sum(sum(abs(g1new-g1))))/(sum(sum(abs(g1new+g1))));
 g1 = g1new;
 end
 sigma1 = beta'*g1*beta;
%Calculate for beta
err = 100;
 while(err>errmax)
 g2new = inv((E(k)+zplus)*eye(NW*NU*2)-alpha-beta*g2*beta');
 err = (sum(sum(abs(g2new-g2))))/(sum(sum(abs(g2new+g2))));
 g2 = g2new;
 end
 sigma2 = beta*g2*beta';
%Calculate self energy matrices
 E1 = kron(diag([1 zeros(1,NL-1)]),sigma1);
 E2 = kron(diag([zeros(1,NL-1) 1]),sigma2);
%Calculate broadening
 G1 = 1i*(E1-E1'); G2 = 1i*(E2-E2');
 
 G=inv((E(k)*eye(NW*(NL)*NU*2))-H-E1-E2);
 T(k) = real(trace(G1*G*G2*G'));
 
 I(count) = I(count) + ((high-low)/length(E))*(f1-f2)*T(k);
 fprintf(1,'|');
end

end

h=plot(Ez_vector,log10(I*quantum));
set(h,'linewidth',3)
set(h,'MarkerSize',10);
hold on

set(gca,'fontsize',40,'linewidth',3,'TickLength',[0.02,0.01]);
%ylim([-13 -10.2]);
%xlabel('V_{GS} (V)','fontsize',50);
%ylabel('log (I_D)','fontsize',50);