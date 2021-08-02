clear;
clc;
sin_30 = 0.5;
cos_30 = sqrt(3)/2;
t = -1.6;
yso = 0.05*t*3*sqrt(3);
soi = yso/(3*sqrt(3));               
l = 0.4;
Ez = 1.5*(yso/l);
deltaz = l*Ez;
yr = 0*0.1*t;

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

%% ZNR dimensions
NL = 2;
NW = 20;

%% alpha
alpha = kron(eye(NW),alphaU);
alpha = alpha + kron(diag(ones(1,NW-1),1),betaN) + kron(diag(ones(1,NW-1),-1),betaN');

%% beta
beta = kron(eye(NW),betaE);
beta = beta + kron(diag(ones(1,NW-1),1),betaNE) + kron(diag(ones(1,NW-1),-1),betaSE);

%% calculation of bandstructure
sHk = @(kx) alpha + exp(1i*kx)*beta + exp(-1i*kx)*beta';
kx_arr = 0*pi:0.01*pi:2*pi;
Ek = [];
for kx = kx_arr
    val = eig(sHk(kx));
    Ek = [Ek val];
end
sz = size(Ek);
sz = sz(1);
figure
for ii = 1:sz
    plot(kx_arr./pi,Ek(ii,:),'LineWidth',1.5)
    hold on;
end
line([0, 2], [0, 0], 'Color', [0,0,0]); 
ylim([-1.6 1.6])
% plot([-1 1],[E-pot E-pot]);
% hold off;

%% calculation of spin-resolved bandstructure
% H = @(kx) alpha + exp(1i*kx)*beta + exp(-1i*kx)*beta';
% B = zeros(2*4*NW,100);
% kx_arr = linspace(0,2*pi,100);
% for i = (1:100)
%     clc;
%     MUT = H(kx_arr(i));
%     Ekx = eig(MUT(1:2:2*4*NW,1:2:2*4*NW));
%     B(1:4*NW,i) = Ekx;
%    
%     Ekx = eig(MUT(2:2:2*4*NW,2:2:2*4*NW));
%     B(4*NW+1:2*4*NW,i) = Ekx;
% 
% end
% figure
% for i = (1:4*NW)
%     plot(kx_arr(1:100)./(2*pi),B(i,1:100),'Color', [1,0,0]);
%     hold on;
% end
% for i = (4*NW+1:2*4*NW)
%     plot(kx_arr(1:100)./(2*pi),B(i,1:100),'Color', [0,1,0]);
%     hold on;
% end
% xlim([0,1]);
% ylim([-5,5]);
% line([0, 1], [0, 0], 'Color', [0,0,0]);   
% % ylabel('E');
% % xlabel('kx');
% hold off;



