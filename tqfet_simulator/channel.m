function [alpha,beta,H] = channel(NL,NW,params)
%CHANNEL Summary of this function goes here
%   Detailed explanation goes here

%% parameters
sin_30 = 0.5;
cos_30 = sqrt(3)/2;
t = params(1);
yso = params(2);
soi = yso/(3*sqrt(3));
l = params(3);
Ez_factor = params(4);
Ez = Ez_factor*(yso/l);
deltaz = l*Ez;
yr_factor = params(5);
yr = yr_factor*deltaz;
NU = 4;

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


end

