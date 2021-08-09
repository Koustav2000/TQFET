function [Ek] = bandstructure(alpha,beta,kx_arr)
%BANDSTRUCTURE Summary of this function goes here
%   Detailed explanation goes here
sHk = @(kx) alpha + exp(1i*kx)*beta + exp(-1i*kx)*beta';
Ek = [];
for kx = kx_arr
    val = eig(sHk(kx));
    Ek = [Ek val];
end
sz = size(Ek);
sz = sz(1);
figure
for ii = 1:sz
    plot(kx_arr./pi,Ek(ii,:),'b','LineWidth',1.5)
    hold on;
end
line([0, 2], [0, 0], 'Color', [0,0,0]); 
ylim([-1.6 1.6])
end

