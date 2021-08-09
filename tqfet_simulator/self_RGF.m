function [sigma1,sigma2] = self_RGF(E,alpha,beta,tol,NW)
%NEGF Summary of this function goes here
%   Detailed explanation goes here
%Define the matrices for NEGF
errmax = tol(1);
zplus = 1i*tol(2);
NU = 4;
g1 = inv(E*eye(NW*NU*2)-alpha);
g2 = inv(E*eye(NW*NU*2)-alpha);
%Calculate surface GFs self consistently
%Calculate for beta'
 err = 100;
 while(err>errmax)
 g1new = inv((E+zplus)*eye(NW*NU*2)-alpha-beta'*g1*beta);
 err = (sum(sum(abs(g1new-g1))))/(sum(sum(abs(g1new+g1))));
 g1 = g1new;
 end
 sigma1 = beta'*g1*beta;
%Calculate for beta
err = 100;
 while(err>errmax)
 g2new = inv((E+zplus)*eye(NW*NU*2)-alpha-beta*g2*beta');
 err = (sum(sum(abs(g2new-g2))))/(sum(sum(abs(g2new+g2))));
 g2 = g2new;
 end
 sigma2 = beta*g2*beta';
end

