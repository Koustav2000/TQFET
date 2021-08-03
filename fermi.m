function [f]=fermi(E,mu,T)
q = 1.60217662 * 10.^-19;
k = 1.38064852 * 10.^-23;
k1 = k/q;
f = 1/(1+exp((E-mu)/(k1*T)));
end
