function [I] = current(Earr,Tarr,mu1,mu2,Temp)
%CURRENT Summary of this function goes here
%   Detailed explanation goes here
f1 = @(E) fermi(E,mu1,Temp);
f2 = @(E) fermi(E,mu2,Temp);
f1arr = arrayfun(f1,Earr);
f2arr = arrayfun(f2,Earr);
sz = size(Earr);
sz = sz(2);
I = 0;
for i = 1:sz-1
    I = I + (Earr(i+1)-Earr(i))*(f2arr(i)-f1arr(i))*Tarr(i);
end

end

