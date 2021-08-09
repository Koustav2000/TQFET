function [Gr,T] = NEGF(E,sigma1,sigma2,H,NW,NL)
%NEGF Summary of this function goes here
%   Detailed explanation goes here
%Calculate self energy matrices
    NU = 4;
 E1 = kron(diag([1 zeros(1,NL-1)]),sigma1);
 E2 = kron(diag([zeros(1,NL-1) 1]),sigma2);
%Calculate broadening
 T1 = 1i*(E1-E1'); T2 = 1i*(E2-E2');
 
Gr=inv((E*eye(NW*(NL)*NU*2))-H-E1-E2);
Ga = Gr';
T = real(trace(T1(1:NW*NU*2,1:NW*NU*2)*Gr(1:NW*NU*2,1+(NL-1)*NW*NU*2:NL*NW*NU*2)*T2(1+(NL-1)*NW*NU*2:NL*NW*NU*2,1+(NL-1)*NW*NU*2:NL*NW*NU*2)*Ga(1+(NL-1)*NW*NU*2:NL*NW*NU*2,1:NW*NU*2)));
% T = real(trace(T1*Gr*T2*Ga));
end

