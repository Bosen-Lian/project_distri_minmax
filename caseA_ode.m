function xdot=caseA_ode(t,x)
global A;
global B1 B2 B3 B4;
global u1 u2 u3 u4;
global Kl1 Kl12 Kl13 Kl14 Kl2 Kl21 Kl23 Kl24 Kl3 Kl31 Kl32 Kl34 Kl4 Kl41 Kl42 Kl43;
global R11 R12 R13 R14 R21 R22 R23 R24 R33 R31 R32 R34 R44 R41 R42 R43;
global Q1 Q2 Q3 Q4 gamma1 gamma3 gamma2 gamma4;


x=[x(1);x(2)];


xdot=[A*x+B1*u1+B2*u2+B3*u3+B4*u4
             -2*(u1+Kl1*x)'*R11*x  % 3,4-th
     2*gamma1^2*(u2+Kl12*x)'*R12*x
     2*gamma1^2*(u3+Kl13*x)'*R13*x
     2*gamma1^2*(u4+Kl14*x)'*R14*x
             -2*(u2+Kl2*x)'*R22*x  % 11,12-th
     2*gamma2^2*(u1+Kl21*x)'*R21*x 
     2*gamma2^2*(u3+Kl23*x)'*R23*x 
     2*gamma2^2*(u4+Kl24*x)'*R24*x  
             -2*(u3+Kl3*x)'*R33*x  % 19,20-th
     2*gamma3^2*(u1+Kl31*x)'*R31*x 
     2*gamma3^2*(u2+Kl32*x)'*R32*x 
     2*gamma3^2*(u4+Kl34*x)'*R34*x 
             -2*(u4+Kl4*x)'*R44*x  % 27,28-th
     2*gamma4^2*(u1+Kl41*x)'*R41*x 
     2*gamma4^2*(u2+Kl42*x)'*R42*x 
     2*gamma4^2*(u3+Kl43*x)'*R43*x 
     x'*Q1*x+(Kl1*x)'*R11*Kl1*x-gamma1^2*(Kl12*x)'*R12*Kl12*x-gamma1^2*(Kl13*x)'*R13*Kl13*x-gamma1^2*(Kl14*x)'*R14*Kl14*x % 35-th
     x'*Q2*x+(Kl2*x)'*R22*Kl2*x-gamma2^2*(Kl21*x)'*R21*Kl21*x-gamma2^2*(Kl23*x)'*R23*Kl23*x-gamma2^2*(Kl24*x)'*R24*Kl24*x
     x'*Q3*x+(Kl3*x)'*R33*Kl3*x-gamma3^2*(Kl31*x)'*R31*Kl31*x-gamma3^2*(Kl32*x)'*R32*Kl32*x-gamma3^2*(Kl34*x)'*R34*Kl34*x
     x'*Q4*x+(Kl4*x)'*R44*Kl4*x-gamma4^2*(Kl41*x)'*R41*Kl41*x-gamma4^2*(Kl42*x)'*R42*Kl42*x-gamma4^2*(Kl43*x)'*R43*Kl43*x];
     

end