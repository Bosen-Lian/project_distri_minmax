close all;
clc;
clear all;

global  A B1 B2 B3 u1 u2 u3 Kl1 Kl12 Kl13 Kl2 Kl21 Kl23 Kl3 Kl31 Kl32...
    R11 R12 R13...
    R21 R22 R23...
    R31 R32 R33...
    Q1 Q2 Q3 gamma1 gamma2 gamma3

A=[-0.0665    8       0         0    ;
    0       -3.663   3.663      0    ;
   -6.86       0    -13.736  -13.736 ;       
   0.6         0      0         0    ];

B1=[0;0;13.736;0];
B2=[-8;0;0;0];
B3=[0;0;-1;0];

R11=1;
R12=1;
R13=2;

R21=1;
R22=2;
R23=2;

R31=2;
R32=3;
R33=1;

gamma1=5;
gamma2=6;
gamma3=7;

Q1=1*eye(4);
Q2=2*eye(4);
Q3=0.5*eye(4);
 
P1e=care(A,[B1,B2,B3],Q1,[R11,0,0;0,-gamma1^2*R12,0;0,0,-gamma1^2*R13]);
K1e=inv(R11)*B1'*P1e;
K12e=-1/gamma1^2*inv(R12)*B2'*P1e;
K13e=-1/gamma1^2*inv(R13)*B3'*P1e;

P2e=care(A,[B1,B2,B3],Q2,[-gamma2^2*R21,0,0;0,R22,0;0,0,-gamma2^2*R23]);
K2e=inv(R22)*B2'*P2e;
K21e=-1/gamma2^2*inv(R21)*B1'*P2e;
K23e=-1/gamma2^2*inv(R23)*B3'*P2e;
      
P3e=care(A,[B1,B2,B3],Q3,[-gamma3^2*R31,0,0;0,-gamma3^2*R32,0;0,0,R33]);
K3e=inv(R33)*B3'*P3e;
K31e=-1/gamma3^2*inv(R31)*B1'*P3e;
K32e=-1/gamma3^2*inv(R32)*B2'*P3e;
     
x=[1 -1  1 -1]';
i=1;
ii=0;

Kl1=[1 2 0 0];
Kl12=[0.1 0.1 0.01 0.1];
Kl13=[0.0001 0.0001 0.0001 0.0001];

Kl2=[-1 0 0 -1];
Kl21=[-0.001 0 0 -0.001];
Kl23=[-0.001 0 0 -0.001];
 

Kl3=[0 -1 0 0];
Kl31=[0 -0.001 0 0];
Kl32=[0 -0.001 0 0];

%  Kl1=[1 2 0 0];
%   Kl12=[0.1 0.1 0 0];
%   Kl13=[0.1 0.1 0 0];
% 
%   Kl2=[-1 0 0 -1];
%   Kl21=[-0.1 0 0 -0.1];
%   Kl23=[-0.1 0 0 -0.1];
% 
%   Kl3=[0 -1 0 0];
%   Kl31=[0 -0.1 0 0];
%   Kl32=[0 -0.1 0 0];
  
K1=[2 3 0 0];
K2=[-2 0 0 -2];
K3=[0 -2 0 0];


eKl1=[];
eKl12=[];
eKl13=[];
eKl2=[];
eKl21=[];
eKl23=[];
eKl3=[];
eKl31=[];
eKl32=[];
dKl1=0.5;
dKl2=0.5;
dKl3=0.5;

KK1=[];
KK2=[];
KK3=[];

kk=0;
T=0.01;

iii=0;

itrr=200;

tic
for k=1:itrr
    
      Kl1last=Kl1;
      Kl2last=Kl2;
      Kl3last=Kl3;
     
    e1=0.00111*rand(1);
    e2=0.00111*rand(1);
    e3=0.00111*rand(1);
    
     u1=-K1*x(:,k)+e1;
     u2=-K2*x(:,k)+e2;
     u3=-K3*x(:,k)+e3;
     
    tspanxd=[0 T];
    X(:,k)=[x(:,k)',zeros(1,39)]';
    [tl,dl]= ode45(@threeplayers_power,tspanxd,X(:,k));
    x(:,k+1)=[dl(length(tl),1:4)]; 
    
    

    iii=iii+1;
    dxx(iii,:)=[x(1,k+1)^2,2*x(1,k+1)*x(2,k+1),2*x(1,k+1)*x(3,k+1),2*x(1,k+1)*x(4,k+1),x(2,k+1)^2,2*x(2,k+1)*x(3,k+1),2*x(2,k+1)*x(4,k+1),x(3,k+1)^2,2*x(3,k+1)*x(4,k+1),x(4,k+1)^2]-[x(1,k)^2,2*x(1,k)*x(2,k),2*x(1,k)*x(3,k),2*x(1,k)*x(4,k),x(2,k)^2,2*x(2,k)*x(3,k),2*x(2,k)*x(4,k),x(3,k)^2,2*x(3,k)*x(4,k),x(4,k)^2];
              
     lk1(iii,:)=dl(length(tl),5:8);
    lk12(iii,:)=dl(length(tl),9:12);
    lk13(iii,:)=dl(length(tl),13:16);
      lk2(iii,:)=dl(length(tl),17:20);
    lk21(iii,:)=dl(length(tl),21:24);
    lk23(iii,:)=dl(length(tl),25:28);
      lk3(iii,:)=dl(length(tl),29:32);
    lk31(iii,:)=dl(length(tl),33:36);
    lk32(iii,:)=dl(length(tl),37:40);
    
    
    fai1(iii,:)=[dxx(iii,:),lk1(iii,:),lk12(iii,:),lk13(iii,:)];
    fai2(iii,:)=[dxx(iii,:),lk2(iii,:),lk21(iii,:),lk23(iii,:)];
    fai3(iii,:)=[dxx(iii,:),lk3(iii,:),lk31(iii,:),lk32(iii,:)];
 
    
    r1(iii)=-1*dl(length(tl),41);
    r2(iii)=-1*dl(length(tl),42);  
    r3(iii)=-1*dl(length(tl),43);
  
    
    a1=rank(fai1);
    a2=rank(fai2);
    a3=rank(fai3);
    

  if dKl1(end)>0.0001 || dKl2(end)>0.0001 || dKl3(end)>0.0001  
    
     if a1==22 && a2==22 && a3==22 
        
                kk=kk+1;
       
         p1=fai1\r1';
         p2=fai2\r2';
         p3=fai3\r3';
      
         
         P1=[p1(1) p1(2) p1(3) p1(4)
             p1(2) p1(5) p1(6) p1(7)
             p1(3) p1(6) p1(8) p1(9)
             p1(4) p1(7) p1(9) p1(10)];
         Kl1=[ p1(11)  p1(12) p1(13)  p1(14)];
         Kl12=[p1(15)  p1(16) p1(17)  p1(18)];
         Kl13=[p1(19)  p1(20) p1(21)  p1(22)];
         
         
         KK1=[KK1; Kl1];
         
         P2=[p2(1) p2(2) p2(3) p2(4)
             p2(2) p2(5) p2(6) p2(7)
             p2(3) p2(6) p2(8) p2(9)
             p2(4) p2(7) p2(9) p2(10)];
         Kl2=[ p2(11)  p2(12) p2(13)  p2(14)];
         Kl21=[p2(15)  p2(16) p2(17)  p2(18)];
         Kl23=[p2(19)  p2(20) p2(21)  p2(22)];
                  
         KK2=[KK2; Kl2];
         
         
         P3=[p3(1) p3(2) p3(3) p3(4)
             p3(2) p3(5) p3(6) p3(7)
             p3(3) p3(6) p3(8) p3(9)
             p3(4) p3(7) p3(9) p3(10)];
         Kl3=[ p3(11)  p3(12) p3(13)  p3(14)];
         Kl31=[p3(15)  p3(16) p3(17)  p3(18)];
         Kl32=[p3(19)  p3(20) p3(21)  p3(22)];
         
         KK3=[KK3; Kl3];
         
                  
         eKl1(kk)=norm(Kl1-K1e);
        eKl12(kk)=norm(Kl12-K12e);
        eKl13(kk)=norm(Kl13-K13e);
       
         eKl2(kk)=norm(Kl2-K2e);
        eKl21(kk)=norm(Kl21-K21e);
        eKl23(kk)=norm(Kl23-K23e);
        
         eKl3(kk)=norm(Kl3-K3e);
        eKl31(kk)=norm(Kl31-K31e);
        eKl32(kk)=norm(Kl32-K32e);
        
         
         eP1(kk)=norm(P1-P1e);
         eP2(kk)=norm(P2-P2e);
         eP3(kk)=norm(P3-P3e);
        
        
        
      dKl1(kk)=norm(Kl1-Kl1last);
      dKl2(kk)=norm(Kl2-Kl2last);
      dKl3(kk)=norm(Kl3-Kl3last);
     
        
         a1=0;
         a2=0;
         a3=0;
         a4=0;    
         
         lk1=[];
         lk12=[];
         lk13=[];
             
         lk2=[];
         lk21=[];
         lk23=[];
      
          lk3=[];
         lk31=[];
         lk32=[];
     
         dxx=[];  
         fai1=[];
         fai2=[];
         fai3=[];
         
         r1=[];
         r2=[];
         r3=[];
         
         p1=[];
         p2=[];
         p3=[];
        
         
         iii=0;
            
        
         
     end
        
   end

        
     if dKl1(end)<=0.0001 &&  dKl2(end)<=0.0001 && dKl3(end)<=0.0001 
         
        break;
     end      
          
end
toc

    

figure  
plot(eP1,'-d','LineWidth',2);
hold on
plot(eP2,'-d','LineWidth',2);
hold on
plot(eP3,'-d','LineWidth',2);
xlabel('Time iteration $h$');legend(' $\Vert P_1^h-P_1 \Vert$',' $\Vert P_2^h-P_2 \Vert$',' $\Vert P_3^h-P_3 \Vert$');
set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'XLabel'),'Interpreter','latex');
set(get(gca,'Legend'),'FontSize',26);
set(legend,'Interpreter','latex');

   
 figure  
plot(eKl1,'-d','LineWidth',2);
hold on
plot(eKl2,'-d','LineWidth',2);
hold on
plot(eKl3,'-d','LineWidth',2);
xlabel('Time iteration $h$');legend(' $\Vert K_1^h-K_{1}^* \Vert$',' $\Vert K_{2}^h-K_{2}^* \Vert$',' $\Vert K_{3}^h-K_{3}^* \Vert$');
set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'XLabel'),'Interpreter','latex');
set(get(gca,'Legend'),'FontSize',26);
set(legend,'Interpreter','latex');


 figure  
plot(eKl12,'-d','LineWidth',2);
hold on
plot(eKl13,'-d','LineWidth',2);
xlabel('Time iteration $h$');legend(' $\Vert K_{12}^h-K_{12}^* \Vert$',' $\Vert K_{13}^h-K_{13}^* \Vert$');
set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'XLabel'),'Interpreter','latex');
set(get(gca,'Legend'),'FontSize',26);
set(legend,'Interpreter','latex');


     figure  
plot(eKl21,'-d','LineWidth',2);
hold on
plot(eKl23,'-d','LineWidth',2);
xlabel('Time iteration $h$');legend(' $\Vert K_{21}^h-K_{21}^* \Vert$',' $\Vert K_{23}^h-K_{23}^* \Vert$');
set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'XLabel'),'Interpreter','latex');
set(get(gca,'Legend'),'FontSize',26);
set(legend,'Interpreter','latex');

figure  
plot(eKl31,'-d','LineWidth',2);
hold on
plot(eKl32,'-d','LineWidth',2);
xlabel('Time iteration $h$');legend('$\Vert K_{31}^h-K_{31}^* \Vert$','$\Vert K_{32}^h-K_{32}^* \Vert$');
set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'XLabel'),'Interpreter','latex');
set(get(gca,'Legend'),'FontSize',26);
set(legend,'Interpreter','latex');


     
uu1=[];
uu2=[];
uu3=[];


xx = [rand(1,4)]';

itlrr=150;
Tl=0.08;
 for l=1:itlrr

    tspand=[0*Tl,1*Tl];

 if l<=kk
            
    u1 = -KK1(l,:)*xx(:,l);
    u2 = -KK2(l,:)*xx(:,l);
    u3 = -KK3(l,:)*xx(:,l);
    
    
 else
     
    u1 = -KK1(end,:)*xx(:,l);
    u2 = -KK2(end,:)*xx(:,l);
    u3 = -KK3(end,:)*xx(:,l);
         
 end   
    
    uu1=[uu1 u1];
    uu2=[uu2 u2];
    uu3=[uu3 u3];
    
     
    [t,d]= ode23('dynamic_power',tspand,xx(:,l));
    xx(:,l+1)=[d(length(t),1:4)]; 
    
 end






% figure
% t=0:T:itlrr*T;
% subplot(4,1,1); 
% plot(t,xx(1,:),'r','Linewidth',2);
% legend('$x_1$');
% set(get(gca,'XLabel'),'FontSize',24);
% set(get(gca,'Legend'),'FontSize',26);
% set(legend,'Interpreter','latex');
% subplot(4,1,2);
% plot(t,xx(2,:),'Linewidth',2);
% xlabel('Time (s)');
% legend('$x_2$');
% set(get(gca,'XLabel'),'FontSize',24);
% set(get(gca,'Legend'),'FontSize',26);
% set(legend,'Interpreter','latex');
% subplot(4,1,3);
% plot(t,xx(3,:),'Linewidth',2);
% xlabel('Time (s)');
% legend('$x_3$');
% set(get(gca,'XLabel'),'FontSize',24);
% set(get(gca,'Legend'),'FontSize',26);
% set(legend,'Interpreter','latex');
% subplot(4,1,4);
% plot(t,xx(4,:),'Linewidth',2);
% xlabel('Time (s)');
% legend('$x_4$');
% set(get(gca,'XLabel'),'FontSize',24);
% set(get(gca,'Legend'),'FontSize',26);
% set(legend,'Interpreter','latex');
% 
% 
% 
% 
% figure
%   tt=T:T:itlrr*T;
%   subplot(3,1,1); 
%   plot(tt,uu1,'-.','Linewidth',2);
%   legend('$u_1$');
% set(get(gca,'XLabel'),'FontSize',24);
% set(get(gca,'XLabel'),'Interpreter','latex');
% set(get(gca,'Legend'),'FontSize',26);
% set(legend,'Interpreter','latex');
% subplot(3,1,2); 
%  plot(tt,uu2,'-.','Linewidth',2);
%  legend('$u_2$');
% set(get(gca,'XLabel'),'FontSize',24);
% set(get(gca,'XLabel'),'Interpreter','latex');
% set(get(gca,'Legend'),'FontSize',26);
% set(legend,'Interpreter','latex');
% subplot(3,1,3); 
%  plot(tt,uu3,'-.','Linewidth',2);
%  legend('$u_3$');
% set(get(gca,'XLabel'),'FontSize',24);
% set(get(gca,'XLabel'),'Interpreter','latex');
% set(get(gca,'Legend'),'FontSize',26);
% set(legend,'Interpreter','latex');



figure
plot(xx(1,:),'Linewidth',2)
hold on
plot(xx(2,:),'Linewidth',2)
hold on
plot(xx(3,:),'Linewidth',2)
hold on
plot(xx(4,:),'Linewidth',2)
legend('x1','x2')

figure
plot(uu1,'Linewidth',2)
hold on
plot(uu2,'Linewidth',2)
hold on
plot(uu3,'Linewidth',2)
legend('u1','u2','u3')


