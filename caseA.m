close all;
clc;
clear all;



global  A B1 B2 B3 B4 u1 u2 u3 u4 Kl1 Kl12 Kl13 Kl14 Kl2 Kl21 Kl23 Kl24 Kl3 Kl31 Kl32 Kl34...
    Kl4 Kl41 Kl42 Kl43...
    R11 R12 R13 R14...
    R21 R22 R23 R24...
    R31 R32 R33 R34...
    R41 R42 R43 R44 Q1 Q2 Q3 Q4 gamma1 gamma3 gamma2 gamma4







A=[0         1;
    31.5397  0];

B1=[0;-4];
B2=[1;-1];
B3=[-1;1.7];
B4=[1;0];


R11=1;
R12=2;
R13=3;
R14=1;


R21=2;
R22=1;
R23=3;
R24=1;

R31=2;
R32=3;
R33=1;
R34=1;

R41=1;
R42=2;
R43=1;
R44=1;


gamma1=5;
gamma2=6;
gamma3=7;
gamma4=6;


 Q1=1*eye(2);
 Q2=2*eye(2);
 Q3=0.5*eye(2);
 Q4=0.25*eye(2);


 
  P1e=care(A,[B1,B2,B3,B4],Q1,[R11,0,0,0;0,-gamma1^2*R12,0,0;0,0,-gamma1^2*R13,0;0,0,0,-gamma1^2*R14])
         K1e=inv(R11)*B1'*P1e;
       K12e=-1/gamma1^2*inv(R12)*B2'*P1e;
       K13e=-1/gamma1^2*inv(R13)*B3'*P1e;
       K14e=-1/gamma1^2*inv(R14)*B4'*P1e;
 
 
  P2e=care(A,[B1,B2,B3,B4],Q2,[-gamma2^2*R21,0,0,0;0,R22,0,0;0,0,-gamma2^2*R23,0;0,0,0,-gamma2^2*R24])
 
        K2e=inv(R22)*B2'*P2e;
       K21e=-1/gamma2^2*inv(R21)*B1'*P2e;
       K23e=-1/gamma2^2*inv(R23)*B3'*P2e;
       K24e=-1/gamma2^2*inv(R24)*B4'*P2e; 
  
  
  P3e=care(A,[B1,B2,B3,B4],Q3,[-gamma3^2*R31,0,0,0;0,-gamma3^2*R32,0,0;0,0,R33,0;0,0,0,-gamma3^2*R34])
 

       K3e=inv(R33)*B3'*P3e;
       K31e=-1/gamma3^2*inv(R31)*B1'*P3e;
       K32e=-1/gamma3^2*inv(R32)*B2'*P3e;
       K34e=-1/gamma3^2*inv(R34)*B4'*P3e;
       
   P4e=care(A,[B1,B2,B3,B4],Q4,[-gamma4^2*R41,0,0,0;0,-gamma4^2*R42,0,0;0,0,-gamma4^2*R43,0;0,0,0,R44])
 
       K4e=inv(R44)*B4'*P4e;
       K41e=-1/gamma4^2*inv(R41)*B1'*P4e;
       K42e=-1/gamma4^2*inv(R42)*B2'*P4e;
       K43e=-1/gamma4^2*inv(R43)*B3'*P4e;  
 
 
 
 

x=[1 -1]';


i=1;
ii=0;




 Kl1=[-20  -2];
% Kl12=[-0.6  -0.1];
% Kl13=[0.2  0.05];
% Kl14=[-1  -0.1];
Kl12=[0  0];
Kl13=[0  0];
Kl14=[0  0];


Kl2=[13  2];
Kl21=[0.1  0.01];
Kl23=[0.1  0.01];
Kl24=[-0.2 -0.1];
 

Kl3=[-11  -2];
Kl31=[0.1  0.001];
Kl32=[-0.1  -0.01];
Kl34=[-0.1 -0.09];

Kl4=[7  1.2];
Kl41=[0.1  0.01];
Kl42=[-0.1  -0.01];
Kl43=[0.1 0.01];


K1=[-27 -4];
K2=[17 3];
K3=[-18 -3];
K4=[11 2];



eKl1=[];
eKl12=[];
eKl13=[];
eKl14=[];
eKl2=[];
eKl21=[];
eKl23=[];
eKl24=[];
eKl3=[];
eKl31=[];
eKl32=[];
eKl34=[];
eKl4=[];
eKl41=[];
eKl42=[];
eKl43=[];
dKl1=0.5;
dKl2=0.5;
dKl3=0.5;
dKl4=0.5;

KK1=[];
KK2=[];
KK3=[];
KK4=[];



kk=0;
T=0.01;

iii=0;

itlr=80;


tic
for k=1:itlr
    
      Kl1last=Kl1;
      Kl2last=Kl2;
      Kl3last=Kl3;
      Kl4last=Kl4;


    e1=0.00111*rand(1);
    e2=0.00111*rand(1);
    e3=0.00111*rand(1);
    e4=0.00111*rand(1);


    
     u1=-K1*x(:,k)+e1;
     u2=-K2*x(:,k)+e2;
     u3=-K3*x(:,k)+e3;
     u4=-K4*x(:,k)+e4;

    tspanxd=[0 T];
    X(:,k)=[x(:,k)',zeros(1,36)]';
    [tl,dl]= ode45(@caseA_fode,tspanxd,X(:,k));
    x(:,k+1)=[dl(length(tl),1:2)]; 
    
    

    iii=iii+1
    dxx(iii,:)=[x(1,k+1)^2,2*x(1,k+1)*x(2,k+1),x(2,k+1)^2]-[x(1,k)^2,2*x(1,k)*x(2,k),x(2,k)^2];
              
     lk1(iii,:)=dl(length(tl),3:4);
    lk12(iii,:)=dl(length(tl),5:6);
    lk13(iii,:)=dl(length(tl),7:8);
    lk14(iii,:)=dl(length(tl),9:10);
     lk2(iii,:)=dl(length(tl),11:12);
    lk21(iii,:)=dl(length(tl),13:14);
    lk23(iii,:)=dl(length(tl),15:16);
    lk24(iii,:)=dl(length(tl),17:18);
     lk3(iii,:)=dl(length(tl),19:20);
    lk31(iii,:)=dl(length(tl),21:22);
    lk32(iii,:)=dl(length(tl),23:24);
    lk34(iii,:)=dl(length(tl),25:26);
     lk4(iii,:)=dl(length(tl),27:28);
    lk41(iii,:)=dl(length(tl),29:30);
    lk42(iii,:)=dl(length(tl),31:32);
    lk43(iii,:)=dl(length(tl),33:34);
    
    fai1(iii,:)=[dxx(iii,:),lk1(iii,:),lk12(iii,:),lk13(iii,:),lk14(iii,:)];
    fai2(iii,:)=[dxx(iii,:),lk2(iii,:),lk21(iii,:),lk23(iii,:),lk24(iii,:)];
    fai3(iii,:)=[dxx(iii,:),lk3(iii,:),lk31(iii,:),lk32(iii,:),lk34(iii,:)];
    fai4(iii,:)=[dxx(iii,:),lk4(iii,:),lk41(iii,:),lk42(iii,:),lk43(iii,:)];
    
    r1(iii)=-1*dl(length(tl),35);
    r2(iii)=-1*dl(length(tl),36);  
    r3(iii)=-1*dl(length(tl),37);
    r4(iii)=-1*dl(length(tl),38);
    
    a1=rank(fai1)
    a2=rank(fai2)
    a3=rank(fai3);
    a4=rank(fai4);

  if dKl1(end)>0.0001 || dKl2(end)>0.0001 || dKl3(end)>0.0001  || dKl4(end)>0.0001 
    
     if a1==11 && a2==11 && a3==11 && a4==11
        
                kk=kk+1;
          
                 
           

         p1=fai1\r1';
         p2=fai2\r2';
         p3=fai3\r3';
         p4=fai4\r4';
         
         P1=[p1(1) p1(2)
             p1(2) p1(3)];
         Kl1=[ p1(4)  p1(5)];
         Kl12=[ p1(6)  p1(7)];
         Kl13=[ p1(8)  p1(9)];
         Kl14=[ p1(10)  p1(11)];
         
         KK1=[KK1; Kl1];
         
         P2=[p2(1) p2(2)
             p2(2) p2(3)];
         Kl2=[ p2(4)  p2(5)];
         Kl21=[ p2(6)  p2(7)];
         Kl23=[ p2(8)  p2(9)];
         Kl24=[ p2(10)  p2(11)];
         
         KK2=[KK2; Kl2];
         
         
         P3=[p3(1) p3(2)
             p3(2) p3(3)];
         Kl3=[ p3(4)  p3(5)];
         Kl31=[ p3(6)  p3(7)];
         Kl32=[ p3(8)  p3(9)];
         Kl34=[ p3(10)  p3(11)];
         
         KK3=[KK3; Kl3];
         
         P4=[p4(1) p4(2)
             p4(2) p4(3)];
         Kl4=[ p4(4)  p4(5)];
         Kl41=[ p4(6)  p4(7)];
         Kl42=[ p4(8)  p4(9)];
         Kl43=[ p4(10)  p4(11)];
         
   KK4=[KK4; Kl4];
         
         eKl1(kk)=norm(Kl1-K1e);
        eKl12(kk)=norm(Kl12-K12e);
        eKl13(kk)=norm(Kl13-K13e);
        eKl14(kk)=norm(Kl14-K14e);
         eKl2(kk)=norm(Kl2-K2e);
        eKl21(kk)=norm(Kl21-K21e);
        eKl23(kk)=norm(Kl23-K23e);
        eKl24(kk)=norm(Kl24-K24e);
         eKl3(kk)=norm(Kl3-K3e);
        eKl31(kk)=norm(Kl31-K31e);
        eKl32(kk)=norm(Kl32-K32e);
        eKl34(kk)=norm(Kl34-K34e);
         eKl4(kk)=norm(Kl4-K4e);
        eKl41(kk)=norm(Kl41-K41e);
        eKl42(kk)=norm(Kl42-K42e);
        eKl43(kk)=norm(Kl43-K43e);
         eP1(kk)=norm(P1-P1e);
         eP2(kk)=norm(P2-P2e);
         eP3(kk)=norm(P3-P3e);
         eP4(kk)=norm(P4-P4e);
        
        
      dKl1(kk)=norm(Kl1-Kl1last);
      dKl2(kk)=norm(Kl2-Kl2last);
      dKl3(kk)=norm(Kl3-Kl3last);
      dKl4(kk)=norm(Kl4-Kl4last);
        
         a1=0;
         a2=0;
         a3=0;
         a4=0;         
         lk1=[];
         lk12=[];
         lk13=[];
         lk14=[];       
         lk2=[];
         lk21=[];
         lk23=[];
         lk24=[];
          lk3=[];
         lk31=[];
         lk32=[];
         lk34=[];  
         lk4=[];
         lk41=[];
         lk42=[];
         lk43=[]; 
         dxx=[];  
         fai1=[];
         fai2=[];
         fai3=[];
         fai4=[];
         r1=[];
         r2=[];
         r3=[];
         r4=[];
         p1=[];
         p2=[];
         p3=[];
         p4=[];
         
         iii=0;
            
        
         
     end
        
   end

        
     if dKl1(end)<=0.0001 &&  dKl2(end)<=0.0001 && dKl3(end)<=0.0001 &&  dKl4(end)<=0.0001
         
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
hold on
plot(eP4,'-d','LineWidth',2);
xlabel('Time iteration $h$');legend(' $\Vert P_1^h-P_1 \Vert$',' $\Vert P_2^h-P_2 \Vert$',' $\Vert P_3^h-P_3 \Vert$',' $\Vert P_4^h-P_4 \Vert$');
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
hold on
plot(eKl4,'-d','LineWidth',2);
xlabel('Time iteration $h$');legend(' $\Vert K_1^h-K_{1}^* \Vert$',' $\Vert K_{2}^h-K_{2}^* \Vert$',' $\Vert K_{3}^h-K_{3}^* \Vert$',' $\Vert K_{4}^h-K_{4}^* \Vert$');
set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'XLabel'),'Interpreter','latex');
set(get(gca,'Legend'),'FontSize',26);
set(legend,'Interpreter','latex');



     figure  
plot(eKl12,'-d','LineWidth',2);
hold on
plot(eKl13,'-d','LineWidth',2);
hold on
plot(eKl14,'-d','LineWidth',2);
xlabel('Time iteration $h$');legend(' $\Vert K_{12}^h-K_{12}^* \Vert$',' $\Vert K_{13}^h-K_{13}^* \Vert$',' $\Vert K_{14}^h-K_{14}^* \Vert$');
set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'XLabel'),'Interpreter','latex');
set(get(gca,'Legend'),'FontSize',26);
set(legend,'Interpreter','latex');


     figure  
plot(eKl21,'-d','LineWidth',2);
hold on
plot(eKl23,'-d','LineWidth',2);
hold on
plot(eKl24,'-d','LineWidth',2);
xlabel('Time iteration $h$');legend(' $\Vert K_{21}^h-K_{21}^* \Vert$',' $\Vert K_{23}^h-K_{23}^* \Vert$',' $\Vert K_{24}^h-K_{24}^* \Vert$');
set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'XLabel'),'Interpreter','latex');
set(get(gca,'Legend'),'FontSize',26);
set(legend,'Interpreter','latex');

       figure  
plot(eKl31,'-d','LineWidth',2);
hold on
plot(eKl32,'-d','LineWidth',2);
hold on
plot(eKl34,'-d','LineWidth',2);
xlabel('Time iteration $h$');legend('$\Vert K_{31}^h-K_{31}^* \Vert$','$\Vert K_{32}^h-K_{32}^* \Vert$','$\Vert K_{34}^h-K_{34}^* \Vert$');
set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'XLabel'),'Interpreter','latex');
set(get(gca,'Legend'),'FontSize',26);
set(legend,'Interpreter','latex');


       figure  
plot(eKl41,'-d','LineWidth',2);
hold on
plot(eKl42,'-d','LineWidth',2);
hold on
plot(eKl43,'-d','LineWidth',2);
xlabel('Time iteration $h$');legend(' $\Vert K_{41}^h-K_{41}^* \Vert$',' $\Vert K_{42}^h-K_{42}^* \Vert$',' $\Vert K_{43}^h-K_{43}^* \Vert$');
set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'XLabel'),'Interpreter','latex');
set(get(gca,'Legend'),'FontSize',26);
set(legend,'Interpreter','latex');



% uu1=[];
% uu2=[];
% uu3=[];
% uu4=[];
% 
% xx = [rand(1,2)]';
% 
% itlrr=100;
%  for l=1:itlrr
% 
%     tspand=[0*T,1*T];
% 
%  if l<=kk
%             
%     u1 = -KK1(l,:)*xx(:,l);
%     u2 = -KK2(l,:)*xx(:,l);
%     u3 = -KK3(l,:)*xx(:,l);
%     u4 = -KK4(l,:)*xx(:,l);
%     
%  else
%      
%     u1 = -KK1(end,:)*xx(:,l);
%     u2 = -KK2(end,:)*xx(:,l);
%     u3 = -KK3(end,:)*xx(:,l);
%     u4 = -KK4(end,:)*xx(:,l);
%      
%  end   
%     
%     uu1=[uu1 u1];
%     uu2=[uu2 u2];
%     uu3=[uu3 u3];
%     uu4=[uu4 u4];
%   
%      
%     [t,d]= ode23('dynamic',tspand,xx(:,l));
%     xx(:,l+1)=[d(length(t),1:2)]; 
%     
%  end
% 
% 
% 
% 
% 
% 
% figure
%   t=0:T:itlrr*T;
% subplot(2,1,1); 
% plot(t,xx(1,:),'r','Linewidth',2);
% legend('$x_1$');
% set(get(gca,'XLabel'),'FontSize',24);
% set(get(gca,'Legend'),'FontSize',26);
% set(legend,'Interpreter','latex');
% subplot(2,1,2);
% plot(t,xx(2,:),'Linewidth',2);
% xlabel('Time (s)');
% legend('$x_2$');
% set(get(gca,'XLabel'),'FontSize',24);
% set(get(gca,'Legend'),'FontSize',26);
% set(legend,'Interpreter','latex');
% 
% 
% 
% 
% figure
%   tt=T:T:itlrr*T;
%   subplot(4,1,1); 
%   plot(tt,uu1,'-.','Linewidth',2);
%   legend('$u_1$');
% set(get(gca,'XLabel'),'FontSize',24);
% set(get(gca,'XLabel'),'Interpreter','latex');
% set(get(gca,'Legend'),'FontSize',26);
% set(legend,'Interpreter','latex');
% % title('Input vs Time');
% subplot(4,1,2); 
%  plot(tt,uu2,'-.','Linewidth',2);
%  legend('$u_2$');
% set(get(gca,'XLabel'),'FontSize',24);
% set(get(gca,'XLabel'),'Interpreter','latex');
% set(get(gca,'Legend'),'FontSize',26);
% set(legend,'Interpreter','latex');
% subplot(4,1,3); 
%  plot(tt,uu3,'-.','Linewidth',2);
%  legend('$u_3$');
% set(get(gca,'XLabel'),'FontSize',24);
% set(get(gca,'XLabel'),'Interpreter','latex');
% set(get(gca,'Legend'),'FontSize',26);
% set(legend,'Interpreter','latex');
% subplot(4,1,4); 
%  plot(tt,uu4,'-.','Linewidth',2);
% xlabel('Time (s)');
% legend('$u_4$');
% set(get(gca,'XLabel'),'FontSize',24);
% set(get(gca,'XLabel'),'Interpreter','latex');
% set(get(gca,'Legend'),'FontSize',26);
% set(legend,'Interpreter','latex');
% 
% 
% 
% 
