close all;
clc;
clear all;



 global A B1 B2 B3 B4 ;


A=[0         1;
    31.5397  0];

B1=[0;-4];
B2=[1;-1];
B3=[-1;1.7];
B4=[1;0];



%% expert parameters
  



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
       
       
       
       
 
 %eig(A-B1*K1e-B2*K2e-B3*K3e-B4*K4e)
 


 K1=[-20  -2.1];
K12=[0  0];
K13=[0  0];
K14=[0  0];



K2=[15  2];
K21=[0  0];
K23=[0  0];
K24=[0  0]; 

K3=[-12  -2];
K31=[0  0];
K32=[0  0];
K34=[0  0];

K4=[7  1.2];
K41=[0  0];
K42=[0  0];
K43=[0  0];


 
 
 iter=10;


tic
for i=1:iter
   
    
    %% use expert information to update kernel matrix Q1, Q2
     
       Abar1(:,:,i)=A-B1*K1(:,:,i)-B2*K12(:,:,i)-B3*K13(:,:,i)-B4*K14(:,:,i);
       S1(:,:,i)=K1(:,:,i)'*R11*K1(:,:,i)-gamma1^2*K12(:,:,i)'*R12*K12(:,:,i)+Q1-gamma1^2*K13(:,:,i)'*R13*K13(:,:,i)-gamma1^2*K14(:,:,i)'*R14*K14(:,:,i);
       P1(:,:,i)=lyap(Abar1(:,:,i)',S1(:,:,i));
       
       Abar2(:,:,i)=A-B2*K2(:,:,i)-B1*K21(:,:,i)-B3*K23(:,:,i)-B4*K24(:,:,i);
       S2(:,:,i)=K2(:,:,i)'*R22*K2(:,:,i)-gamma2^2*K21(:,:,i)'*R21*K21(:,:,i)+Q2-gamma2^2*K23(:,:,i)'*R23*K23(:,:,i)-gamma2^2*K24(:,:,i)'*R24*K24(:,:,i);
       P2(:,:,i)=lyap(Abar2(:,:,i)',S2(:,:,i));       
       
       
       Abar3(:,:,i)=A-B3*K3(:,:,i)-B1*K31(:,:,i)-B2*K32(:,:,i)-B4*K34(:,:,i);
       S3(:,:,i)=K3(:,:,i)'*R33*K3(:,:,i)-gamma3^2*K31(:,:,i)'*R31*K31(:,:,i)+Q3-gamma3^2*K32(:,:,i)'*R32*K32(:,:,i)-gamma3^2*K34(:,:,i)'*R34*K34(:,:,i);
       P3(:,:,i)=lyap(Abar3(:,:,i)',S3(:,:,i));       
       
       
       Abar4(:,:,i)=A-B4*K4(:,:,i)-B1*K41(:,:,i)-B2*K42(:,:,i)-B3*K43(:,:,i);
       S4(:,:,i)=K4(:,:,i)'*R44*K4(:,:,i)-gamma4^2*K41(:,:,i)'*R41*K41(:,:,i)+Q4-gamma4^2*K42(:,:,i)'*R42*K42(:,:,i)-gamma4^2*K43(:,:,i)'*R43*K43(:,:,i);
       P4(:,:,i)=lyap(Abar4(:,:,i)',S4(:,:,i));       
       
                 
       K1(:,:,i+1)=inv(R11)*B1'*P1(:,:,i);
       K12(:,:,i+1)=-1/gamma1^2*inv(R12)*B2'*P1(:,:,i);
       K13(:,:,i+1)=-1/gamma1^2*inv(R13)*B3'*P1(:,:,i);
       K14(:,:,i+1)=-1/gamma1^2*inv(R14)*B4'*P1(:,:,i);

       K2(:,:,i+1)=inv(R22)*B2'*P2(:,:,i);
       K21(:,:,i+1)=-1/gamma2^2*inv(R21)*B1'*P2(:,:,i);
       K23(:,:,i+1)=-1/gamma2^2*inv(R23)*B3'*P2(:,:,i);
       K24(:,:,i+1)=-1/gamma2^2*inv(R24)*B4'*P2(:,:,i);
  
       K3(:,:,i+1)=inv(R33)*B3'*P3(:,:,i);
       K31(:,:,i+1)=-1/gamma3^2*inv(R31)*B1'*P3(:,:,i);
       K32(:,:,i+1)=-1/gamma3^2*inv(R32)*B2'*P3(:,:,i);
       K34(:,:,i+1)=-1/gamma3^2*inv(R34)*B4'*P3(:,:,i);
 
       K4(:,:,i+1)=inv(R44)*B4'*P4(:,:,i);
       K41(:,:,i+1)=-1/gamma4^2*inv(R41)*B1'*P4(:,:,i);
       K42(:,:,i+1)=-1/gamma4^2*inv(R42)*B2'*P4(:,:,i);
       K43(:,:,i+1)=-1/gamma4^2*inv(R43)*B3'*P4(:,:,i);
       
%        p1(:,i)=norm(P1(:,:,i)-P1e);
%        p2(:,i)=norm(P2(:,:,i)-P2e);
%        p3(:,i)=norm(P3(:,:,i)-P3e);
%        p4(:,i)=norm(P4(:,:,i)-P4e);
%        
       k1(:,i)=norm(K1(:,:,i)-K1e);
       k12(:,i)=norm(K12(:,:,i)-K12e);
       k13(:,i)=norm(K13(:,:,i)-K13e);
       k14(:,i)=norm(K14(:,:,i)-K14e);
       
       k2(:,i)=norm(K2(:,:,i)-K2e);
       k21(:,i)=norm(K21(:,:,i)-K21e);
       k23(:,i)=norm(K23(:,:,i)-K23e);
       k24(:,i)=norm(K24(:,:,i)-K24e);
       
       k3(:,i)=norm(K3(:,:,i)-K3e);
       k31(:,i)=norm(K31(:,:,i)-K31e);
       k32(:,i)=norm(K32(:,:,i)-K32e);
       k34(:,i)=norm(K34(:,:,i)-K34e);
       
       k4(:,i)=norm(K4(:,:,i)-K4e);
       k41(:,i)=norm(K41(:,:,i)-K41e);
       k42(:,i)=norm(K42(:,:,i)-K42e);
       k43(:,i)=norm(K43(:,:,i)-K43e);

       


       p1(:,i)=norm(P1(:,:,i)-P1e);
       p2(:,i)=norm(P2(:,:,i)-P2e);
       p3(:,i)=norm(P3(:,:,i)-P3e);
       p4(:,i)=norm(P4(:,:,i)-P4e);      
            
%        k1(:,i)=norm(K1(:,:,i));
%        k12(:,i)=norm(K12(:,:,i));
%        k13(:,i)=norm(K13(:,:,i));
%        k14(:,i)=norm(K14(:,:,i));
%        
%        k2(:,i)=norm(K2(:,:,i));
%        k21(:,i)=norm(K21(:,:,i));
%        k23(:,i)=norm(K23(:,:,i));
%        k24(:,i)=norm(K24(:,:,i));
%        
%        k3(:,i)=norm(K3(:,:,i));
%        k31(:,i)=norm(K31(:,:,i));
%        k32(:,i)=norm(K32(:,:,i));
%        k34(:,i)=norm(K34(:,:,i));
%        
%        k4(:,i)=norm(K4(:,:,i));
%        k41(:,i)=norm(K41(:,:,i));
%        k42(:,i)=norm(K42(:,:,i));
%        k43(:,i)=norm(K43(:,:,i));
       
       %barA(:,i)=eig(A-B1*K1(:,:,i)-B2*K2(:,:,i)-B3*K3(:,:,i)-B4*K4(:,:,i));
       
       
       
       
       
       
              
end 
       


toc






     figure  
plot(p1,'-*','LineWidth',2);
hold on
plot(p2,'-*','LineWidth',2);
hold on
plot(p3,'-*','LineWidth',2);
hold on
plot(p4,'-*','LineWidth',2);
xlabel('Time iteration $h$');legend(' $\Vert P_1^h-P_1 \Vert$',' $\Vert P_2^h-P_2 \Vert$',' $\Vert P_3^h-P_3 \Vert$',' $\Vert P_4^h-P_4 \Vert$');
set(get(gca,'XLabel'),'FontSize',18);
set(get(gca,'XLabel'),'Interpreter','latex');
set(get(gca,'Legend'),'FontSize',18);
set(legend,'Interpreter','latex');

   


     figure  
plot(k1,'-*','LineWidth',2);
hold on
plot(k2,'-*','LineWidth',2);
hold on
plot(k3,'-*','LineWidth',2);
hold on
plot(k4,'-*','LineWidth',2);
xlabel('Time iteration $h$');legend(' $\Vert K_1^h-K_{1}^* \Vert$',' $\Vert K_{2}^h-K_{2}^* \Vert$',' $\Vert K_{3}^h-K_{3}^* \Vert$',' $\Vert K_{4}^h-K_{4}^* \Vert$');
set(get(gca,'XLabel'),'FontSize',18);
set(get(gca,'XLabel'),'Interpreter','latex');
set(get(gca,'Legend'),'FontSize',18);
set(legend,'Interpreter','latex');


   
     figure  
plot(k12,'-*','LineWidth',2);
hold on
plot(k13,'-*','LineWidth',2);
hold on
plot(k14,'-*','LineWidth',2);
xlabel('Time iteration $h$');legend(' $\Vert K_{12}^h-K_{12}^* \Vert$',' $\Vert K_{13}^h-K_{13}^* \Vert$',' $\Vert K_{14}^h-K_{14}^* \Vert$');
set(get(gca,'XLabel'),'FontSize',18);
set(get(gca,'XLabel'),'Interpreter','latex');
set(get(gca,'Legend'),'FontSize',18);
set(legend,'Interpreter','latex');

   
     figure  
plot(k21,'-*','LineWidth',2);
hold on
plot(k23,'-*','LineWidth',2);
hold on
plot(k24,'-*','LineWidth',2);
xlabel('Time iteration $h$');legend(' $\Vert K_{21}^h-K_{21}^* \Vert$',' $\Vert K_{23}^h-K_{23}^* \Vert$',' $\Vert K_{24}^h-K_{24}^* \Vert$');
set(get(gca,'XLabel'),'FontSize',18);
set(get(gca,'XLabel'),'Interpreter','latex');
set(get(gca,'Legend'),'FontSize',18);
set(legend,'Interpreter','latex');

   
       figure  
plot(k31,'-*','LineWidth',2);
hold on
plot(k32,'-*','LineWidth',2);
hold on
plot(k34,'-*','LineWidth',2);
xlabel('Time iteration $h$');legend('$\Vert K_{31}^h-K_{31}^* \Vert$','$\Vert K_{32}^h-K_{32}^* \Vert$','$\Vert K_{34}^h-K_{34}^* \Vert$');
set(get(gca,'XLabel'),'FontSize',18);
set(get(gca,'XLabel'),'Interpreter','latex');
set(get(gca,'Legend'),'FontSize',18);
set(legend,'Interpreter','latex');


       figure  
plot(k41,'-*','LineWidth',2);
hold on
plot(k42,'-*','LineWidth',2);
hold on
plot(k43,'-*','LineWidth',2);
xlabel('Time iteration $h$');legend(' $\Vert K_{41}^h-K_{41}^* \Vert$',' $\Vert K_{42}^h-K_{42}^* \Vert$',' $\Vert K_{43}^h-K_{43}^* \Vert$');
set(get(gca,'XLabel'),'FontSize',18);
set(get(gca,'XLabel'),'Interpreter','latex');
set(get(gca,'Legend'),'FontSize',18);
set(legend,'Interpreter','latex');










% uu1=[];
% uu2=[];
% uu3=[];
% uu4=[];
% 
% xx = rand(1,2)';
% iter=80;
% T=0.01;
%  for l=1:iter
% 
%     tspand=[0*T,1*T];
% 
%  
%             
%     u1 = -K1(:,:,end)*xx(:,l);
%     u2 = -K2(:,:,end)*xx(:,l);
%     u3 = -K3(:,:,end)*xx(:,l);% xx = rand(1,2)';
% iter=80;
% T=0.01;
%  for l=1:iter
% 
%     tspand=[0*T,1*T];
% 
%  
%             
%     u1 = -K1(:,:,end)*xx(:,l);
%     u2 = -K2(:,:,end)*xx(:,l);
%     u3 = -K3(:,:,end)*xx(:,l);
%     u4 = -K4(:,:,end)*xx(:,l);
%     
%     uu1=[uu1 u1];
%     uu2=[uu2 u2];
%     uu3=[uu3 u3];
%     uu4=[uu4 u4];
%   
%      
%     [t,d]= ode45('dynamic',tspand,xx(:,l));
%     xx(:,l+1)=[d(length(t),1:2)]; 
%     
%  end
% 
% figure
% plot(xx(1,:),'Linewidth',2)
% hold on
% plot(xx(2,:),'Linewidth',2)
% legend('x1','x2')

%     u4 = -K4(:,:,end)*xx(:,l);
%     
%     uu1=[uu1 u1];
%     uu2=[uu2 u2];
%     uu3=[uu3 u3];
%     uu4=[uu4 u4];
%   
%      
%     [t,d]= ode45('dynamic',tspand,xx(:,l));
%     xx(:,l+1)=[d(length(t),1:2)]; 
%     
%  end
% 
% figure
% plot(xx(1,:),'Linewidth',2)
% hold on
% plot(xx(2,:),'Linewidth',2)
% legend('x1','x2')



