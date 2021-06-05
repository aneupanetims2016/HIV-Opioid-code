

clear all
close all
clc

a6 = 1;
a7 = 1;
a8 = 1;
a9 = 1;


initial_cond = [7779681  42450 79912 10 5 31384];
terminal_cond = [ 0 0 0 0 0 0];
tforward= 0:0.04:40;
tbackwards = 40:-0.04:0;
options = odeset('RelTol',1e-9,'AbsTol',1e-12);

params = [0.374770124132209,0.715206465517946,0.197263367461474,1.00388980628920,0.210020488926058];


% S0N0ratio  = initial_cond (1)/(initial_cond(1)+initial_cond(2)+initial_cond(3)+initial_cond(4)+initial_cond(5)+initial_cond(6));
% beta = 4.3/(S0N0ratio *(params (1)/(1/14) + (1-0.42)*params (4)/1 + 0.42/(params(5) + 0.1765*params(5))));


mu= 1/79 ;

Lambda= 7779681*mu ;

   
beta_u= params(1) ;
beta_v= params(2) ;
delta=0.0010;
q_u=10.0454;
gamma_v=0.499993;
delta_2= 0.9987;
mu_u= params(3);
mu_v=params(5);
q_v=params(4);
alpha= 24.9993;
gamma_2= 0.99921;
mu_a=0.0879;
q  = 0.03;
mu_1 = 46.1262;
mu_2 = 8.8916;



R1 = beta_u/(mu+mu_u+delta)

R2 = beta_v/(mu+mu_v+gamma_v)


u = zeros(length(tforward),4);
u2 = ones(length(tforward),4);
[t_H, y_H] = ode15s(@(t,y)Model_HIV_Control(t,y,tforward,u),tforward,initial_cond,options);
[t_O, y_O] = ode15s(@(t,y)Model_HIV_Control(t,y,tforward,u2),tforward,initial_cond,options);

figure
plot(t_H,y_H(:,1),'-.','LineWidth',2.5)
 hold on
plot(t_O,y_O(:,1),'g','LineWidth',2.5)
legend('with u=0 ','with u=1') 
 title('Number of susceptible Individuals','fontweight','normal','fontsize',18)
 ylabel('U(t)','fontweight','normal','fontsize',18)
 xlabel('Years','fontweight','normal','fontsize',18)
 xticks([1 6 11 16 21 26 31 36 40 ])
 yticks([2000000 4000000 6000000 8000000 ])
 yticklabels({'2x10^6', '4x10^6', '6x10^6' ,'8x10^6'})
 set(gca, 'YGrid', 'on', 'XGrid', 'off','LineWidth', 2,'fontsize',14)
hold off


figure
plot(t_H,y_H(:,2),'-.','LineWidth',2.5)
 hold on
plot(t_O,y_O(:,2),'g','LineWidth',2.5)
legend('with u=0 ','with u=1','Location','northwest') 
 title('Number of opioid addicted Individuals','fontweight','normal','fontsize',18)
 ylabel('U(t)','fontweight','normal','fontsize',18)
 xlabel('Years','fontweight','normal','fontsize',18)
 xticks([1 6 11 16 21 26 31 36 40 ])
 yticks([200000 400000 600000 800000 1000000 1200000 1400000])
 yticklabels({'2x10^5', '4x10^5', '6x10^5' ,'8x10^5' ,'10x10^5','12x10^5','14x10^5'})
 set(gca, 'YGrid', 'on', 'XGrid', 'off','LineWidth', 2,'fontsize',14)
hold off

figure
plot(t_H,y_H(:,3),'-.','LineWidth',2.5)
 hold on
plot(t_O,y_O(:,3),'g','LineWidth',2.5)
legend('with u=0','with u=1') 
 title('Number of HIV Positive','fontweight','normal','fontsize',18)
 ylabel('V(t)','fontweight','normal','fontsize',18)
 xlabel('Years','fontweight','normal','fontsize',18)
 xticks([1 6 11 16 21 26 31 36 40 ])
 yticks([500000 1000000 1500000 2000000 2500000 3000000 3500000 ])
 yticklabels({'0.5x10^6', '1x10^6', '1.5x10^6' ,'2x10^6', '2.5x10^6','3x10^6','3.5x10^6'})
 set(gca, 'YGrid', 'on', 'XGrid', 'off','LineWidth', 2,'fontsize',14)
hold off

figure
plot(t_H,y_H(:,4),'-.','LineWidth',2.5)
 hold on
plot(t_O,y_O(:,4),'g','LineWidth',2.5)
legend('with u=0 ','with u=1') 
 title('Number of addicted and first stage infected','fontweight','normal','fontsize',18)
 ylabel('I_1(t)','fontweight','normal','fontsize',18)
 xlabel('Years','fontweight','normal','fontsize',18)
 xticks([1 6 11 16 21 26 31 36 40 ])
%    yticks([50 100 150 200 250 ])
%    yticklabels({'50', '100', '150' ,'200', '250'})
 set(gca, 'YGrid', 'on', 'XGrid', 'off','LineWidth', 2,'fontsize',14)
hold off

figure
plot(t_H,y_H(:,5),'-.','LineWidth',2.5)
 hold on
plot(t_O,y_O(:,5),'g','LineWidth',2.5)
legend('with u=0 ','with u=1') 
 title('Number of addicted and second stage infected','fontweight','normal','fontsize',18)
 ylabel('I_2(t)','fontweight','normal','fontsize',18)
 xlabel('Years','fontweight','normal','fontsize',18)
 xticks([1 6 11 16 21 26 31 36 40 ])
%     yticks([100 200 300 400 500 ])
%     yticklabels({'100', '200', '300' ,'400', '500'})
 set(gca, 'YGrid', 'on', 'XGrid', 'off','LineWidth', 2,'fontsize',14)
hold off

figure
plot(t_H,y_H(:,6),'-.','LineWidth',2.5)
 hold on
plot(t_O,y_O(:,6),'g','LineWidth',2.5)
legend('with u=0','with u=1') 
 title('Number of AIDS','fontweight','normal','fontsize',18)
 ylabel('A(t)','fontweight','normal','fontsize',18)
 xlabel('Years','fontweight','normal','fontsize',18)
 xticks([1 6 11 16 21 26 31 36 40 ])
 yticks([50000 100000 150000 200000 250000 ])
 yticklabels({'0.5x10^5', '10x10^5', '1.5x10^5' ,'2x10^5', '2.5x10^5'})
 set(gca, 'YGrid', 'on', 'XGrid', 'off','LineWidth', 2,'fontsize',14)
hold off

%     return

y = zeros(length(tforward),6);
lambda = zeros(length(tforward),6);
test = -1;
delta = 0.001;

 while(test <0)
   
ut = u;
yold = y;
lambda_old = lambda;
[t, y] = ode15s(@(t,y)Model_HIV_Control(t,y,tforward,ut),tforward,initial_cond,options);

[t_b, lambda] = ode15s(@(t,lambda)HIV_Adjoint_Eqn(t,lambda,tforward,y,ut),tbackwards,terminal_cond,options);


S = y(:,1);
U = y(:,2);
V = y(:,3);
I_1 = y(:,4);
I_2 = y(:,5);
A =  y(:,6);

N = (S+U+V+I_1+I_2+A);


tempu1= (delta.*(lambda(:,2)-lambda(:,1)).*U)./(2*a6);
tempu2=  q_u.*beta_v.*((lambda(:,4)-lambda(:,2)).*(V+I_1+I_2).*U)./(2*a7*N);
tempu3 = gamma_v.*(V.*(lambda(:,6)-lambda(:,3)))./(2*a8)+gamma_2.*(I_2.*(lambda(:,6)-lambda(:,5)))./(2*a8);
tempu4 = (delta_2.*((lambda(:,5)-lambda(:,3)).*I_2))./(2*a9);



u11 = min(10,max(0,tempu1));
u21 = min(1,max(0,tempu2));
u31 = min(1,max(0,tempu3));
u41 = min(1,max(0,tempu4));

u1 = 0.5*(u11 + ut(:,1));
u2 = 0.5*(u21 + ut(:,2));
u3 = 0.5*(u31 + ut(:,3));
u4 = 0.5*(u41 + ut(:,4));

u = [u1, u2,u3,u4];

  temp1 = delta*norm(u) - norm(u-ut);
  temp2 = delta*norm(y) - norm(y-yold);
  temp3 = delta*norm(lambda) - norm(lambda - lambda_old);
  
  test = min([temp1 temp2 temp3])
  
 end


[t_r, y_r] = ode15s(@(t,y)Model_HIV(y),tforward,initial_cond);



 
 
 figure 
 plot(t,y(:,1),'g','LineWidth',2.5)
 hold on
 plot(t,y_r(:,1),'--','LineWidth',2.5)
 legend('with controls','with out controls') 
 title('Number of Susceptible Individuals','fontweight','normal','fontsize',18)
 ylabel('S (t)','fontweight','normal','fontsize',18)
 xlabel('Years','fontweight','normal','fontsize',18)
 yticks([2000000 4000000 6000000 8000000])
 yticklabels({'2x10^6', '4x10^6', '6x10^6' ,'8x10^6'})
 xticks([1 6 11 16 21 26 31 36 40 ])
 xticklabels({'1','6','11','16','21','26','31','36','40'})
 set(gca, 'YGrid', 'on', 'XGrid', 'off','LineWidth', 2,'fontsize',14)
hold off 
figure
plot(t,y(:,2),'g','LineWidth',2.5)
 hold on
plot(t,y_r(:,2),'--','LineWidth',2.5)
legend('with controls','with out controls','location','northwest') 
 title('Number of opioid addicted Individuals','fontweight','normal','fontsize',18)
 ylabel('U(t)','fontweight','normal','fontsize',18)
 xlabel('Years','fontweight','normal','fontsize',18)
 yticks([200000 400000 600000 800000 1000000 1200000 1400000])
 yticklabels({'2x10^5','4x10^5','6x10^5','8x10^5','10x10^5','14x10^5'})
 xticks([1 6 11 16 21 26 31 36 40 ])
 xticklabels({'1','6','11','16','21','26','31','36','40'})
 set(gca, 'YGrid', 'on', 'XGrid', 'off','LineWidth', 2,'fontsize',14)
hold off

figure
plot(t,y(:,3),'g','LineWidth',2.5)
 hold on
plot(t,y_r(:,3),'--','LineWidth',2.5)
legend('with controls','with out controls') 
 title('Number of HIV positive Individuals','fontweight','normal','fontsize',18)
 ylabel('V(t)','fontweight','normal','fontsize',18)
 xlabel('Years','fontweight','normal','fontsize',18)
%  yticks([200000  400000 600000  ])
%  yticklabels({'2x10^5', '4x10^5', '6x10^5'})
 xticks([1 6 11 16 21 26 31 36 40 ])
 xticklabels({'1','6','11','16','21','26','31','36','40'})
 set(gca, 'YGrid', 'on', 'XGrid', 'off','LineWidth', 2,'fontsize',14)
hold off
% 
figure
plot(t,y(:,4),'g','LineWidth',2.5)
 hold on
plot(t,y_r(:,4),'--','LineWidth',2.5)
legend('with controls','with out controls','location','northwest') 
 title('Number of addicted with first stage infected Individuals','fontweight','normal','fontsize',15)
 ylabel('I_1(t)','fontweight','normal','fontsize',18)
 xlabel('Years','fontweight','normal','fontsize',18)
%   yticks([50 100 150 200 250])
%    yticklabels({'50', '100', '150' ,'200' ,'250'})
 xticks([1 6 11 16 21 26 31 36 40 ])
 xticklabels({'1','6','11','16','21','26','31','36','40'})
 set(gca, 'YGrid', 'on', 'XGrid', 'off','LineWidth', 2,'fontsize',12)
hold off
% 
figure
plot(t,y(:,5),'g','LineWidth',2.5)
 hold on
plot(t,y_r(:,5),'--','LineWidth',2.5)
legend('with controls','with out controls','location','northwest') 
 title('Number of addicted and second stage infected Individuals','fontweight','normal','fontsize',15)
 ylabel('I_2(t)','fontweight','normal','fontsize',18)
 xlabel('Years','fontweight','normal','fontsize',18)
%  yticks([100 200 300 400 500])
%  yticklabels({'100', '200', '300' ,'400' ,'500'})
 xticks([1 6 11 16 21 26 31 36 40 ])
 xticklabels({'1','6','11','16','21','26','31','36','40'})
 set(gca, 'YGrid', 'on', 'XGrid', 'off','LineWidth', 2,'fontsize',12)
hold off

figure
plot(t,y(:,6),'g','LineWidth',2.5)
 hold on
plot(t,y_r(:,6),'--','LineWidth',2.5)
legend('with controls','with out controls') 
 title('Number of AIDS Individuals','fontweight','normal','fontsize',18)
 ylabel('A(t)','fontweight','normal','fontsize',18)
 xlabel('Years','fontweight','normal','fontsize',18)
 yticks([200000 400000 600000 800000 1000000])
 yticklabels({'2x10^{5}', '4x10^{5}', '6x10^{5}' ,'8x10^{5}' ,'10x10^{5}'})
 xticks([1 6 11 16 21 26 31 36 40 ])
 xticklabels({'1','6','11','16','21','26','31','36','40'})
set(gca, 'YGrid', 'on', 'XGrid', 'off','LineWidth', 2,'fontsize',14)
hold off

figure
plot(tforward,u(:,1),'r-.','LineWidth',2.5)
title('Treatment for addicted individuals (u_1)','fontweight','normal','fontsize',18)
xlabel('Time(year)','fontweight','normal','fontsize',18)
ylabel('u_1(t)','fontweight','normal','fontsize',18)
legend('0\leq u_1 \leq 1') 
%yticks([0.2 0.4 0.6 0.8 1])
xlabel('Years','fontweight','normal','fontsize',18)
% yticks([0.00002 0.00004 0.00006 0.00008 0.0001 0.00012 0.00014 0.00016 0.00018 0.00025 ])   
%  yticklabels({'0.2x10^{-4}','0.4x10^{-4}','0.6x10^{-4}','0.8x10^{-4}', '1x10^{-4}', '1.2x10^{-4}','1.4x10^{-4}','1.6x10^{-4}','1.8x10^{-4}', '2.5X10^{-4}'})
xlabel('Years','fontweight','normal','fontsize',18)
xticks([1 6 11 16 21 26 31 36 40 ])
xticklabels({'1','6','11','16','21','26','31','36','40'})
set(gca, 'YGrid', 'on', 'XGrid', 'off','LineWidth', 2,'fontsize',14)
hold off

 figure
 plot(tforward,u(:,2),'r-.','LineWidth',2.5)
 title('Education and awareness for HIV spreading (u_2)','fontweight','normal','fontsize',12)
xlabel('Time(year)','fontweight','normal','fontsize',18)
ylabel('u_2(t)','fontweight','normal','fontsize',18)
legend('0\leq u_2 \leq 1','location','northwest') 
%yticks([0.2 0.4 0.6 0.8 1])
xlabel('Years','fontweight','normal','fontsize',18)
 ylim([0 1])  
% yticklabels({'0','1'})
xlabel('Years','fontweight','normal','fontsize',18)
xticks([1 6 11 16 21 26 31 36 40 ])
xticklabels({'1','6','11','16','21','26','31','36','40'})
set(gca, 'YGrid', 'on', 'XGrid', 'off')
hold off

 figure
 plot(tforward,u(:,3),'r-.','LineWidth',2.5)
title('HAART for HIV positive(u_3)','fontweight','normal','fontsize',18)
legend('0\leq u_3 \leq 1') 
xlabel('Time(year)','fontweight','normal','fontsize',18)
ylabel('u_3(t)','fontweight','normal','fontsize',18)
xlabel('Years','fontweight','normal','fontsize',18)
yticks([0.2 0.4 0.6 0.8 1 ])
yticklabels({'0.2','0.4', '0.6','0.8', '1'})
xlabel('Years','fontweight','normal','fontsize',18)
xticks([1 6 11 16 21 26 31 36 40 ])
xticklabels({'1','6','11','16','21','26','31','36','40'})
set(gca, 'YGrid', 'on', 'XGrid', 'off')
hold off

 figure
 plot(tforward,u(:,4),'r-.','LineWidth',2.5)
 title('Treatment for addiction for Infected-addicted individuals','fontweight','normal','fontsize',14)
xlabel('Time(year)','fontweight','normal','fontsize',18)
ylabel('u_4(t)','fontweight','normal','fontsize',18)
legend('0\leq u_4 \leq 1') 
 yticks([0.2 0.4 0.6 0.8 1  ])
 yticklabels({'0.2','0.4', '0.6','0.8', '1'})
xlabel('Years','fontweight','normal','fontsize',18)
xticks([1 6 11 16 21 26 31 36 40 ])
xticklabels({'1','6','11','16','21','26','31','36','40'})
set(gca, 'YGrid', 'on', 'XGrid', 'off')
hold off


function dy = Model_HIV_Control(t,y,timespan,utimedep)


params = [0.374770124132209,0.715206465517946,0.197263367461474,1.00388980628920,0.210020488926058];





mu= 1/79 ;

Lambda= 7779681*mu ;

   
beta_u= params(1) ;
beta_v= params(2) ;
delta=0.0010;
q_u=10.0454;
gamma_v=0.499993;
delta_2= 0.9987;
mu_u= params(3);
mu_v=params(5);
q_v=params(4);
alpha= 24.9993;
gamma_2= 0.99921;
mu_a=0.0879;
q  = 0.03;
mu_1 = 46.1262;
mu_2 = 8.8916;

dy = zeros(6,1);

% controls 

u1timedep = utimedep(:,1);
u2timedep = utimedep(:,2);
u3timedep = utimedep(:,3);
u4timedep = utimedep(:,4);

   u1t = interp1(timespan,u1timedep,t);
   u2t = interp1(timespan,u2timedep,t);
   u3t = interp1(timespan,u3timedep,t);
   u4t = interp1(timespan,u4timedep,t);
   
   
   

S  = y(1);
U  = y(2);
V = y(3);
I_1 = y(4);
I_2  = y(5);
A  = y(6);

N  = y(1) + y(2) + y(3) + y(4) + y(6);


dy(1) = Lambda- (beta_u*(U+I_1+I_2)*S./N) - (beta_v*(V+I_1+I_2)*S./N )- mu*S + delta*u1t*U ;
dy(2) = (beta_u*(U+I_1+I_2)*S./N)-(q_u*(1-u2t)*beta_v*(V+I_1+I_2)*U./N)- (mu+mu_u + delta*u1t)*U;
dy(3) = (beta_v*(V+I_1+I_2)*S./N)-(q_v*beta_u*(U+I_1+I_2)*V./N)- (mu+mu_v+gamma_v*(1-u3t))*V+ delta_2*u4t*I_2;
dy(4) =  (q_u*(1-u2t)*beta_v*(V+I_1+I_2)*U./N) + (q_v*beta_u*(U+I_1+I_2)*V./N)-(mu+alpha+mu_1)*I_1;
dy(5) = alpha*I_1-(mu+gamma_2*(1-u3t)+mu_2 + delta_2*u4t)*I_2;
dy(6) = (gamma_v*(1-u3t)*V) + (gamma_2*(1-u3t)*I_2)-((mu+mu_a)*A);

end

function dlambda = HIV_Adjoint_Eqn(t,lambda,timespan,ytimedep,utimedep)

params = [0.374770124132209,0.715206465517946,0.197263367461474,1.00388980628920,0.210020488926058];




mu= 1/79 ;
Lambda= 7779681*mu ;   
beta_u= params(1) ;
beta_v= params(2) ;
delta=0.0010;
q_u=10.0454;
gamma_v=0.499993;
delta_2= 0.9987;
mu_u= params(3);
mu_v=params(5);
q_v=params(4);
alpha= 24.9993;
gamma_2= 0.99921;
mu_a=0.0879;
q  = 0.03;
mu_1 = 46.1262;
mu_2 = 8.8916;

dlambda = zeros(6,1);

S = ytimedep(:,1);
U = ytimedep(:,2);
V = ytimedep(:,3);
I_1 = ytimedep(:,4);
I_2 = ytimedep(:,5);
A = ytimedep(:,6);


   St = interp1(timespan,S,t);
   Ut = interp1(timespan,U,t);
   Vt = interp1(timespan,V,t);
   I_1t = interp1(timespan,I_1,t);
   I_2t = interp1(timespan,I_2,t);
   At = interp1(timespan,A,t);
   
   Nt = (St+Ut+Vt+I_1t+I_2t+At);
  
% controls 

u1timedep = utimedep(:,1);
u2timedep = utimedep(:,2);
u3timedep = utimedep(:,3);
u4timedep = utimedep(:,4);

   u1t = interp1(timespan,u1timedep,t);
   u2t = interp1(timespan,u2timedep,t);
   u3t = interp1(timespan,u3timedep,t);
   u4t = interp1(timespan,u4timedep,t);
   
   
   dlambda(1) =((lambda(1)-lambda(2))*beta_u*(Ut+I_1t+I_2t)*(Nt-St)./Nt^2)+...
         ((lambda(1)-lambda(3))*beta_v*(Vt+I_1t+I_2t)*(Nt-St)./Nt^2)+...
         ((lambda(4)-lambda(2))*q_u*(1-u2t)*beta_v*(Vt+I_1t+I_2t)*Ut./Nt^2)+...
         ((lambda(4)-lambda(3))*q_v*beta_u*(Ut+I_1t+I_2t)*Vt./Nt^2)+...
         (lambda(1)*mu);   
                                  
  
   dlambda(2) = -1 + ((lambda(1)-lambda(2))*beta_u*St*(St+Vt+At)./Nt^2)+...
               ((lambda(3)-lambda(1))*beta_v*St*(Vt+I_1t+I_2t)./Nt^2)+...
               ((lambda(2)-lambda(4))*q_u*(1-u2t)*beta_v*(Vt+I_1t+I_2t)*(Nt-Ut)./Nt^2)+...
               ((lambda(3)-lambda(4))*q_v*beta_u*Vt*(St+Vt+At)./Nt^2)+...
               ((lambda(2)-lambda(1))*delta*u1t)+...
               (lambda(2)*(mu+mu_u));
          
   dlambda(3) = -0.01 + ((lambda(2)-lambda(1))*beta_u*St*(Ut+I_1t+I_2t)./Nt^2)+...
               ((lambda(1)-lambda(3))*beta_v*St*(St+Ut+At)./Nt^2)+...
               ((lambda(2)-lambda(4))*q_u*(1-u2t)*beta_v*Ut*(St+Ut+At)./Nt^2)+...
               ((lambda(3)-lambda(4))*q_v*beta_u*(Ut+I_1t+I_2t)*(Nt-Vt)./Nt^2)+...
               ((lambda(3)-lambda(6))*gamma_v*(1-u3t))+...
               (lambda(3)*(mu+mu_v));
           
   dlambda(4) = -1 + ((lambda(1)-lambda(2))* beta_u*St* (St+Vt+At)./Nt^2) +...
               ((lambda(1)-lambda(3))*beta_v*St*(St+Ut+At)./Nt^2)+...
               ((lambda(2)-lambda(4))*q_u*(1-u2t)*beta_v*Ut*(St+Ut+At)./Nt^2)+...
               ((lambda(3)-lambda(4))*q_v*beta_u*Vt* (St+Vt+At)./Nt^2)+...
               ((lambda(4)-lambda(5))*alpha)+...
               (lambda(4)*(mu + mu_1));
           
   dlambda(5) = -1 + ((lambda(1)-lambda(2))* beta_u*St* (St+Vt+At)./Nt^2) +...
              ((lambda(1)-lambda(3))*beta_v*St*(St+Ut+At)./Nt^2)+...
              ((lambda(2)-lambda(4))*q_u*(1-u2t)*beta_v*Ut*(St+Ut+At)./Nt^2)+...
              ((lambda(3)-lambda(4))*q_v*beta_u*Vt* (St+Vt+At)./Nt^2)+...
              ((lambda(5)-lambda(6))*gamma_2*(1-u3t))+...
              (lambda(5)*(mu+mu_2))+...
              ((lambda(5)-lambda(3))*delta_2*u4t);  
   dlambda(6) = -0.89 + ((lambda(2)-lambda(1))*beta_u*St*(Ut+I_1t+I_2t)./Nt^2)+...
               ((lambda(3)-lambda(1))*beta_v*St*(Vt+I_1t+I_2t)./Nt^2)+...
               ((lambda(4)-lambda(2))*q_u*(1-u2t)*beta_v*(Vt+I_1t+I_2t)*Ut./Nt^2)+...
               ((lambda(4)-lambda(3))*q_v*beta_u*(Ut+I_1t+I_2t)*Vt./Nt^2)+...
               (lambda(6)*(mu+mu_a));
end

function dy = Model_HIV(y)


params = [0.374770124132209,0.715206465517946,0.197263367461474,1.00388980628920,0.210020488926058];


mu= 1/79 ;
Lambda= 7779681*mu ;   
beta_u= params(1) ;
beta_v= params(2) ;
delta=0.0010;
q_u=10.0454;
gamma_v=0.499993;
delta_2= 0.9987;
mu_u= params(3);
mu_v=params(5);
q_v=params(4);
alpha= 24.9993;
gamma_2= 0.99921;
mu_a=0.0879;
q  = 0.03;
mu_1 = 46.1262;
mu_2 = 8.8916;

dy = zeros(6,1);

S  = y(1);
U  = y(2);
V = y(3);
I_1 = y(4);
I_2 = y(5);
A = y(6);

N  = y(1) + y(2) + y(3) + y(4) + y(5)+ y(6);



   dy(1) = Lambda - beta_u*(U+I_1+I_2)*S./N - beta_v*(V+I_1+I_2)*S./N- mu*S + delta*U;   
   dy(2) = beta_u*(U+I_1+I_2)*S./N - q_u*beta_v*(V+I_1+I_2)*U./N - (mu + mu_u + delta)*U;
   dy(3) = beta_v*(V+I_1+I_2)*S./N - q_v*beta_u*(U+I_1+I_2)*V./N - (mu + mu_v + gamma_v)*V + delta_2*I_2;
   dy(4) = q_u*beta_v*(V+I_1+I_2)*U./N + q_v*beta_u*(U+I_1+I_2)*V./N - (mu + alpha + mu_1)*I_1;
   dy(5) = alpha*I_1 - (mu + gamma_2+ mu_2 + delta_2)*I_2;
   dy(6) = gamma_v*V + gamma_2*I_2 - (mu + mu_a)*A;

end
