
% this is the code which we use for parameter

clear all
close all
clc

%t=0 is 1999

global AIDSCases HIVCases HIVDeaths DrugDeaths...
       tforward  tmeasure_AIDScases tmeasure_OPIOIDdeaths...
       tmeasure_HIVcases tmeasure_HIVdeaths  

AIDSCases = [38285;  36922; 36726; 37317; 36220; 34261; 32790; 31984;...
31384; 30187; 27401; 25620; 24684; 23656; 19313; 18590; 18375; 17749;17113]; % data starts in 2000


% AIDSDeaths = [18318; 18472; 18397; 18287; 17898; 17573; 17186; 16352;...
% 15435; 14923; 13754; 13393; 13087; 12992; 13067; 12804; 12983; 12824; 12146]; % data starts in 2000

HIVCases = [47247; 44716; 43020; 41265; 40512; 39232; 40005; 39817;...
39569; 38351; 37428]; % data starts in 2008


HIVDeaths = [18525; 18043; 16742; 16300; 16018; 15908; 16145; 15860; 16395;...
16358; 15483]; % data starts in 2008

DrugDeaths = [8050; 8407; 9496; 11920; 12940; 13756; 14918; 17545; 18516;...
19582; 20422; 21089; 22784; 23166; 25052; 28647; 33091; 42249; 47600; 46802]; %data starts 1999


dt = 0.1;
tdata_AIDScases = 1:1:19;

tdata_OPIOIDdeaths = 0:1:19;
tdata_HIVcases = 9:1:19;
tdata_HIVdeaths = 9:1:19;

tforward = (0:dt:19); 

tmeasure_AIDScases = 11:1/dt:length(tforward);
tmeasure_OPIOIDdeaths = 1:1/dt:length(tforward);
tmeasure_HIVcases = 91:1/dt:length(tforward);
tmeasure_HIVdeaths = 91:1/dt:length(tforward); 
  
  
  %params= [ beta_u, beta_v, mu_u,q_v, mu_v]

    
%results in       

params = [0.374770124132209,0.715206465517946,0.197263367461474,1.00388980628920,0.210020488926058];

lb = [0    0    0    1    0  ]; 
ub = [Inf Inf  Inf  Inf  Inf ];
  
  
   
   q  = 0.03;
   gamma_2 = 0.99921;
   gamma_v = 0.499993;
   alpha = 24.9993;
   mu_2 = 8.8916;
   mu_1 = 46.1262;
   
   mu_a = 0.0879;
   delta_2 = 0.9987;
   q_u= 10.0454;
   delta = 0.0010; 

  for i = 1:3
  [params,fval] =  fminsearchbnd(@err_in_data,params,lb,ub,optimset('Display','iter'));
  end

%initial_cond = [params(13)  params(14) params(15) 10 5 31384];
initial_cond = [7779681 42450 79912 10 5 31384];

 [~, y_r] = ode15s(@(t,y)HIV_Opioid_model(y,params),tforward,initial_cond);

 AidsCases_p = gamma_v*y_r(:,3) + gamma_2*y_r(:,5);

HivCases_p = params(2)*(y_r(:,3)+y_r(:,4)+y_r(:,5)).*y_r(:,1)./(y_r(:,1) +y_r(:,2)+y_r(:,3)+...
     y_r(:,4)+y_r(:,5)+y_r(:,6))+...
     q_u*params(2)*(y_r(:,3)+y_r(:,4)+y_r(:,5)).*y_r(:,2)./(y_r(:,1)+y_r(:,2)+y_r(:,3)+...
     y_r(:,4)+y_r(:,5)+y_r(:,6));

HivCases_end = params(2)*(y_r(tmeasure_HIVcases(:),3)+y_r(tmeasure_HIVcases(:),4)+y_r(tmeasure_HIVcases(:),5)).*...
     y_r(tmeasure_HIVcases(:),1)./(y_r(tmeasure_HIVcases(:),1)+y_r(tmeasure_HIVcases(:),2)+y_r(tmeasure_HIVcases(:),3)+...
     y_r(tmeasure_HIVcases(:),4)+y_r(tmeasure_HIVcases(:),5)+y_r(tmeasure_HIVcases(:),6))+...
     q_u*params(2)*(y_r(tmeasure_HIVcases(:),3)+y_r(tmeasure_HIVcases(:),4)+y_r(tmeasure_HIVcases(:),5)).*...
     y_r(tmeasure_HIVcases(:),2)./(y_r(tmeasure_HIVcases(:),1)+y_r(tmeasure_HIVcases(:),2)+y_r(tmeasure_HIVcases(:),3)+...
     y_r(tmeasure_HIVcases(:),4)+y_r(tmeasure_HIVcases(:),5)+y_r(tmeasure_HIVcases(:),6));
 

  AidsCases_end = gamma_v*y_r(tmeasure_AIDScases(:),3) + gamma_2*y_r(tmeasure_AIDScases(:),5);



 HivDeaths_end = params(5)*y_r(tmeasure_HIVdeaths(:),3) +...
     (1-q)*mu_1*y_r(tmeasure_HIVdeaths(:),4) + (1-q)*mu_2*y_r(tmeasure_HIVdeaths(:),5);
 
 Odeath_end = params(3)*y_r(tmeasure_OPIOIDdeaths(:),2) + ...
     q*mu_1*y_r(tmeasure_OPIOIDdeaths(:),4) + q*mu_2*y_r(tmeasure_OPIOIDdeaths(:),5);
 

  
   w_Aidscases = 1/((mean(AIDSCases)^2)*length(AIDSCases));
  
  w_Hivcases = 1/((mean(HIVCases)^2)*length(HIVCases));
  w_Hivdeaths = 1/((mean(HIVDeaths)^2)*length(HIVDeaths));
  w_Odeaths = 1/((mean(DrugDeaths)^2)*length(DrugDeaths));

   error_in_AidsCases = sum((AidsCases_end - AIDSCases).^2)*w_Aidscases;
   
  
  error_in_HivCases = sum((HivCases_end - HIVCases).^2)*w_Hivcases;
  error_in_HivDeaths = sum((HivDeaths_end - HIVDeaths).^2)*w_Hivdeaths; 
  error_in_ODeaths = sum((Odeath_end - DrugDeaths).^2)*w_Odeaths;
  
  dimc = [0.6 0.6 0.6];

figure(1)
bar(tdata_AIDScases, AIDSCases,'FaceColor',dimc)

%plot(tforward,params(8)*y_r(:,3) + params(11)*y_r(:,5),'LineWidth',2.5)
% hold on 
% H=gca;
% H.LineWidth=2;
%plot(tdata_AIDScases, AIDSCases, 'r.', 'MarkerSize',20)
bar(tdata_AIDScases, AIDSCases,'FaceColor',dimc)
hold on 
plot(tforward,gamma_v*y_r(:,3) + gamma_2*y_r(:,5),'LineWidth',2.5)
ylabel('AIDS Cases','fontweight','normal','fontsize',18);
title({'Number of AIDS diagnoses in US per year', '2000-2018'},'fontweight','normal','fontsize',18)
xlabel('Time(years)','fontweight','normal','fontsize',18)
xline(4,'--k',{'2003'})
xline(7,'--k',{'2006'})
xline(10,'--k',{'2009'})
xline(13,'--k',{'2012'})
xline(16,'--k',{'2015'})
xline(19,'--k',{'2018'})
xticks([ 1 4 7 10 13 16 19])
xticklabels({'2000','2003','2006','2009','2012','2015','2018'})
yticks([ 20000 30000 40000 50000])
yticklabels({'2x10^4', '3x10^4', '4x10^4', '5x10^4'})
set(gca, 'YGrid', 'on', 'XGrid', 'off')
H=gca;
H.LineWidth=2;
hold off




figure(2)
% plot(tforward,HivCases_p,'LineWidth',2.5)
% hold on 
% H=gca;
% H.LineWidth=2;
% plot(tdata_HIVcases, HIVCases, 'r.', 'MarkerSize',20)
bar(tdata_HIVcases, HIVCases, 'FaceColor',dimc)
hold on
plot(tforward,HivCases_p,'LineWidth',2.5)
ylabel('HIV Cases','fontweight','normal','fontsize',18);
title({'Number of HIV diagnoses in US per year', '2008-2018'},'fontweight','normal','fontsize',18)
xlabel('Time(years)','fontweight','normal','fontsize',18)
xline(4,'--k',{'2003'})
xline(7,'--k',{'2006'})
xline(10,'--k',{'2009'})
xline(13,'--k',{'2012'})
xline(16,'--k',{'2015'})
xline(19,'--k',{'2018'})
xticks([ 1 4 7 10 13 16 19])
xticklabels({'2000','2003','2006','2009','2012','2015','2018'})
yticks([ 30000 40000 50000 60000 70000])
yticklabels({'3x10^4', '4x10^4', '5x10^4', '6x10^4','7x10^4'})
set(gca, 'YGrid', 'on', 'XGrid', 'off')
H=gca;
H.LineWidth=2;
hold off


figure(3)
bar(tdata_HIVdeaths, HIVDeaths, 'FaceColor',dimc)
hold on 
plot(tforward,params(5)*y_r(:,3)+(1-q)*mu_1*y_r(:,4) + (1-q)*mu_2*y_r(:,5),'LineWidth',2.5)
ylabel('HIV Deaths','fontweight','normal','fontsize',18);
title({'Number of deaths due to HIV in US per year', '2008-2018'},'fontweight','normal','fontsize',18)
xlabel('Time(years)','fontweight','normal','fontsize',18)
xline(4,'--k',{'2003'})
xline(7,'--k',{'2006'})
xline(10,'--k',{'2009'})
xline(13,'--k',{'2012'})
xline(16,'--k',{'2015'})
xline(19,'--k',{'2018'})
xticks([ 1 4 7 10 13 16 19])
xticklabels({'2000','2003','2006','2009','2012','2015','2018'})
yticks([ 10000 15000 20000 25000 30000])
yticklabels({'1x10^4', '1.5x10^4', '2x10^4', '1.5x10^4','3x10^4'})
set(gca, 'YGrid', 'on', 'XGrid', 'off')
H=gca;
H.LineWidth=2;
hold off

figure(4)
bar(tdata_OPIOIDdeaths, DrugDeaths, 'FaceColor',dimc)
hold on 
plot(tforward,params(3)*y_r(:,2)+ q*mu_1*y_r(:,4) + q*mu_2*y_r(:,5),'LineWidth',2.5)
ylabel('Drug Deaths','fontweight','normal','fontsize',18);
title({'Number of deaths due to drug overdose in US per year', '1999-2018'},'fontweight','normal','fontsize',16)
xlabel('Time(years)','fontweight','normal','fontsize',18)
xline(4,'--k',{'2003'})
xline(7,'--k',{'2006'})
xline(10,'--k',{'2009'})
xline(13,'--k',{'2012'})
xline(16,'--k',{'2015'})
xline(19,'--k',{'2018'})
xticks([ 1 4 7 10 13 16 19])
xticklabels({'2000','2003','2006','2009','2012','2015','2018'})
yticks([ 10000 20000 30000 40000 50000 60000 70000 80000])
yticklabels({'1x10^4', '2x10^4', '3x10^4', '4x10^4','5x10^4','6x10^4','7x10^4','8x10^4'})
set(gca, 'YGrid', 'on', 'XGrid', 'off')
H=gca;
H.LineWidth=2;
hold off
% 
% figure(6)
% plot(tforward, y_r(:,1))
% title('Susceptible')
% figure(7)
% plot(tforward, y_r(:,2))
% title('Opioid users')
% figure(8)
% plot(tforward, y_r(:,3))
% title('HIV')
% 
% figure(9)
% plot(tforward, y_r(:,4))
% title('first infected I_1') 
% 
% figure(10)
% plot(tforward, y_r(:,5))
% title('second infected I_2')
% figure(11)
% plot(tforward, y_r(:,6))
% title('AIDS')
%  

display('Parameters after data fitting: \n ')
 

 fprintf('beta_u = %g\n', params(1));
 fprintf('beta_v =%g\n',  params(2));   
 fprintf('delta = %g\n', delta);
 fprintf('q_u = %g\n', q_u);
 fprintf('mu_u =%g\n',  params(3));
 fprintf('q_v = %g\n', params(4));
 fprintf('mu_v =%g\n',  params(5));
 fprintf('gamma_v =%g\n',  gamma_v);
 fprintf('delta_2 = %g\n',  delta_2);
 fprintf('alpha  = %g\n', alpha);
 fprintf('gamma_2 = %g\n', gamma_2);
 fprintf('mu_a = %g\n', mu_a);
 fprintf('mu_1 = %g\n', mu_1);
 fprintf('mu_2 = %g\n', mu_2);
 fprintf('q(fixed) = %g\n', q);
 fprintf('R_u = %g\n', params(1)/(1/79 + params(3) + delta));
 fprintf('R_v = %g\n', params(2)/(1/79 + params(5) + gamma_v));
%  R_u = params(1)/(1/79 + params(5) + params(3));
%  R_v = params(2)/(1/79 + params(7) + params(8));

function error_in_data = err_in_data(k)

global AIDSCases  HIVCases HIVDeaths DrugDeaths...
       tforward  tmeasure_AIDScases tmeasure_OPIOIDdeaths...
       tmeasure_HIVcases tmeasure_HIVdeaths 
   
   %initial_cond = [k(15)  k(13) k(14) 0 0 31384];
   %initial_cond = [k(13)  k(14) k(15) 10 5 31384];
   %initial_cond = [ 6601176 35729 36653 0 0 31384];
   %initial_cond = [6601176 35729 36653 10 5 31384];
  initial_cond = [7779681 42450 79912 10 5 31384];
   
   q=0.03;

  gamma_2 = 0.99921;
  gamma_v = 0.499993;
  alpha = 24.9993; 
  
  mu_2= 8.8916;
  mu_1 = 46.1262;
  mu_a = 0.0879;
  delta_2 = 0.9987;
  q_u=10.0454;
  delta=0.0010;
 
 [~,y] = ode15s(@(t,y)HIV_Opioid_model(y,k),tforward,initial_cond);
 
 N = (y(tmeasure_HIVcases(:),1) + y(tmeasure_HIVcases(:),2) + y(tmeasure_HIVcases(:),3)+...
      y(tmeasure_HIVcases(:),4) + y(tmeasure_HIVcases(:),5) + y(tmeasure_HIVcases(:),6));
  

  ModelAIDSCases = gamma_v*y(tmeasure_AIDScases(:),3) + gamma_2*y(tmeasure_AIDScases(:),5);

 
 ModelHIVCases = k(2)*((y(tmeasure_HIVcases(:),3) + y(tmeasure_HIVcases(:),4) + y(tmeasure_HIVcases(:),5)).*y(tmeasure_HIVcases(:),1))./N+...
      q_u*k(2)*((y(tmeasure_HIVcases(:),3) + y(tmeasure_HIVcases(:),4) + y(tmeasure_HIVcases(:),5)).*y(tmeasure_HIVcases(:),2))./N;
  
 ModelHIVDeaths = k(5)*y(tmeasure_HIVdeaths(:),3) + (1-q)*mu_1*y(tmeasure_HIVdeaths(:),4) + ...
     (1-q)*mu_2*y(tmeasure_HIVdeaths(:),5);
 
 OPIOIDdeaths = k(3)*y(tmeasure_OPIOIDdeaths(:),2) + q*mu_1*y(tmeasure_OPIOIDdeaths(:),4) +...
     q*mu_2*y(tmeasure_OPIOIDdeaths(:),5) ;

 
 w_aidscases = 1/((mean(AIDSCases)^2)*length(AIDSCases));
 w_HIVcases = 1/((mean(HIVCases)^2)*length(HIVCases));
 w_HIVdeaths = 1/((mean(HIVDeaths)^2)*length(HIVDeaths));
 w_Opdeaths = 1/((mean(DrugDeaths)^2)*length(DrugDeaths));
 
 error_in_data = sum((ModelAIDSCases - AIDSCases).^2)*w_aidscases+...                 
                 sum((ModelHIVCases - HIVCases).^2)*w_HIVcases+...
                 sum((ModelHIVDeaths - HIVDeaths).^2)*w_HIVdeaths+...
                 sum((OPIOIDdeaths - DrugDeaths).^2)*w_Opdeaths; 

%   error_in_data =10*sum((ModelHIVCases - HIVCases).^2)*w_HIVcases+...
%                  1000*sum((ModelHIVDeaths - HIVDeaths).^2)*w_HIVdeaths+...
%                  100*sum((OPIOIDdeaths - DrugDeaths).^2)*w_Opdeaths; 
                
% % 
%  error_in_data = sum((ModelAIDSCases - AIDSCases).^2)*w_aidscases+...                 
%                  sum((ModelHIVDeaths - HIVDeaths).^2)*w_HIVdeaths+...
%                  sum((OPIOIDdeaths - DrugDeaths).^2)*w_Opdeaths; 
  
end

function dy = HIV_Opioid_model(y,params)

   dy=zeros(6,1);
   
  
   mu = 1/79;
    
   Lambda = 7933442*mu;
  
   gamma_2 = 0.99921;
   gamma_v = 0.499993;
   alpha = 24.9993;
   mu_2= 8.8916;
   mu_1 = 46.1262;
   mu_a = 0.0879;
   delta_2 = 0.9987;
   q_u=10.0454;
   delta=0.0010;
  
   beta_u = params(1);
   beta_v = params(2);   
   delta = delta;
   q_u = q_u;
   mu_u = params(3);
   q_v = params(4);
   mu_v = params(5);
  
   delta_2 = delta_2;
   

   mu_a = mu_a;
   mu_1 = mu_1;
   mu_2 = mu_2;

   alpha = alpha;
   gamma_v = gamma_v;
   gamma_2 = gamma_2;
   
   
   S = y(1);
   U = y(2);
   V = y(3);
   I_1 = y(4);
   I_2 = y(5);
   A = y(6);
   N = y(1) + y(2) + y(3) + y(4) + y(5) + y(6);
  

   dy(1) = Lambda - beta_u*(U+I_1+I_2)*S./N - beta_v*(V+I_1+I_2)*S./N- mu*S + delta*U;   
   dy(2) = beta_u*(U+I_1+I_2)*S./N - q_u*beta_v*(V+I_1+I_2)*U./N - (mu + mu_u + delta)*U;
   dy(3) = beta_v*(V+I_1+I_2)*S./N - q_v*beta_u*(U+I_1+I_2)*V./N - (mu + mu_v + gamma_v)*V + delta_2*I_2;
   dy(4) = q_u*beta_v*(V+I_1+I_2)*U./N + q_v*beta_u*(U+I_1+I_2)*V./N - (mu + alpha + mu_1)*I_1;
   dy(5) = alpha*I_1 - (mu + gamma_2+ mu_2 + delta_2)*I_2;
   dy(6) = gamma_v*V + gamma_2*I_2 - (mu + mu_a)*A;

    
end


