clear all
close all
clc

global tforward  tmeasure_AIDScases tmeasure_OPIOIDdeaths tmeasure_HIVcases tmeasure_HIVdeaths initial_cond 
   
numiter = 1000; 



          
      %params= [ beta_u, beta_v, mu_u,q_v, mu_v]    
true_params = [0.374770124132209,0.715206465517946,0.197263367461474,1.00388980628920,0.210020488926058];          
           
X = zeros(length(true_params),numiter); 

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

initial_cond = [7779681 42450 79912 10 5 31384];


[~, y_trp] = ode15s(@(t,y)HIV_Opioid_model(y,true_params),tforward,initial_cond);

noiselevel = [0, 0.01, 0.05, 0.1, 0.2, 0.3];

total_ARE =  zeros(length(noiselevel), length(true_params));

total_ARE_Table = {'beta_u', 'beta_v','mu_u','q_v','mu_v'};

for noisei = 1:6
    
rng default
noiselev = noiselevel(noisei)


        for i= 1:numiter
            i
            
             AIDSCases_trp =  gamma_v*y_trp(tmeasure_AIDScases(:),3) + gamma_2*y_trp(tmeasure_AIDScases(:),5);
             
             HivCases_trp =  true_params(2)*(y_trp(tmeasure_HIVcases(:),3) + y_trp(tmeasure_HIVcases(:),4) + y_trp(tmeasure_HIVcases(:),5)).*...
     y_trp(tmeasure_HIVcases(:),1)./(y_trp(tmeasure_HIVcases(:),1)+y_trp(tmeasure_HIVcases(:),2)+y_trp(tmeasure_HIVcases(:),3)+...
     y_trp(tmeasure_HIVcases(:),4)+y_trp(tmeasure_HIVcases(:),5)+y_trp(tmeasure_HIVcases(:),6))+...
     q_u*true_params(2)*(y_trp(tmeasure_HIVcases(:),3)+y_trp(tmeasure_HIVcases(:),4)+y_trp(tmeasure_HIVcases(:),5)).*...
     y_trp(tmeasure_HIVcases(:),2)./(y_trp(tmeasure_HIVcases(:),1)+y_trp(tmeasure_HIVcases(:),2)+y_trp(tmeasure_HIVcases(:),3)+...
     y_trp(tmeasure_HIVcases(:),4)+y_trp(tmeasure_HIVcases(:),5)+y_trp(tmeasure_HIVcases(:),6));
 
            HivDeaths_trp = true_params(5)*y_trp(tmeasure_HIVdeaths(:),3) + (1-q)*mu_1*y_trp(tmeasure_HIVdeaths(:),4) + (1-q)*mu_2*y_trp(tmeasure_HIVdeaths(:),5);

            Drugdeath_trp = true_params(3)*y_trp(tmeasure_OPIOIDdeaths(:),2) + q*mu_1*y_trp(tmeasure_OPIOIDdeaths(:),4) + q*mu_2*y_trp(tmeasure_OPIOIDdeaths(:),5);

  
             AIDSCases =  (noiselev*(AIDSCases_trp).*randn(length(tmeasure_AIDScases),1)) + AIDSCases_trp;
             HIVCases =  (noiselev*(HivCases_trp).*randn(length(tmeasure_HIVcases),1)) + HivCases_trp;
             HIVDeaths =  (noiselev*(HivDeaths_trp).*randn(length(tmeasure_HIVdeaths),1)) + HivDeaths_trp;
             DrugDeaths = (noiselev*(Drugdeath_trp).*randn(length(tmeasure_OPIOIDdeaths),1)) + Drugdeath_trp;
             
             k = true_params;
            lb = [0    0    0    1    0  ];  
             
             k = fminsearchbnd(@(k)err_in_data(k, AIDSCases,HIVCases,HIVDeaths, DrugDeaths),k,lb);
             
             X(:,i) = k';
             
        end
        
        arescore = zeros(1,length(true_params));
    for i = 1:length(true_params)
        arescore(i) = 100*sum(abs(true_params(i) - X(i,:))/abs(true_params(i)))/numiter;
    end
    
    total_ARE(noisei,:) = round(arescore,1);
    total_ARE_Table(noisei+1,:) = num2cell(total_ARE(noisei,:));

end

function error_in_data = err_in_data(k, AIDSCases,HIVCases,HIVDeaths, DrugDeaths)

global tforward  tmeasure_AIDScases tmeasure_OPIOIDdeaths...
       tmeasure_HIVcases tmeasure_HIVdeaths initial_cond 
   
  
   
   
 [~,y] = ode15s(@(t,y)HIV_Opioid_model(y,k),tforward,initial_cond);
 
 N = (y(tmeasure_HIVcases(:),1) + y(tmeasure_HIVcases(:),2) + y(tmeasure_HIVcases(:),3)+...
      y(tmeasure_HIVcases(:),4) + y(tmeasure_HIVcases(:),5) + y(tmeasure_HIVcases(:),6));
  

 ModelAIDSCases = 0.5*y(tmeasure_AIDScases(:),3) + 0.25*y(tmeasure_AIDScases(:),5);

 
 ModelHIVCases = k(2)*((y(tmeasure_HIVcases(:),3) + y(tmeasure_HIVcases(:),4) + y(tmeasure_HIVcases(:),5)).*y(tmeasure_HIVcases(:),1))./N+...
     10.0454*k(2)*((y(tmeasure_HIVcases(:),3) + y(tmeasure_HIVcases(:),4) + y(tmeasure_HIVcases(:),5)).*y(tmeasure_HIVcases(:),2))./N;
  
 ModelHIVDeaths = k(5)*y(tmeasure_HIVdeaths(:),3) + (1-0.03)*46.1262*y(tmeasure_HIVdeaths(:),4) + (1-0.03)*8.8916*y(tmeasure_HIVdeaths(:),5);
 
 OPIOIDdeaths = k(3)*y(tmeasure_OPIOIDdeaths(:),2) + 0.03 * 46.1262*y(tmeasure_OPIOIDdeaths(:),4) + 0.03*8.8916*y(tmeasure_OPIOIDdeaths(:),5) ;

 
 w_aidscases = 1/((mean(AIDSCases)^2)*length(AIDSCases));
 w_HIVcases = 1/((mean(HIVCases)^2)*length(HIVCases));
 w_HIVdeaths = 1/((mean(HIVDeaths)^2)*length(HIVDeaths));
 w_Opdeaths = 1/((mean(DrugDeaths)^2)*length(DrugDeaths));
 
 error_in_data = sum((ModelAIDSCases - AIDSCases).^2)*w_aidscases+...                 
                 sum((ModelHIVCases - HIVCases).^2)*w_HIVcases+...
                 1000*sum((ModelHIVDeaths - HIVDeaths).^2)*w_HIVdeaths+...
                 sum((OPIOIDdeaths - DrugDeaths).^2)*w_Opdeaths;                

  
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
   gamma_v = gamma_v;
   delta_2 = delta_2;
   gamma_2 = gamma_2;
   mu_a = mu_a;
   mu_1 = mu_1;
   mu_2 = mu_2;


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