%%% Run simulation to predict the inhibition of AC. The input of Gi is
%%% based on the experiment result and simuation of activation of MOR. All
%%% the parameters are derived from literature( Science 283, 381-387 (1999)).

clear

%% Step1: Define Constants
%%%%%%%%%% Model Parameters
global kAC KmAC kPDE KmPDE kf1 kb1
kAC = 900; % 1/sec/mM
KmAC = 20; % mM
KmPDE = 19.842; % mM
kPDE = 22.222; % 1/sec/mM
kf1 = 3e-2; % 1/sec/mM
kb1 = 8e-2; % 1/sec/mM

%% Step2 Define Simulation Time
 tsim = 0:0.1:60; % s
  
 %%%%%%%%%% Step3 Initial Condition
 Gi_AC20 = 0;
 cAMP = 3.9; % mM
 cAMP0 = 3.9;
 AC180 = 0.02; % mM
 Gi0 = [1,0.5,0]; % mM
 
 Gi_AC2 = zeros(300,length(Gi0));
 cAMP = zeros(300,length(Gi0));
 AC18 = zeros(300,length(Gi0));
 Gi = zeros(300,length(Gi0));
 
 %%%%%%%%%% Step4 Run Simulation
 
 for i=1:length(Gi0)
 param_i = [Gi_AC20,cAMP0,AC180,Gi0(i)];
 [time,param] = ode15s(@dydt_Signaling,tsim,param_i);
 Gi_AC2(1:length(time),i) = param(:,1);
 cAMP(1:length(time),i) = param(:,2);
 AC18(1:length(time),i) = param(:,3);
 Gi(1:length(time),i) = param(:,4);
 end
 
 %%%%%%%%%% Step5 Plot Results
 
 colors = repmat('krgbmcy',1,300);
 figure
 handel1 = gcf();
 figure
 handel2 = gcf();
 figurelegend = char('Wild Type','MOR-siRNA injected','Control');
 for i=1:length(Gi0)
 figure(handel1)
 hold on
 plot(time,cAMP(1:length(time),i)/cAMP0,colors(i),'LineWidth',2.25)
 hold off
 figure(handel2)
 hold on
 plot(time,AC18(1:length(time),i)/AC180,colors(i),'LineWidth',2.25)
 hold off
 end
 
figure(handel1)
ylim([0.6,1])
xlabel('Time/s','FontSize',14)
ylabel('Relative level of Cellular cAMP','FontSize',14)
legend(figurelegend)
figure(handel2)
ylim([0.6,1])
xlabel('Time/s','FontSize',14)
ylabel('Relative level of AC','FontSize',14)
legend(figurelegend)
 