%%% This files contains codes of stochastic model of GABA release.

clear

%% Step1: Model Input & Initial Condition

Gi = [0,100,200]; % #

rand('seed',0) % We set the seed to debug and reproduce our work.

% set the inital four state respectively

state_C0 = zeros(4e4,length(Gi));
state_C1 = zeros(4e4,length(Gi));
state_C2 = zeros(4e4,length(Gi));
state_C3 = zeros(4e4,length(Gi));

tjump(1) = 0; % set the inital time

%% Step2: Run simulations

for i = 1:length(Gi)
    C0 = 100; % #
    C1 = 100; % #
    C2 = 100; % #
    C3 = 0; % #
    G = Gi(i);

    for k=1:40000
    kf1 = 2e3*C0; % #/s
    kb1 = 1e3*C1; % #/s
    kf2 = 5e2*C1;  % #/s
    kb2 = 1e3*C2; % #/s
    kf3 = 1e1*G*C2; % #/s
    kb3 = 5e2*C3; % #/s

    time_kf1 = -log(rand)./kf1; % waiting time for one transition from vesicle to dockingto happen
    time_kb1 = -log(rand)./kb1; % waiting time for one transition from docking to vesicle to happen
    time_kf2 = -log(rand)./kf2; % waiting time for one transition from docking to release to happen
    time_kb2 = -log(rand)./kb2; % waiting time for one transition from release to docking to happen
    time_kf3 = -log(rand)./kf3; % waiting time for one transition from docking to inhibited to happen
    time_kb3 = -log(rand)./kb3; % waiting time for one transition from inhibited to docking to happen

  %%% The minium simulated time interval above is the exact reaction that 
  %%% might happen in a certain period of time.  
  
    time = min([time_kf1,time_kb1,time_kf2,time_kb2,time_kf3,time_kb3]); 

  %%% Update the results according to happened reactions
  
        if time == time_kf1 % transition from vesicle to docking happens
            C0 = C0 - 1;
            C1 = C1 + 1;
        elseif time == time_kb1 % transition from docking to vesicle happens
            C1 = C1 - 1;
            C0 = C0 + 1;
        elseif time == time_kf2 % transition from docking to release happens
            C1 = C1 - 1;
            C2 = C2 + 1;
        elseif time == time_kb2 % transition from release to docking happens
            C1 = C1 + 1;
            C2 = C2 - 1;
        elseif time == time_kf3 % transition from docking to inhibited happens
            C2 = C2 - 1;
            C3 = C3 + 1;
        elseif time == time_kb3 % transition from inhibited to docking happens
            C3 = C3 - 1;
            C2 = C2 + 1;
        end
    
%%% save the resutls
state_C0(k,i) = C0;
state_C1(k,i) = C1;
state_C2(k,i) = C2;
state_C3(k,i) = C3;
tjump(k) = time;
tlast(k) = sum(tjump(1:k)); % update the total reaction time

    end
    
 % iterate the process
 
end

%% Step3: Save/Plot results

for i = 1:length(Gi)
figure
hold on
plot(tlast,state_C0(:,i),'r','LineWidth',1.5)
plot(tlast,state_C1(:,i),'b','LineWidth',1.5)
plot(tlast,state_C2(:,i),'k','LineWidth',1.5)
plot(tlast,state_C3(:,i),'y','LineWidth',1.5)
figurelegend = char('synthesized','docked','released','inhibited');
legend(figurelegend,'Location','EastOutside')
xlabel('Time/s','FontSize',14)
ylabel('Number of GABA vesicles','FontSize',14)
hold off
end

% errorbar
mean = [mean(state_C2(:,1)),mean(state_C3(:,1));mean(state_C2(:,2)),mean(state_C3(:,2));...
    mean(state_C2(:,3)),mean(state_C3(:,3))];
std = [std(state_C2(:,1)),std(state_C3(:,1));std(state_C2(:,2)),std(state_C3(:,2));...
    std(state_C2(:,3)),std(state_C3(:,3))];
gpname = {['Control'],['MOR-siRNA injected'],['Wild Type']};
leg = {['GABA released'],['GABA inhibited']};
figure
barweb(mean,std,1,gpname);
legend(leg,'Location','EastOutside')
xticklabel_rotate([1 2 3],45,gpname)
ylabel('Number of GABA vesicles')
