clear

%%%%%%%%%% Step1: Define Constants
%%%%%%%%%% Model Parameters
 
 %%%% Circulation/extracellular transport
    global AR 
    global partitionbrain partitionliver partitionlung partitionspleen partitionother
    global kblooddis kbloodbind ktransblood km
    global kbindbrain kbindlung kbindliver kbindspleen kbindother
    global kintbrain kintliver kintlung kintspleen kintother
    global kelimtbrain kelimtliver kelimtlung kelimtspleen kelimtother
    global Qbrain  Qlung Qliver Qspleen Qother Qc
    global kescendvec kunpackend kc 
    
    kblooddis = 1e-2/60;
    kbloodbind = 1e-4/60;
    ktransblood = 4e-1;
  
    AR = 0; % #
    km = 2.4e-9; % #/min
    
    partitionbrain = 1e-1;
    partitionliver = 3.5e-1;
    partitionlung = 1.8e-2;
    partitionspleen = 7e-2;
    partitionother = 1e-3;
    
    % min^-1
    kbindbrain = 1;
    kbindlung = 2e-1;
    kbindliver = 1e2;
    kbindspleen = 1;
    kbindother = 1e-1;
    
    % min^-1
    kintbrain = 1e-2;
    kintliver = 3e-2;
    kintlung = 1e2;
    kintspleen = 1e2;
    kintother = 2e1;
    
    % min^-1
    kelimtbrain = 0.6e-1;
    kelimtlung = 0.2e-1;
    kelimtliver = 0.6e-1;
    kelimtspleen = 0.4e-2;
    kelimtother = 0.5e-2;
    
    Qc = 5.58; % L/min
    Qbrain = 0.14*Qc; 
    Qliver = 0.25*Qc;
    Qlung = Qc;
    Qspleen = 0.10*Qc;
    Qother = (1-0.14-0.25-0.1)*Qc;
    
    kescendvec = 1e-1/60;
    kunpackend = 1e-3/60;
    kc = 1e-15*6.02e23; % fmol/ug
    
%%%%%%%%%% Step2 Define Simulation Time
    tlast = 0:0.01:400; %% min

%%%%%%%%%% Step3 Initial Condition
    Bcf =  3.5e4; % ug/L
    Bcf_origin = 3.5e4;
    Bcb = 0;
    Etbrain = 0;
    Etlung = 0;
    Etliver = 0;
    Etspleen = 0;
    Etother = 0;
    Tbbrain = 0;
    Tblung = 0;
    Tbliver = 0;
    Tbspleen = 0;
    Tbother = 0;
    brain = 0;
    lung = 0;
    liver = 0;
    spleen = 0;
    other = 0;
    Enc = 0;
    
    statevar_i = [Bcf,Bcb,Etbrain,Etlung,Etliver,Etspleen,...
        Etother,Tbbrain,Tblung,Tbliver,Tbspleen,Tbother,...
        brain,lung,liver,spleen,other,Enc];
    
%%%%%%%%%% Step4 Run Simulation
    
    [time,statevars] = ode15s(@dydt_transport,tlast,statevar_i);

    Bcf = statevars(:,1);
    Bcb = statevars(:,2);
    Etbrain = statevars(:,3);
    Etlung = statevars(:,4);
    Etliver = statevars(:,5);
    Etspleen = statevars(:,6);
    Etother = statevars(:,7);
    Tbbrain = statevars(:,8);
    Tblung = statevars(:,9);
    Tbliver = statevars(:,10);
    Tbspleen = statevars(:,11);
    Tbother = statevars(:,12);
    brain = statevars(:,13);
    lung = statevars(:,14);
    liver = statevars(:,15);
    spleen = statevars(:,16);
    other = statevars(:,17);
    Enc = statevars(:,18);
   
%%%%%%%%%% Step5 Plot Results

    colors = repmat('krgbmcy',1,300);
    tissue = char('blood','brain','lung','liver','spleen','other'); 
    figure
    hold on
    handle = gcf();
    plot(time,statevars(:,1)/Bcf_origin*100,colors(1),'LineWidth',2.25)
    figurelegend{1} = tissue(1,:);
    for i=2:6
    plot(time, statevars(:,i+11)/Bcf_origin*100,colors(i),'LineWidth',2.25)
    figurelegend{i} = tissue(i,:);
    end
    legend(figurelegend,'Location','Eastoutside')
    xlabel('Time/min','FontSize',14);
    ylabel('Percentage of total exosomes','FontSize',14);
    
    
    %%%% Barplot
    %%% Note that the time is provided by ODE15s solver, the time is
    %%% slightly adjusted.
    
    %%% group time = 1min
    i_1 = find(round(time)==1);
    i_1 = i_1(1);
    group_1 = [brain(i_1) lung(i_1) liver(i_1) spleen(i_1) other(i_1)];
    
    %%% group time = 5min
    i_5 = find(round(time)==5);
    i_5 = i_5(1);
    group_5 = [brain(i_5) lung(i_5) liver(i_5) spleen(i_5) other(i_5)];
    
    %%% group time = 10min
    i_10 = find(round(time)==10);
    i_10 = i_10(1);
    group_10 = [brain(i_10) lung(i_10) liver(i_10) spleen(i_10) other(i_10)];
    
    %%% group time = 30min
    i_30 = find(round(time)==30);
    i_30 = i_30(1);
    group_30 = [brain(i_30) lung(i_30) liver(i_30) spleen(i_30) other(i_30)];
    
    %%% group time = 60min
    i_60 = find(round(time)==60);
    i_60 = i_60(1);
    group_60 = [brain(i_60) lung(i_60) liver(i_60) spleen(i_60) other(i_60)];
    
     %%% group time = 240min
    i_240 = find(round(time)==240);
    i_240 = i_240(1);
    group_240 = [brain(i_240) lung(i_240) liver(i_240) spleen(i_240) other(i_240)];
    
    brain_bar = [group_1(1) group_5(1) group_10(1) group_30(1) group_60(1) group_240(1)];
    lung_bar = [group_1(2) group_5(2) group_10(2) group_30(2) group_60(2) group_240(2)];
    liver_bar = [group_1(3) group_5(3) group_10(3) group_30(3) group_60(3) group_240(3)];
    spleen_bar = [group_1(4) group_5(4) group_10(4) group_30(4) group_60(4) group_240(4)];
    other_bar = [group_1(5) group_5(5) group_10(5) group_30(5) group_60(5) group_240(5)];
    
    y = [brain_bar;lung_bar;liver_bar;spleen_bar;other_bar]/Bcf_origin*100;
  
    figure
    bar(y)
    set(gca,'FontSize',12, 'XTick',[1 2 3 4 5]); 
    xticklabel_rotate([1 2 3 4 5],45,{'Brain','Lung','Liver','Spleen','Other'});
    ylabel('Percentage of total exosomes','FontSize',14);
    legend('1min','5min','10min','30min','60min','240min','Location','Eastoutside')

  