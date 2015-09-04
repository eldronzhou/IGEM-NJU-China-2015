%%% This module simulates the knockdown kinetics of RNA inference after
%%% delivery of exosome to central nerve system.This work is based on the
%%% paper written by Bartlett and Davis.
%%% (Nucleic Acids Res 34, 322-333, doi:10.1093/nar/gkj439 (2006))

%%% In order to run thissimulation, you need to run delivery module first 
%%% and set simulation time exactly the same as this file(tsim). Do not 
%%% clear the result when running last module. The input value of Enc is 
%%% based on the output of results in last module.

%% Step1: Define Constants
%%%%%%%%%% Model Parameters
  global kescendvec kdisRISC kformRISC rtot kcleavage kdegRISC kformRISCm
  global kdisRISCm kdeginna kformmRNA kdegmRNA kformprot kdegprot 
        
  kescendvec = 1e-2/60; % min^-1
  kdisRISC = 1e-9/60;   % min^-1
  kformRISC = 2e-19/60; % min^-1
  rtot = 1.9e15;    % #
  kcleavage = 7.2/60;   % min^-1
  kdegRISC = 7.7e-2/60; % min^-1
  kformRISCm = 1.1e14/60;   % #
  kdisRISCm = 1/60; % min^-1
  kdeginna = 2.9e-2/60; % min^-1
  kformmRNA = 5.2*3e13/60;  % min^-1
  kdegmRNA = 2/60;  % min^-1
  kformprot = 5.2e2/60; % min^-1
  kdegprot = 3.5e-1/60; % min^-1
 
  
%%% Step2 Define Simulation Time
   tsim = 0:0.05:24*15*60; %% min
    
%% Step3 Initial Condition
  R = 0;
  Cna = 0;
  C = 0;
  M0 = kformmRNA/kdegmRNA; % steady level #
  P0 = kformprot/kdegprot*kformmRNA/kdegmRNA; % steady level #
  
  param_i = [R,Cna,C,M0,P0];
  
%% Step4 Run Simulation
 
  at = time; % synchronize simulation time 
  a = Enc; % set input of Enc from workplace
  [time,param] = ode15s(@(t,param) dydt_RNAi(t,at,a,param),tsim,param_i);

  
  R = param(:,1);
  Cna = param(:,2);
  C = param(:,3);
  M = param(:,4); 
  P = param(:,5); 
  
%% Step5 Plot Results

  figure
  plot(time/60,P/P0,'LineWidth',2.25)
  xlabel('Time/hours','FontSize',14);
  ylabel('Relative level of MOR Protein','FontSize',14);
  ylim([0,1]);
  figure
  plot(time/60,M/M0,'LineWidth',2.25)
  xlabel('Time/hours','FontSize',14);
  ylabel('Relative level of MOR mRNA','FontSize',14);
  ylim([0,1]);
  
  %%% Then change the input quantity of exosomes in last module, update and
  %%% save the simulation results then make the plot.
 