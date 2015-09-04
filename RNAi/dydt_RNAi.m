%%% This files contains ODEs of RNA inference. 

function RNAi=dydt_RNAi(t,at,a,param)
%% Model Parameters
   global kescendvec kdisRISC kformRISC rtot kcleavage kdegRISC kformRISCm
   global kdisRISCm kdeginna kformmRNA kdegmRNA kformprot kdegprot 

%% Model Variables
    R = param(1);
    Cna = param(2);
    C = param(3);
    M = param(4);
    P = param(5);
    
%% Model Equation
    a = interp1(at,a,t); % Reading Enc from external workspace(simulation results in the previous module)
    
    dR = kescendvec*a - kdisRISC*R + kformRISC*(rtot+kdisRISC*R-R-C)*Cna + ...
        kcleavage*C - kdegRISC*(R+C) - kformRISCm*R*M;
    
    dCna = kdisRISC*R - kformRISC*(rtot+kdisRISC*R-R-C)*Cna - kdeginna*Cna;
    
    dC = kformRISCm*R*M - kdisRISCm*C - kdegRISC*(R+C) - kcleavage*C;
    
    dM = kformmRNA + kdisRISCm*C - kdegmRNA*M - kformRISCm*R*M;
    
    dP = kformprot*M - kdegprot*P;
    
    RNAi = [dR;dCna;dC;dM;dP];
    
return