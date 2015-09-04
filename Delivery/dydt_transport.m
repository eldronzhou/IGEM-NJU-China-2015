%%% This file contains ODEs of pharmacokinetic models. 

function deriv = dydt_transport(t,statevar)
%% Model Parameters

    global AR 
    global partitionbrain partitionliver partitionlung partitionspleen partitionother
    global kblooddis kbloodbind ktransblood km
    global kbindbrain kbindlung kbindliver kbindspleen kbindother
    global kintbrain kintliver kintlung kintspleen kintother
    global kelimtbrain kelimtliver kelimtlung kelimtspleen kelimtother
    global Qbrain  Qlung Qliver Qspleen Qother Qc
    global kescendvec kc
    
%% Model Variables

    Bcf = statevar(1);
    Bcb = statevar(2);
    Etbrain = statevar(3);
    Etlung = statevar(4);
    Etliver = statevar(5);
    Etspleen = statevar(6);
    Etother = statevar(7);
    Tbbrain = statevar(8);
    Tblung = statevar(9);
    Tbliver = statevar(10);
    Tbspleen = statevar(11);
    Tbother = statevar(12);
    brain = statevar(13);
    lung = statevar(14);
    liver = statevar(15);
    spleen = statevar(16);
    other = statevar(17);
    Enc = statevar(18);
  
%% Model Equations
    if Bcf<0
        Bcf = 0;
    end
    
    dBcf = kblooddis*Bcb - kbloodbind*Bcf - ...
        (Etbrain + Etlung + Etliver + ...
        Etspleen + Etother);
    
    dBcb = kbloodbind*Bcf - kblooddis*Bcb;
   
    dEtbrain = ktransblood*partitionbrain*Qbrain/Qc*Bcf - ...
        kbindbrain*Etbrain - km*AR*Etbrain;
     
    dEtlung= ktransblood*partitionlung*Qlung/Qc*Bcf - kbindlung*Etlung ;
     
    dEtliver = ktransblood*partitionliver*Qliver/Qc*Bcf - kbindliver*Etliver;
     
    dEtspleen = ktransblood*partitionspleen*Qspleen/Qc*Bcf - kbindspleen*Etspleen ;
     
    dEtother = ktransblood*partitionother*Qother/Qc*Bcf - kbindother*Etother ;
     
    dTbbrain = kbindbrain*Etbrain  - kintbrain*Tbbrain + km*AR*Etbrain;
    
    dTblung = kbindlung*Etlung - kintlung*Tblung ;
    
    dTbliver = kbindliver*Etliver - kintliver*Tbliver;
    
    dTbspleen = kbindspleen*Etspleen - kintspleen*Tbspleen ;
    
    dTbother = kbindother*Etother - kintother*Tbother ;
     
    dbrain = kintbrain*Tbbrain - kelimtbrain*brain;
    
    dlung = kintlung*Tblung - kelimtlung*lung;
    
    dliver = kintliver*Tbliver - kelimtliver*liver;
    
    dspleen = kintspleen*Tbspleen - kelimtspleen*spleen;
    
    dother = kintother*Tbother - kelimtother*other;
    
    dEnc = kc*kintbrain*Tbbrain - kescendvec*Enc ;
     
    deriv = [dBcf;dBcb;dEtbrain;dEtlung;dEtliver;dEtspleen;...
        dEtother;dTbbrain;dTblung;dTbliver;dTbspleen;dTbother;...
        dbrain;dlung;dliver;dspleen;dother;dEnc];
    
return