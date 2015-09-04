%%% This files contains ODEs to simulate inhibition of AC in response to
%%% activation of MOR.

function Signaling = dydt_Signaling(t,param)
 global kAC KmAC kPDE KmPDE kf1 kb1

 Gi_AC2 = param(1);
 cAMP = param(2);
 AC18 = param(3);
 Gi = param(4);

 ATP = 2;
 PDE = 0.45;
 
ReactionFlux1 = kAC*AC18*ATP/(KmAC+ATP);
ReactionFlux2 = kPDE*PDE*cAMP/(KmPDE+cAMP);
ReactionFlux3 = kf1*AC18*Gi - kb1*Gi_AC2;

dGi_AC2 = (ReactionFlux3);
dcAMP = (ReactionFlux1 - ReactionFlux2);
dAC18 = (-ReactionFlux3);
dGi = (-ReactionFlux3);

Signaling=[dGi_AC2;dcAMP;dAC18;dGi];

return