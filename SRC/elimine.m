function [tilde_AA, tilde_LL] = elimine(AA, LL, Refneu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% elimination :
% Effectue la pseudo-élimination.
%
% SYNOPSIS [tilde_AA, tilde_LL] = elimination(AA, LL, Refneu)
%          
% INPUT * AA, LL, Refneu : les deux matrices du système ainsi que le tableau
%                          indiquant si un élément est sur le bord d'omega.
%
% OUTPUT - [tilde_AA, tilde_LL]: les matrices obtenues après pseudo-élimination.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:length(Refneu)
   if Refneu(i) == 1
       LL(i) = 0;
       AA(i,:) = 0;
       AA(:,i) = 0;
       AA(i,i) = 1;
    end;
end;
       
tilde_AA = AA;
tilde_LL = LL;
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                     fin de la fonction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2020