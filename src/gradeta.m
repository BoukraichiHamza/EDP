%% Fonction qui calcule le gradient des fonctions de bases sur un �l�ment P1
%
% Param�tres :
% Tx : abscisses des points du triangles
% Ty : ordonn�es des points du triangles
% i : indice de la fonction de base
% Retour :
% grad : gradient de la fonction de base
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function grad = gradeta(Tx,Ty,i)
grad =  [Ty(mod(i,3)+1)-Ty(mod(i+1,3)+1);Tx(mod(i+1,3)+1)-Tx(mod(i,3)+1)];
end

