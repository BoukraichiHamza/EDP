function [M,alpha] = raideur_triangle(Tx,Ty)
%% Fonction qui calcule la matrice de Raideur sur des �l�ments triangles
% Paramétres :
% Tx : abscisses des points du triangles
% Ty : ordonnées des points du triangles
% Retour :
% M : matrice de raideur
% alpha : valeur absolue du determinant de la Jacobienne de la fonction de
% passage.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    M = zeros(3,3);
    matAlpha = [ Tx(2) - Tx(1) , Tx(3) - Tx(1) ; Ty(2) - Ty(1) , Ty(3) - Ty(1)];
    alpha = det(matAlpha);
    for i = 1:3
        gradni = 1/alpha*gradeta(Tx,Ty,i);
        for j = 1:3
            gradnj = 1/alpha*gradeta(Tx,Ty,j);
            M(i,j) = alpha/2 * gradni' * gradnj;
        end
    end
end