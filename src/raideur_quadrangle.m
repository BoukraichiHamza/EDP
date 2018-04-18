function [M,alpha] = raideur_triangle(Tx,Ty)
%% Fonction qui calcule la matrice de Raideur sur des �l�ments quadrangles
%
% Paramétres :
% Tx : abscisses des points du quadrangles
% Ty : ordonnées des points du quadrangles
% Retour :
% M : matrice de raideur
% alpha : valeur absolue du determinant de la Jacobienne de la fonction de
% passage.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Jphi = [ Tx(2) - Tx(1) , Tx(4) - Tx(1) ; Ty(2) - Ty(1), Ty(4) - Ty(1)];
alpha = det(Jphi);
O = inv(Jphi'*Jphi);
a = O(1,1);
b = O(2,1);
c = O(2,2);

M = zeros(4,4);
M(1,1) = 2 * a + 3 *b + 2*c;
M(1,2) = -2*a + c;
M(1,3) = - a  - 3*b - c;
M(1,4) = a - 2*c;
M(2,2) = 2*a-3*b + 2*c;
M(2,3) = a - 2 *c;
M(2,4) = -a + 3*b - c;
M(3,3) = 2 * a + 3 *b + 2*c;
M(3,4) = -2*a + c;
M(4,4) = 2 * a - 3 *b + 2*c;
M = M' + M - diag(M).*eye(4);
M = abs(alpha)/6 * M;
end