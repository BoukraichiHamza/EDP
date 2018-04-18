%% Fonction qui calcul une solution approch�es sur les �l�ments triangles
% Paramétres :
% n : taille du maillage
% f : fonction du second membre
% Retour :
% u : solution approchée
% coordinates : position des éléments des triangles
% dirichlet : indices des éléments appartenant à la frontière de dirichlet
% nnR : nombres de non zeros de la matrice R = chol(A)
% nnA : nombres de non zeros de la matrice A
% nnP : nombres de non zeros de la matrice P = chol(A permutée)
function [u,coordinates,dirichlet,nnR,nnA,nnP] =  cal_sol_maillage_carre(n,f)
%% Récupération des données du maillage
[coordinates, elements3, elements4, dirichlet, neumann] = maillage_carre(n);
m = size(coordinates,1);
l = size(elements3,1);

%% Calcul des matrices de raideur et assemblage Triangle
%Matrice A
A = sparse(m,m);
b = zeros(m,1);
for j = 1:l
    T = elements3(j,:);
    C = coordinates(T,:);
    [M,alpha] = raideur_triangle(C(:,1),C(:,2));
    A(T,T) = A(T,T) + M;
    bary = sum(C)/3;
    b(T) = b(T) + alpha*f(bary)/6;
end


%% Suite Second membre b
%Conditions de dirichlet
for j=1:m
    b(j) = b(j) - A(j,dirichlet)*u_d(coordinates(dirichlet,:));
end
%% Resolution du systeme linéaire
Dl = setdiff((1:m),dirichlet);
u = zeros(m,1);
u(dirichlet) = u_d(coordinates(dirichlet));
% Resolution
nnR = nnz(chol(A(Dl,Dl)));
nnA = nnz(A(Dl,Dl));
B = A(Dl,Dl);
p = symamd(B);
nnP = nnz(chol(B(p,p)));
v = A(Dl,Dl)\b(Dl);
u(Dl) = v;

    