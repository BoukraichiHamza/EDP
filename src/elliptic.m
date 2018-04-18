%% Booleans pour le calcul
use_maillage_carre = false;
solve_chol = false;

%% Récupération des données du maillage
if (use_maillage_carre)
    n = 40;
    [coordinates, elements3, elements4, dirichlet, neumann] = maillage_carre(n);
else
    load('coordinates.dat');
    load('dirichlet.dat');
    load('elements3.dat');
    load('elements4.dat');
    load('neumann.dat');
end
m = size(coordinates,1);
l = size(elements3,1);
k = size(elements4,1);
nbarretes = size(neumann,1);

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

%% Calcul des matrices de raideur et assemblage Quadrangle
%Matrice A
for j = 1:k
    T = elements4(j,:);
    C = coordinates(T,:);
    [M,alpha] = raideur_quadrangle(C(:,1),C(:,2));
    A(T,T) = A(T,T) + M;
    bary = sum(C)/4;
    b(T) = b(T) + alpha*f(bary)/4;
end

%% Suite Second membre b
%Conditions de dirichlet
for j=1:m
    b(j) = b(j) - A(j,dirichlet)*u_d(coordinates(dirichlet,:));
end
%Conditions de Newman
for j=1:nbarretes
    Arr = neumann(j,:);
    c = coordinates(Arr,:);
    p1 = c(:,1);
    p2 = c(:,2);
    b(Arr) = b(Arr) + sqrt((p1(1)-p2(1))^2 + (p1(2) - p2(2))^2) * g ( (p1 + p2 )/ 2)/2;
end
%% Resolution du systeme linéaire
Dl = setdiff((1:m),dirichlet);
u = zeros(m,1);
u(dirichlet) = u_d(coordinates(dirichlet));
v = zeros(size(Dl,1),size(Dl,1));
if (solve_chol)
    % Resolution avec factorisation de cholesky
    B = A(Dl,Dl);
    p = symamd(B);
    R = chol(B(p,p));
    v(p) = R\(R'\b(Dl(p)));
else
    % Resolution directe
    v = A(Dl,Dl)\b(Dl);
end
% Affichage
u(Dl) = v;
show(elements3,elements4,coordinates,u);