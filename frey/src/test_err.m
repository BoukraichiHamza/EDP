% Fichier pour calculer et afficher l'erreur entre la solution exacte et la
% solution calculée sur les éléments triangles et aussi de comparer le
% caractère creux de la matrice A et sa factorisation de Cholesky

%% Calcul de l'erreur et le nombre de non z�ros des matrices
close all
n = 20;
err = zeros(n,1);
nnsR = zeros(n,1);
nnsA = zeros(n,1);
nnsP = zeros(n,1);

h = zeros(n,1);
taillem = zeros(n,1);
for i = 4:n
    [u,coord,dirichlet,nnR,nnA,nnP] = cal_sol_maillage_carre(i,@(x) fs(x));
    uexx = uex(coord);
    uexx(dirichlet) = 0;
    err(i) = norm(u - uexx,2)/i;
    nnsR(i) = nnR;
    nnsA(i) = nnA;
    nnsP(i) = nnP;
    taillem(i) = i^2;
    h(i) = 1 / sqrt(length(u));
end

%% Affichage de l'erreur
figure(1)
title("L'erreur en fonction de h.")
hold on
xlabel("log(h)")
hold on
ylabel("log(err(h)")
hold on
plot(log10(h(4:n)),log10(err(4:n)),'k')
hold on
plot(log10(h(4:n)),log10(h(4:n).^2),'g')
hold on
loglog(log10(h(4:n)),log10(h(4:n).^3),'r')
hold on
legend('err','h^2','h^3')

%% Affichage du nombre de non z�ros
figure(2)
title("Le nombre de non z�ro de la matrice en fonction de sa taille.")
hold on
xlabel("Taille de la matrice")
hold on
ylabel("Nombre de non z�ros de la matrice")
hold on
plot(taillem(4:n),nnsR(4:n),'-k')
hold on
plot(taillem(4:n),nnsA(4:n),'-r')
hold on
plot(taillem(4:n),nnsP(4:n),'-o')
legend('R','A','P = Facto de Cholesky de A permut�e')

