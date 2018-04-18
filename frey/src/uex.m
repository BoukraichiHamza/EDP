function y = uex (X)
n = size(X,1);
y = zeros(n,1);
for i = 1:n
y(i) = sin ( pi * X(i,1)) * sin ( pi * X(i,2));
end