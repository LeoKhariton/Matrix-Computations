clc; clearvars;

n = 6;
m = 3;
A = rand(n);
B = rand(m);
C = rand(n, m);
Y = zeros(n, m);

%решение встроенной функцией matlab
disp('Решение встроенной функцией: X =');
disp(sylvester(A,B,C));%начиная с версии 2014
disp(sylv(A,B,C));

[Q_hess, A_hess] = hess(A);
[Q_shur, B_shur] = schur(B');

B_shur = B_shur';
CC = Q_hess' * C * Q_shur;

% disp('Матрица А в правой верхней форме Хессенберга:');
% disp(A_hess);
% disp('Матрица B в левой нижней форме Шура:');
% disp(B_shur);

k = m;
while k >= 1
    if k == 1 || abs(B_shur(k-1,k)) <= 1e-15
        s = B_shur(k+1:m,k)'*Y(:,k+1:m)';
        Y(:,k) = inv(A_hess+B_shur(k,k)*eye(n))*(CC(:,k)-s');
        k=k-1;
    else
        s1 = B_shur(k+1:m,k-1)'*Y(1:n,k+1:m)';
        s2 = B_shur(k+1:m,k)'*Y(1:n,k+1:m)';
        YY = inv([A_hess+B_shur(k-1,k-1)*eye(n), B_shur(k,k-1)*eye(n);
                  B_shur(k-1,k)*eye(n),          A_hess+B_shur(k,k)*eye(n)])*[CC(:,k-1)-s1';
                                                                              CC(:,k)-s2'];
        Y(:,k-1) = YY(1:n);
        Y(:,k) = YY(n+1:n*2);
        k = k-2;
    end
end

disp('Решение X =');
X = Q_hess * Y * Q_shur';
disp(X);
disp('Невязка:');
disp(A*X + X*B - C);
