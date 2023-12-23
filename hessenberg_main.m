function hessenberg_main()
    clc; clearvars;
    n = 5;
    A = rand(n);
    B = rand(n);
    fprintf('Исходная пара матриц:\n');
    catmat(A, B);
    %приводим матрицу А к форме Хессенберга, а матрицу B к треугольному виду
    fprintf('\nМатрицы были приведены к форме Хессенберга:\n');
    [AA, BB, Q, Z] = hessenberg(A, B);
    catmat(AA, BB);
    fprintf('\nТрансформирующие матрицы Q и Z:\n');
    catmat(Q, Z);
    fprintf('\nПроверка невязок для А и В:\n');
    disp(Q'*A*Z-AA);%проверка невязки для А
    disp(Q'*B*Z-BB);%проверка невязки для В
    fprintf('\nПроверка матриц Q и Z на ортогональность:\n');
    disp(Q*Q'- eye(n));
    disp(Z*Z'- eye(n));
    %решение с помощью встроенной функции matlab
    [AA, BB, Q, Z] = hess(A, B);
    fprintf('\nМатрицы были приведены к форме Хессенберга с помощью встроенной функции:\n');
    catmat(AA, BB);
    fprintf('\nМатрицы Q и Z:\n');
    catmat(Q, Z);
end

function [AA, BB, Q, Z] = hessenberg(A, B)
    [n, ~] = size(A);
    [Q, ~] = qrgivens(B);
    AA = Q' * A;
    BB = Q' * B;
    Z = eye(n);
    for k = 1 : n-2
        for i = n : -1 : k+2
            Rl = eye(n);
            [c, s] = givensrotation(AA(i-1, k), AA(i, k));
            Rl([i-1, i], [i-1, i]) = [c -s; s c];
            Q = Q * Rl;
            AA = Rl'*AA;
            BB = Rl'*BB;
%             catmat(AA,BB);fprintf('\n');

            Rr = eye(n);
            [c, s] = givensrotation(-BB(i, i), BB(i, i-1));
            Rr([i-1, i],[i-1, i]) = [c -s; s c];
            Z = Z * Rr;
            AA = AA*Rr;
            BB = BB*Rr;
%             catmat(AA,BB);fprintf('\n');
        end
    end
end

function catmat(A, B)
    [n,m]=size(A);
    for i = 1:n
        for j = 1:m
            fprintf('%.4f\t\t', A(i, j));
        end
        fprintf('|\t\t');
        for j = 1:n
            fprintf('%.4f\t\t', B(i, j));
        end
        fprintf('\n');
    end
end

function [Q,R] = qrgivens(A)
    [m,n] = size(A);
    Q = eye(m);
    R = A;
    for j = 1:n
        for i = m:-1:(j+1)
            G = eye(m);
            [c,s] = givensrotation(R(i-1,j),R(i,j));
            G([i-1, i],[i-1, i]) = [c -s; s c];
            R = G'*R;
            Q = Q*G;
        end
    end
end

function [c,s] = givensrotation(a,b)
    if b == 0
        c = 1;
        s = 0;
    else
        if abs(b) > abs(a)
            r = a / b;
            s = 1 / sqrt(1 + r^2);
            c = s*r;
        else
            r = b / a;
            c = 1 / sqrt(1 + r^2);
            s = c*r;
        end
    end
end
