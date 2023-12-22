function hessenberg-main()
    clc; clearvars;
    n = 5;
    A1 = rand(n);
    A2 = rand(n);
    fprintf('Исходная пара матриц:\n');
    catmat(A1, A2);
    %приводим матрицу А1 к форме Хессенберга, а матрицу А2 к треугольному виду
    fprintf('\nМатрицы были приведены к форме Хессенберга:\n');
    [AA, BB, Q, Z]=hessenberg(A1, A2);
    catmat(AA, BB);
    fprintf('\nТрансформирующие матрицы Q и Z:\n');
    catmat(Q, Z);
    fprintf('\nПроверка:\n');
    catmat(Q'*A1*Z-AA, Q'*A2*Z-BB);
    Q'*A1*Z-AA
    Q'*A2*Z-BB
%     fprintf('\nПроверка матриц Q и Z на ортогональность:\n');
%     catmat(Q'*Q, Z'*Z);
    fprintf('\nПроверка невязки Q и Z\n');
    catmat(Q*Q'- eye(n), Z*Z'- eye(n));
    Q*Q'- eye(n)
    Z*Z'- eye(n)
    %с помощью встроенной функции matlab
    [AA, BB, Q, Z] = hess(A1, A2);
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
