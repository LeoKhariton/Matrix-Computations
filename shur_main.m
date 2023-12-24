function shur_main()
    clc; clearvars;
    n = 7;
    A = rand(n);
    B = rand(n);
    fprintf('Исходные матрицы:\n');
    catmat(A,B);
    %приводим матрицу А к форме Шура, а матрицу B к треугольному виду
    fprintf('\nМатрицы были приведены к форме Шура:\n');
    [AA, BB, Q, Z]=shur(A, B);
    catmat(AA, BB);
    fprintf('\nТрансформирующие матрицы Q и Z:\n');
    catmat(Q, Z);
    fprintf('\nПроверка:\n');
    disp(Q'*A*Z-AA);
    disp(Q'*B*Z-BB);
    fprintf('\nПри этом матрица Q приводит к форме Шура матрицу AB^-1:\n');
    disp(Q'*A*B^-1*Q);
    fprintf('\nПроверка Q и Z на ортогональность:\n');
    disp(Q*Q'-eye(n));
    disp(Z*Z'-eye(n));
end

function [AA, BB, Q, Z] = shur(A, B)
    [n, ~] = size(A);
    [AA, BB, Q, Z] = hess(A, B);
    if n == 1 || n == 2, return; end
    Q=Q';
    eps = 1e-15;
    
    q = 0;
    while q ~= n 
        for i = 2 : n
            if abs(AA(i, i-1)) <= eps*(abs(AA(i-1, i-1)) + abs(AA(i, i)))
                AA(i, i-1) = 0;
            end
        end
        
        [p, q] = findpq(AA, eps);
        
        if q < n
            A22 = AA(1+p:n-q, 1+p:n-q);
            B22 = BB(1+p:n-q, 1+p:n-q);
            
            if det(B22) == 0
                AA(n-q,n-q-1) = 0;
            else
                [~, ~, Qk, Zk] = qzstep(A22, B22);
                Q = Q * blkdiag(eye(p),Qk,eye(q));
                Z = Z * blkdiag(eye(p),Zk,eye(q));
                AA=blkdiag(eye(p),Qk,eye(q))'*AA*blkdiag(eye(p),Zk,eye(q));
                BB=blkdiag(eye(p),Qk,eye(q))'*BB*blkdiag(eye(p),Zk,eye(q));
            end
        end
    end
end

function [p, q] = findpq(A, eps)
    [n, ~] = size(A);
    
    A33 = zeros(n);
    for i = n-1 : -1 : 3
        if abs(A(i+1, i)) <= eps && abs(A(i,i-1)) > eps && abs(A(i-1,i-2)) > eps
            A33 = A(i+1:end,i+1:end);
            break;
        elseif abs(A(i+1, i)) > eps && abs(A(i,i-1)) > eps && abs(A(i-1,i-2)) > eps
            A33 = [];
            break;
        end
    end
    q = size(A33, 1);
    
    A11=[];
    for i = 1 : n-q-3
        if abs(A(i+1, i)) <= eps && abs(A(i+2,i+1)) > eps && abs(A(i+3,i+2)) > eps
            A11 = A(1:i,1:i);
            break;
        end
    end
    p=size(A11, 1);
end

function [A, B, Q, Z] = qzstep(A, B)
    [n, ~] = size(A);
    C = A*B^-1;
    lambdas = eig(C(n-1:n, n-1:n));
    xyz = (C-lambdas(1)*eye(n))*(C-lambdas(2)*eye(n))*[1; zeros(n-1, 1)];
    
    Q = eye(n);
    Z = eye(n);
    
    for k = 1 : n-2
        Qk = findhouse1(xyz(1:3));
        A = blkdiag(eye(k-1),Qk,eye(n-k-2))*A;
        B = blkdiag(eye(k-1),Qk,eye(n-k-2))*B;
        
        Q = Q * blkdiag(eye(k-1),Qk,eye(n-k-2));
        
        Zk1 = findhouse2(B(k+2,k:k+2)');
        A = A*blkdiag(eye(k-1),Zk1,eye(n-k-2));
        B = B*blkdiag(eye(k-1),Zk1,eye(n-k-2));
        
        Zk2 = findhouse2(B(k+1,k:k+1)');
        A = A*blkdiag(eye(k-1),Zk2,eye(n-k-1));
        B = B*blkdiag(eye(k-1),Zk2,eye(n-k-1));
        
        Z = Z * blkdiag(eye(k-1),Zk1,eye(n-k-2)) * blkdiag(eye(k-1),Zk2,eye(n-k-1));
        
        xyz(1) = A(k+1,k);
        xyz(2) = A(k+2,k);
        if k<n-2
            xyz(3) = A(k+3,k);
        end
    end

    Qn1 = findhouse1(xyz(1:2));
    A=blkdiag(eye(n-2),Qn1)*A;
    B=blkdiag(eye(n-2),Qn1)*B;
    
    Q = Q * blkdiag(eye(n-2),Qn1);

    Zn1 = findhouse2(B(n,n-1:n)');
    A=A*blkdiag(eye(n-2),Zn1);
    B=B*blkdiag(eye(n-2),Zn1);
    
    Z = Z * blkdiag(eye(n-2),Zn1);
end

function Q = findhouse1(x)
    n = length(x);
    x = x/norm(x);
    s = x(2:n)'*x(2:n);
    v = [1; x(2:n)];
    if s==0, beta=0;
    else
        mu = sqrt(x(1)^2 + s);
        if x(1)<=0
            v(1)=x(1)-mu;
        else
            v(1) = -s/(x(1)+mu);
        end
        beta = 2*v(1)^2/(s+v(1)^2);
        v=v/v(1);
    end
    Q = eye(n)-beta*(v*v');
end
function Q = findhouse2(x)
    n = length(x);
    x = x/norm(x);
    s = x(1:n-1)'*x(1:n-1);
    v = [x(1:n-1);1];
    if s==0, beta=0;
    else
        mu = sqrt(x(n)^2 + s);
        if x(n)<=0
            v(n)=x(n)-mu;
        else
            v(n) = -s/(x(n)+mu);
        end
        beta = 2*v(n)^2/(s+v(n)^2);
        v=v/v(n);
    end
    Q = eye(n)-beta*(v*v');
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
