function [x] = ris_sist_sup(A,b)
    n=size(A,1);
    x = zeros(n,1);
    x(n) = b(n) / A(n, n);
    for i = n-1:-1:1
        x(i) = (b(i)-A(i, i+1:n)*x(i+1:n)) / A(i, i);
    end 
end