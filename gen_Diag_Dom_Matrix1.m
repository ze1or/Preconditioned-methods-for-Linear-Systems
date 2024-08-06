function A = gen_Diag_Dom_Matrix1(n)
    A = randn(n,n);
    for i = 1:n
        A(i,i) = sum(abs(A(i,:))) + rand(1);
    end
end
