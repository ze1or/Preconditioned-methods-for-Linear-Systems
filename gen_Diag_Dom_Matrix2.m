function A = gen_Diag_Dom_Matrix2(n)
    a = rand(n)+1;
    A = a*transpose(a);
    dnorm = diag(vecnorm(A,1,2)) + eye(n); 
    A = A + dnorm;
end