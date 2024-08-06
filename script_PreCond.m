% in questo script confrontiamo la velocità di convergenza e 
% di computazione di: jacobi, richardson e gauss-siedel
clc;
clear;
close all;
N = 100; Tol = 10e-6; ITMAX = 10000;
times1 = zeros(N,1);
times2 = zeros(N,1);
times3 = zeros(N,1);

X = 1:N;
for n=1:N
    %A = rand(n)+1; A = A*transpose(A); % simmetrica definita positiva
    A = gen_Diag_Dom_Matrix1(n); % a_ii < 10 , |a_ij| < 2 i =/= j 
    %A = gen_Diag_Dom_Matrix2(n); % a_ii > 100 , a_ij > 10 i =/= j SPD
    %A = gen_Poisson_Mat1D(n);
    b = ones(n,1);
    
    tic
    [k1,resOut1] = metodo_jacobi(A, b, Tol, ITMAX);
    t1 = toc;
    times1(n) = t1;
    
    tic
    lambdas = eig(A);
    L1 = lambdas(end); Ln = lambdas(1);
    alpha = 2/(L1 + Ln);
    [k2,resOut2] = metodo_richardson(A, b, alpha, Tol, ITMAX);
    t2 = toc;
    times2(n) = t2;

    tic
    [k3,resOut3] = metodo_gauss_siedel_simm(A, b, Tol, ITMAX);
    t3 = toc;
    times3(n) = t3;
end

%interpolazione
grade = 4;
pol1 = polyfit(X, times1, grade);
pol2 = polyfit(X, times2, grade);
pol3 = polyfit(X, times3, grade);
% Valutazione del polinomio interpolante
x_interp = linspace(1, N, 1000); % genero 1000 punti tra 1 e N
y_interp1 = polyval(pol1, x_interp);
y_interp2 = polyval(pol2, x_interp);
y_interp3 = polyval(pol3, x_interp);
% Visualizzazione della forma algebrica del polinomio
polynomial_str1 = poly2str(pol1, 'x'); % Converti i coefficienti in una stringa
disp(['Polinomio interpolante1: ', polynomial_str1])
polynomial_str2 = poly2str(pol2, 'x'); % Converti i coefficienti in una stringa
disp(['Polinomio interpolante2: ', polynomial_str2])
polynomial_str3 = poly2str(pol3, 'x'); % Converti i coefficienti in una stringa
disp(['Polinomio interpolante3: ', polynomial_str3])

predicted_time1 = polyval(pol1, 10000);
predicted_time2 = polyval(pol2, 10000);
predicted_time3 = polyval(pol3, 10000);
disp('predicted times:')
disp([predicted_time1 predicted_time2 predicted_time3])

figure;
semilogy(X,times1,'go', x_interp,y_interp1,'g-', ...
         X,times2,'ro', x_interp,y_interp2,'r-', ...
         X,times3,'bo', x_interp,y_interp3,'b-');
legend('J','J interp', 'R','R interp', 'GS','GS interp');
xlabel('n');
ylabel('log(t)');
title('Confronto tempo impiegato per Ax=b')
disp([ ' J: '  num2str(k1)  ' iterazioni in '  num2str(t1)  ' secondi ' newline ...
       ' R: '  num2str(k2)  ' iterazioni in '  num2str(t2)  ' secondi ' newline ...
       'GS: '  num2str(k3)  ' iterazioni in '  num2str(t3)  ' secondi ' ])


figure;
Xits = (1:ITMAX);
plot(Xits, log(resOut1), 'g-', ...
     Xits, log(resOut2), 'r-', ...
     Xits, log(resOut3), 'b-');

xlabel('Iteration Number'); ylabel('log(Residue)');
title('residuo b-Ax');
legend('J','R','GS');

