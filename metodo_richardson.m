function [k, resOut, x] = metodo_richardson(A, b, alpha, Tol, ITMAX)
    % cerchiamo di risolvere Ax=b per x con:
    % A: matrice di dimensioni (n,n) invertibile
    % b: vettore di dimensioni (n,1)
    % alpha: scalare che stabilirà la convergenza del metodo
    %        alpha ottimale = 1/(lambda_1 + lambda_n) 
    %        rispettivamente autovalore massimo e minimo
    % Tol: uno scalare che rappresenta il minimo valore
    %      che la norma del residuo r = b-Ax puo assumere
    % ITMAX: massimo di iterazioni ammissibili
    %        per evitare cicli infiniti
    if nargin == 3
        ITMAX = 10000;
        Tol   = 10e-6;
    end
    n   = size(A,1);
    % M   = (1/alpha)*eye(n);  % matrice diagonale con gli elementi diagonali di A
    % M_1 = alpha*eye(n);
    xk  = rand(n, 1);     % vettore soluzione
    rk  = b - A*xk;

    % vettore: contiene la norma del residuo al tempo k
    res    = zeros(ITMAX,1);
    res(1) = norm(rk);

    xk1 = rand(n,1);
    for k=1:ITMAX
        % xk1 = M_1*(M*xk+(b - M*xk)); 
        xk1 = xk + alpha*rk; 
        rk1 = b - A*xk1;
        res(k) = norm(rk);

        if (res(k) < Tol || norm(xk-xk1) < Tol)
            break
        end
        if res(k) > res(1)*1e+2
            disp(['errore aumenta! ci fermiamo a ',  num2str(k), ' iterazioni (richardson)']);
            break
        end

        xk = xk1;
        rk = rk1;
    end
    resOut = res;
    if nargout == 3
        x{1} = xk;
    end
end