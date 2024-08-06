function [k, resOut, x] = metodo_gauss_siedel_simm(A, b, Tol, ITMAX)
    % cerchiamo di risolvere Ax=b per x con:
    % A: matrice di dimensioni (n,n) invertibile
    % b: vettore di dimensioni (n,1)
    % Tol: uno scalare che rappresenta il minimo valore
    %      che la norma del residuo r = b-Ax puo assumere
    % ITMAX: massimo di iterazioni ammissibili
    %        per evitare cicli infiniti
    if nargin == 2
        ITMAX = 10000;
        Tol   = 10e-6;
    end
    n   = size(A,1);
    E1 = tril(A); % lower triangular
    F1 = A - E1;
    F2 = triu(A); % upper triangular
    E2 = A - F2;
    xk1 = rand(n, 1);
    rk = b - A*xk1;

    % vettore: contiene la norma del residuo al tempo k
    res        = zeros(ITMAX,1);
    res(1)     = norm(rk);

    for k=1:ITMAX
        
        xk2 = ris_sist_inf(E1, b-F1*xk1);
        % xk2 = E1\(b-F1*xk1);
        xk3 = ris_sist_sup(F2, b-E2*xk2);
        % xk3 = F2\(b-E2*xk2);
        rk1 = b - A*xk3;
        res(k) = norm(rk);
        
        if (res(k) < Tol || norm(xk1-xk3) < Tol)
            break
        end
        if res(k) > res(1)*1e+2
            disp(['errore aumenta! ci fermiamo a ',  num2str(k), ' iterazioni (Gauss-siedel)']);
            break
        end

        xk1 = xk3;
        rk = rk1;
    end

    resOut = res;
    if nargout == 3
        x{1} = xk1;
    end
end