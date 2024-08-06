function [k, resOut, x] = metodo_jacobi(A, b, Tol, ITMAX)
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
    D   = diag(diag(A));  % matrice diagonale con gli elementi diagonali di A
    D_1 = eye(n)/D;
    E   = A - D;          % creo E come A con diagonale nulla
    xk  = rand(n, 1);     % vettore soluzione
    % xk1 = zeros(n,1);   % serve se uso il ciclo for
    rk  = b - A*xk;

    % vettore: contiene la norma del residuo al tempo k
    res        = zeros(ITMAX,1);
    res(1)     = norm(rk);

    for k=1:ITMAX
        
        xk1 = D_1*(b - E*xk); 
        % for i = 1:n
        %     xk1(i) = (b(i) - A(i, [1:i-1, i+1:n]) * xk([1:i-1, i+1:n])) / A(i, i);
        % end
        rk1 = b - A*xk1;
        res(k) = norm(rk);

        if (res(k) < Tol || norm(xk-xk1) < Tol)
            break
        end
        if res(k) > res(1)*1e+2
            disp(['errore aumenta! ci fermiamo a ',  num2str(k), ' iterazioni (jacobi)']);
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