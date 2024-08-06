function A = gen_Poisson_Mat1D(n)
    % Create the main diagonal with 4s
    mainDiag = 4 * ones(n, 1);
    
    % Create the sub and super diagonals with -1s
    offDiag = -1 * ones(n-1, 1);
    
    % Pad the off-diagonal vectors to match the main diagonal length
    offDiagPadded = [offDiag; 0];
    
    % Create the sparse matrix using spdiags
    A = spdiags([offDiagPadded, mainDiag, [0; offDiag]], -1:1, n, n);
end
