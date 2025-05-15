
    
%Density Matrix Renormalization Group for the 1D XXZ Heisenberg Model
%H = - \sum{ Sx Sx + Sy Sy + delta Sz Sz }

delta = 0.;   % ZZ anisotropy
m = 10;        %Number of initial states
iter = 100;    %Number of iterations


%Intialization (local operators)
I= eye(2);
Sz = [1/2 0 ; 0 -1/2];
Sp = [0 0 ; 1 0];
Sm = [0 1 ; 0 0];

% Initialization (blocks, assuming reflection symmetry)
Szblock = Sz;
Spblock = Sp;
Smblock = Sm;
Iblock  = I;
Hblock  = zeros(2);
En = 0.;

%%(DRG)

for i = 1:iter
    
    syslenght = 2*i + 2;

    %(block + site) operators
    Hblock = kron(Hblock,I) - delta * kron(Szblock, Sz) + 0.5 * ( kron(Spblock, Sm) + kron(Smblock, Sp) );
    Szblock = kron(Iblock, Sz);
    Spblock = kron(Iblock, Sp);
    Smblock = kron(Iblock, Sm);
    Iblock  = kron(Iblock, I);
    %HSB matrix (superblock)
    Hsb = kron(Hblock, Iblock) + kron(Iblock, Hblock) - delta * kron(Szblock, Szblock) + 0.5 * ( kron(Spblock, Smblock) + kron(Smblock, Spblock) );
    %(Symmetry ensureness)
    Hsb = 0.5 * (Hsb + Hsb');
    
    %Diagonalize Hsb
    LastE = En;
    opts.disp = 0;
    opts.issym = 1;
    opts.real = 1;
    [psi, En] = eigs(Hsb, 1, 'SA', opts);
    %compute energies
    Ebond = (En - LastE) / 2;
    Edens = En / syslenght;
    %reduced density matrix
    [rows,cols] = size(psi);
    dim = sqrt(rows);
    psiMatrix = reshape(psi, dim, dim);
    rho = psiMatrix * psiMatrix';

    %diagonalize rho
    [V, D] = eig(rho);
    %sort descend Rho eigenval and associate eigenvec matrix
    [D, Index] = sort(diag(D), 'descend');  
    V = V(:,Index);
    %truncation and new diagonalizing matrix of Hsb
    nkeep = min(size(D, 1), m);
    Omatrix = V(:, 1:nkeep);
    truncErr = 1 - sum(D(1:nkeep));
    
    %print out results
    fprintf('%d\t%f\t%f\t%f\t%f\n', syslenght, En, Ebond, Edens, truncErr);
    
    %transform the block operators into the truncated basis
    Hblock  = Omatrix' * Hblock * Omatrix;
    Szblock = Omatrix' * Szblock * Omatrix;
    Spblock = Omatrix' * Spblock * Omatrix;
    Smblock = Omatrix' * Smblock * Omatrix;
    Iblock =  Omatrix' * Iblock  * Omatrix;

end
    
    