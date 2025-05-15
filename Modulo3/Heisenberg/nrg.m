
    
%Numerical Renormalization Group for the 1D XXZ Heisenberg Model
%H = - \sum{ Sx Sx + Sy Sy + delta Sz Sz }

delta = 0.;   % ZZ anisotropy
m = 5;        %Number of initial states
iter = 8;    %Number of iterations


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

%%(NRG)

for i = 1:iter
    
    syslenght = 2^i;
    
    %HSB matrix (superblock)
    Hsb = kron(Hblock, Iblock) + kron(Iblock, Hblock) - delta * kron(Szblock, Szblock) + 0.5 * ( kron(Spblock, Smblock) + kron(Smblock, Spblock) );
    %(Symmetry ensureness)
    Hsb = 0.5 * (Hsb + Hsb');
    
    %Diagonalize Hsb
    LastE = En;
    [eigenvec, eigenval] = eig(Hsb);
    %sort ascend Hsb eigenval and associate eigenvec matrix
    [eigenval, Index] = sort(diag(eigenval), 'ascend');
    eigenvec = eigenvec(:,Index);
    %compute energies
    En = eigenval(1);
    Edens  = En / syslenght;
    Ebond = (En - LastE) / (syslenght/2);
    %truncation and new diagonalizing matrix of Hsb
    nkeep = min(size(eigenval, 1), m);
    Omatrix = eigenvec(:, 1:nkeep);

    %print out results
    fprintf('%d\t%f\t%f\t%f\n', syslenght, En, Ebond, Edens);

    %blocks' update for next iteration
    Szblock = kron(Iblock,Szblock);
    Spblock = kron(Iblock,Spblock);
    Smblock = kron(Iblock,Smblock);
    Iblock = kron(Iblock,Iblock);
    
    %transform the block operators into the truncated basis
    Hblock  = Omatrix' * Hsb * Omatrix;
    Szblock = Omatrix' * Szblock * Omatrix;
    Spblock = Omatrix' * Spblock * Omatrix;
    Smblock = Omatrix' * Smblock * Omatrix;
    Iblock =  Omatrix' * Iblock  * Omatrix;


end
