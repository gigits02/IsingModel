
    
%H = - \sum{ delta Sz Sz }

delta = 1.;  
m = 10;        %Number of initial states
iter = 100;    %Number of iterations

I= eye(2);
Sz = [1 0 ; 0 -1];
Sp = [0 0 ; 1 0];
Sm = [0 1 ; 0 0];

H = kron(Sz, Sz);
%(Symmetry ensureness)
H = 0.5 * (H + H');

fprintf('%s\n', 'H:');
disp(H);

%Diagonalize H (LANCZOS)
opts.disp = 0;
opts.issym = 1;
opts.real = 1;
[psi, En] = eigs(H, 1, 'SA', opts);
Edens = En / 2;

%reduced density matrix
[rows,cols] = size(psi);
dim = sqrt(rows);
psiMatrix = reshape(psi, dim, dim);
rho = psiMatrix * psiMatrix';

%show rho
disp(['dimensioni di rho:', num2str(size(rho))]);
fprintf('%s\n', 'rho:');
disp(rho);

%diagonalize rho
[V, D] = eig(rho);

%sort descend Rho eigenval and associate eigenvec matrix
[D, Index] = sort(diag(D), 'descend');  
V = V(:,Index);
disp('autovalori di rho:');
disp(D);
disp('autovettori di rho:');
disp(V);

%compute correlation as <psi(ZZ) psi>
c = psi'*H*psi;
disp('c:');
disp(c);
%compute magnetization as <psi(kron(Z,I) + kron(I,Z))psi>
mz = psi'*(kron(Sz,I) + kron(I, Sz))*psi;
disp('mz:');
disp(mz);


