psi = [1 2 3 4];
psi = psi';
disp('Array originale:');
disp(psi);

[rows,cols] = size(psi);
dim = sqrt(rows);
psiMatrix = reshape(psi, dim, dim);

disp('Matrice ristrutturata:');
disp(psiMatrix);