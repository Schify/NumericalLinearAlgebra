A = mmread("bcsstk16.mtx");
b = rand(size(A, 1), 1);
x = rand(size(A, 1), 1);
xsol = jacobi_method(A, b, x, 20000);

