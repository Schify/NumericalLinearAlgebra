rng(42)
n = 1000;
A = rand(n, 100)*rand(100,n);%low rank
b = rand(n, 1);

alpha = 0.7;
k = 50;


[x, res] = problem_2(A, b, alpha, k)


