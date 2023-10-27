addpath("../..")

 A=mmread("Test.mtx");

 [Q_k,T_k,r,err_ind]=Lanczos_HW1(A, 40, zeros(size(A,1), 1), norm(A, "inf"));