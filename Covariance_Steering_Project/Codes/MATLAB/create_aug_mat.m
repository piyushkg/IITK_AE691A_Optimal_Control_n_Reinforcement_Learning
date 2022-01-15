function [A_out,B_out,D_out] = create_aug_mat(A_k,B_k,D_k,A_kp1,B_kp1,D_kp1)
A_out = A_k*A_kp1;
B_out = [A_kp1*B_k, B_kp1];
D_out = [A_kp1*D_k, D_kp1];
end