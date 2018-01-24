function [proj, projcompl] = lq_ortho(U,Y)
% return projection of Y onto U and onto orthogonal complement of U

km = size(U,1);
kp = size(Y,1);
[Q,L] = qr([U;Y]',0);
Q = Q'; L = L';
L11 = L(1:km,1:km);
L21 = L(km+1:km+kp, 1:km);
L22 = L(km+1:km+kp,km+1:km+kp);
Q1 = Q(1:km,:);
Q2 = Q(km+1:end,:);
proj = L21*Q1;

projcompl = L22*Q2;
end