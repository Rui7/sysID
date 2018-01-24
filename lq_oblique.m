function [proj] = lq_oblique(B,C,A)
% project A onto row space of C along B

kb = size(B,1);
kc = size(C,1);
ka = size(A,1);
[Q,L] = qr([B; C; A]',0);
Q = Q'; L = L';
L11 = L(1:kb,1:kb);
L21 = L(kb+1:kb+kc, 1:kb);
L22 = L(kb+1:kb+kc,kb+1:kb+kc);
L31 = L(kb+kc+1:kb+kc+ka, 1:kb);
L32 = L(kb+kc+1:kb+kc+ka, kb+1:kb+kc);
L33 = L(kb+kc+1:kb+kc+ka, kb+kc+1:kb+kc+ka);
Q1 = Q(1:kb,:);
Q2 = Q(kb+1:kb+kc,:);
Q3 = Q(kb+kc+1:kb+kc+ka,:);
proj = L32*inv(L22)*[L21 L22]*[Q1; Q2];

% check LQ decomposition - all should be zero
%  sum(sum( abs(L11*Q1 - B) >= 1e-06))
%  sum(sum( abs(L21*Q1 + L22*Q2 - C) >= 1e-06))
%  sum(sum( abs(L31*Q1 + L32*Q2 + L33*Q3 - A) >= 1e-06))
end