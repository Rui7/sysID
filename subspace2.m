% (Naive) Implementation of open-loop N4SID 
% Jodie Simkoff
% 
% References:
% [1] Ljung, Lennart. "System identification." Signal analysis and 
%    prediction. Birkhäuser, Boston, MA, 1998. 163-173.
% [2] Trnka, Pavel. Subspace identification methods. Diss. Ph. D. Thesis, 
%    Czech Technical University in Prague, 2007.
% 

clear all
close all
rng('shuffle')

% define system
Td = 1;
sys = ss(0.9*eye(3,3), 2*eye(3,3), eye(3,3), zeros(3,3),Td);

m = size(sys.B,1);
p = size(sys.C,1);

% N4SID parameters
r = 10; 
s1 = r/2;

% simulate data
N = 1000;
udata = randn(N+r-1,m); 
T = [0:N+r-2]';
[ydata] = lsim(sys, udata, T) + 0.2*randn(N+r-1,p); 

% vectorize data 
yvec = []; uvec = [];
for i = 1:N+r-1
    yvec = [yvec; ydata(i,:)'];
    uvec = [uvec; udata(i,:)'];
end
for i = 1:N
    Yr(1:p*r,i) = yvec((i-1)*p+1:(i+r-1)*(p));
    Ur(1:m*r,i) = uvec((i-1)*m+1:(i+r-1)*(m));
end

% form block Hankel matrices
Yp = Yr(1:p*r/2,:); Yf = Yr(p*r/2+1:end,:);
Ypplus = Yr(1:p*(r/2+1),:); Yfminus = Yr(p*(r/2+1)+1:end,:);

Up = Ur(1:p*r/2,:); Uf = Ur(p*r/2+1:end,:);
Upplus = Ur(1:p*(r/2+1),:); Ufminus = Ur(p*(r/2+1)+1:end,:);

% stacked Hankel matrices
Wp = [Up; Yp]; Wpplus = [Upplus; Ypplus];

% project Yf onto past data (Wp) along row space of Uf
Oh = lq_oblique(Uf, Wp, Yf);
Ohplus = lq_oblique(Ufminus, Wpplus, Yfminus);

W1 = eye(p*s1,p*s1);
W2 = eye(N,N);

% perform singular value decomposition of k-step ahead predictors
[U,S,V]= svd(W1*Oh*W2,'econ');
% determine order of system from singular values
ss = diag(S);
n=find(cumsum(ss)>0.85*sum(ss),1);
U1 = U(:,1:n); S1 = S(1:n,1:n); V1 = V(1:n,:);

% form extended observability matrix
Gammah = W1\U1*sqrtm(S1);
Gammahbar = Gammah(1:end-n,:); 
xhat = pinv(Gammah)*Oh;
xhatplus = pinv(Gammahbar)*Ohplus;

% find state-space model parameters via least-squares
lhs = [xhatplus(:,1:N); ydata((r/2)+1:N+(r/2),:)']';
rhs = [xhat(:,1:N); udata((r/2)+1:N+(r/2),:)']';
params = (rhs'*rhs)\rhs'*lhs;
Ahat = params(1:n,1:n);
Bhat = params(1:n,n+1:n+m);
Chat = params(n+1:n+p,1:n);
Dhat = params(n+1:n+p,n+1:n+m);

% verify that state-space model is (approximately) equivalent to true system
% to within a similarity transform
T = Chat/sys.C;
AhatSim = T*Ahat*inv(T)
BhatSim = T*Bhat
ChatSim = Chat*inv(T)
DhatSim = Dhat

% compute noise covariance and kalman filter gain
sigma = 1/(N-(n+m))*(lhs - rhs*params)'*(lhs-rhs*params);
sigma11 = sigma(1:n,1:n);
sigma12 = sigma(1:n,n+1:n+m);
sigma21 = sigma(n+1:n+p,1:n);
sigma22 = sigma(n+1:n+p,n+1:n+m);

R = sigma22;
K = sigma12*inv(sigma22);
