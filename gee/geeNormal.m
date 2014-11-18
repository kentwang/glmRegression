function [b,sigma2,MVb,EVb,R]=geeNormal(Y,Z,s,t,b0)
%
% Y = concatonated tsx1 response vector
% Z = concatonated tsxp matrix of regressor variables
% s = number of subjects
% t = number of observations per subject
% b0 = estimated coefficients from last iteration
%

  y = {};
  X = {};
  muhat = {};
  r = zeros(s, t);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  p=size(Z,2);
  m=0;
  for j=1:s
      y{j}=Y(m+1:m+t,1);
      X{j}=Z(m+1:m+t,:);
      m=m+t;
  endfor
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



  % Estimate mu and A
  for j=1:s
      muhat{j}=X{j} * b0;
  endfor


  for j = 1:s
    r(j, :) = (y{j} - muhat{j})';
  endfor

  % compute scale parameter estimate
  C=0;
  for j=1:s
    C=C+((y{j}-muhat{j})*(y{j}-muhat{j})');
  end
  sigma2=trace(C/(t*s-p));


  %%%%%%%%% Compute working correlation matrix
  % Exchangable
%  alpha=0;
%  for j=1:s
%    for i=1:t-1
%      for m=i+1:t
%        alpha = alpha + (r(j,i)*r(j,m))/((0.5*s*(t-1))*p);
%      end
%    end
%  end
%  alpha=alpha/sigma2;

%   AR(1)
   alpha=0;
   for j=1:s
       for i=1:t-1
           alpha=alpha+(r(j,i)*r(j,i+1)/((t-1)*s-p));
       end
   end
   alpha=alpha/sigma2;

  % Fill in working correlation matrix
  R=eye(t);
  for u=1:t-1
      for v=u+1:t
          % Exchangable 
          % R(u,v)=alpha; R(v,u)=alpha;
          % AR(1)
          R(u,v)=alpha^(abs(u-v)); R(v,u)=alpha^(abs(u-v));
      end
  end

  % Compute RStart
  RStar = [];
  for j=1:s
    RStar = blkdiag(RStar, R);
  endfor

  % Adjust coefficient estimates
  %b=b0+inv(D'*inv(V)*D)*D'*inv(V)*(Y-exp(Z*b0)); This is wrong for normal case
  b = inv(Z'*inv(RStar)*Z)*Z'*inv(RStar)*Y;

  % Compute model-based estimate for variance-covariance matrix of b
  MVb=inv(Z'*inv(RStar)*Z);

  % Compute empirical estimator of variance-covariance matrix of b
  % H=zeros(p);
  H = [];
  for j=1:s
    H = blkdiag(H, (y{j}-muhat{j})*(y{j}-muhat{j})');
  end
  EVb=inv(Z'*inv(RStar)*Z)*Z'*inv(RStar)*H*inv(RStar)*Z*inv(Z'*inv(RStar)*Z);   

endfunction


