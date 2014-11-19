function [b,sigma2,MVb,EVb,R]=fitpoissonGEE(Y,Z,s,t,b0)
%
% Y = concatonated tsx1 response vector
% Z = concatonated tsxp matrix of regressor variables
% s = number of subjects
% t = number of observations per subject
%
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p=size(Z,2);
m=0;
for j=1:s
    y{j}=Y(m+1:m+t,1);
    X{j}=Z(m+1:m+t,:);
    m=m+t;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Estimate mu and A
for j=1:s
    muhat{j}=exp(X{j}*b0);
    A{j}=diag(muhat{j}); % hold on, I think this is the independent case. Not really. The book has a typo
end

% Compute Pearson residuals
for j=1:s
    %
    % Create temporary variables
    mutemp=muhat{j};
    Atemp=A{j};
    ytemp=y{j};
    %
    for i=1:t
        r(j,i)=(ytemp(i)-mutemp(i))/sqrt(Atemp(i,i));
    end
end

% compute scale parameter estimate
C=0;
for j=1:s
    %for i=1:t
     %   C=C+r(j,i)^2;
        C=C+(A{j}^(-1/2)*(y{j}-muhat{j})*(y{j}-muhat{j})'*A{j}^(-1/2));
    %end
end
sigma2=trace(C/(t*s-p));
%sigma2=((C/(t*s-p)));
%sigma2=1;

%%%%%%%%% Compute working correlation matrix
%
% Unspecified
% C=zeros(t,t);
% for j=1:s
%     C=C + A{j}^(-1/2)*(y{j}-muhat{j})*(y{j}-muhat{j})'*A{j}^(-1/2)/(sigma2*(s-p));
% end
%
%
% Exchangable
% alpha=0;
% for j=1:s
%    for i=1:t-1
%        for m=i+1:t
%            alpha = alpha + (r(j,i)*r(j,m))/((0.5*s*(t-1))*p);
%        end
%    end
% end
% alpha=alpha/sigma2;
%
%
% m-dependent
% m=1;
% alpha=0;
% for j=1:s
%     for i=1:t-m
%         alpha=alpha+(r(j,i)*r(j,i+m)/((t-m)*s-p));
%     end
% end
% alpha=alpha/sigma2;
%
%
% AR(1)
% alpha=0;
% for j=1:s
%     for i=1:t-1
%         alpha=alpha+(r(j,i)*r(j,i+1)/((t-1)*s-p));
%     end
% end
% alpha=alpha/sigma2;
%
% Independent
alpha=0;
%
% Fill in working correlation matrix
R=eye(t);
for u=1:t-1
    for v=u+1:t
        %
        % Independent
        R(u,v)=alpha; R(v,u)=alpha;
        %
        % Unspecified
        %R(u,v)=C(u,v); R(v,u)=C(u,v);
        %
        % Exchangable 
        %R(u,v)=alpha; R(v,u)=alpha;
        %
        % m-dependent
        %if abs(u-v)<=m
        %    R(u,v)=alpha; R(v,u)=alpha;
        %else
        %    R(u,v)=0; R(v,u)=0;
        %end
        %
        % AR(1)
        %R(u,v)=alpha^(abs(u-v)); R(v,u)=alpha^(abs(u-v));
    end
end
%
% Display working correlation matrix
%R
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute V
V=[];
for j=1:s
    P{j}=A{j}^(1/2)*R*A{j}^(1/2)*sigma2;
    V=blkdiag(V,P{j});
end

% Compute D
D=[];
for j=1:s
    tempX=X{j};
    tempD=zeros(t,p);
    for i=1:t
        for k=1:p
            tempD(i,k)=tempX(i,k)*exp(tempX(i,:)*b0);
        end
    end
    G{j}=tempD;
    D=[D; tempD];
end

% Adjust coefficient estimates
b=b0+inv(D'*inv(V)*D)*D'*inv(V)*(Y-exp(Z*b0));

% Compute model-based estimate for variance-covariance matrix of b
MVb=inv(D'*inv(V)*D);

% Compute empirical estimator of variance-covariance matrix of b
H=zeros(p);
for j=1:s
    H=H+G{j}'*inv(P{j})*(y{j}-muhat{j})*(y{j}-muhat{j})'*inv(P{j})*G{j};
end
EVb=inv(D'*inv(V)*D)*H*inv(D'*inv(V)*D);

    




