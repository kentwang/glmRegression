Yraw = [63, 89;
        54, 91;
        61, 62;
        50, 80;
        52, 72;
        59, 69;
        48, 73;
        74, 81;
        71, 69;
        54, 88;
        48, 92;
        59, 64];
Y = [Yraw([1, 4, 7, 10, 2, 5, 8, 11, 3, 6, 9, 12], 1); Yraw([1, 4, 7, 10, 2, 5, 8, 11, 3, 6, 9, 12], 2)];

Z = [1, -1, 0;
     1, -1, 0;
     1, -1, 0;
     1, -1, 0;
     1, 0, 0;
     1, 0, 0;
     1, 0, 0;
     1, 0, 0;
     1, 1, 0;
     1, 1, 0;
     1, 1, 0;
     1, 1, 0;
     1, -1, 1;
     1, -1, 1;
     1, -1, 1;
     1, -1, 1;
     1, 0, 1;
     1, 0, 1;
     1, 0, 1;
     1, 0, 1;
     1, 1, 1
     1, 1, 1;
     1, 1, 1;
     1, 1, 1;];
s = 6;
t = 4;
b0 = inverse(Z'*Z)*Z'*Y;
bnew = b0;
%bnew = [1, 0, 0]';
bold = ones(length(b0), 1);
epsilon = 10^-6;
iter = 0;


while(max(abs(bnew - bold)) / sum(abs(bold)) > epsilon)
  iter++;
  printf("Iteration %d\n", iter);
  bold = bnew;
  [bnew,sigma2,MVb,EVb,R] = geeNormal(Y,Z,s,t,bnew,"Exchangeable");
  disp(bnew);
endwhile
