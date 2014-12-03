function dytheta = dboxcox(y, theta)
  if abs(theta) < 0.000001
    dytheta = 0;
  else
    dytheta = (theta*y.^theta.*log(y) - y.^theta + 1)/theta^2;
  endif