function ytheta = boxcox(y, theta)
  if abs(theta) < 0.000001
    ytheta = log(y);
  else
    ytheta = (y.^theta - 1)/theta;
  endif