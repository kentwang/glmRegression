%- Vectorization of nchoosek()
function r = NchooseK(N, K)
  r = zeros(length(N), 1);
  for i = 1:length(r)
    r(i) = nchoosek(N(i), K(i));
  endfor