function [x, y, ymaj] = generateline(n)
x = 1.5:1:n;
x = [ x; x; repmat(nan,1,(n-1)) ];
y = [ 0.5 n+0.5 nan ].';
y = repmat(y,1,(n-1));

ymaj = [ -5 n+5.5 nan]';
ymaj = repmat(ymaj,1,(n-1));