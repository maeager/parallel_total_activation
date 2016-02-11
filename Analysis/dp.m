function [p,q,D,sum] = dp(M)  %#codegen
% [p,q] = dp(M) 
%    Use dynamic programming to find a min-cost path through matrix M.
%    Return state sequence in p,q
% 2003-03-15 dpwe@ee.columbia.edu

% Copyright (c) 2003 Dan Ellis <dpwe@ee.columbia.edu>
% released under GPL - see file COPYRIGHT

[r,c] = size(M);

% costs
D = zeros(r+1, c+1);
D(1,:) = NaN;
D(:,1) = NaN;
D(1,1) = 0;
D(2:(r+1), 2:(c+1)) = M;

% traceback
phi = zeros(r,c);
sum=0;
for i = 1:r; 
  for j = 1:c;
    [dmax, tb] = min([D(i, j), D(i, j+1), D(i+1, j)]);
    D(i+1,j+1) = D(i+1,j+1)+dmax;
    phi(i,j) = tb;
  end
end

% Traceback from top left
i = r; 
j = c;
p = zeros(1,r*c); 
q = zeros(1,r*c); 
p(1) = i; p_count=2;
q(1) = j; q_count=2;
while i > 1 & j > 1
  tb = phi(i,j);
  if (tb == 1)
    i = i-1;
    j = j-1;
  elseif (tb == 2)
    i = i-1;
  elseif (tb == 3)
    j = j-1;
  else    
    error;
  end
  p(p_count) = i;p_count=p_count+1;
  q(q_count) = i;q_count=q_count+1;
  sum = sum+M(i,j);
end

% Strip off the edges of the D matrix before returning
D = D(2:(r+1),2:(c+1));
