function r=sf3(p)
%  sf3.m   2-19-96
%
%  Spatial Error Analysis ToolBox Version 1.0,   October 5, 1997
%  Copyright 1997-1998 by David Y. Hsu  All rights reserved.
%  dhsu@littongcs.com
%
% given p and f(x), find r such that      F(r)-p = 0,
% where F(r) = integral of f(x) from 0 to r.
%      Maxwell pdf used for f(x) here.
%
if p>1 | p<0
  disp( '0 <= p <= 1')
  error( 'input error for sf3(p)')
end

  if p==1
    r=Inf;
  else
    r=1;
    p0=sf2(r);
%%%%%%%%%  searching for a, b  %%%%%%%%%%%%%%%%
for j=1:2:5
  while (p0 > p)
    r=r-(0.1)^j ;
    p0=sf2(r);
  end

  while (p0 < p)
    r=r+(0.1)^j ;
    p0=sf2(r);
  end
end
end
