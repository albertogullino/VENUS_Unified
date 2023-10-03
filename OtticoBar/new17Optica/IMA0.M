function [nba]=ima0(nin)
 if length(nin)>0
  ndum=nin;
  fi=find(imag(ndum)>-1e-2);
  ndum(fi)=real(ndum(fi));
  nba=ndum;
 else
  nba=nin;
 end
