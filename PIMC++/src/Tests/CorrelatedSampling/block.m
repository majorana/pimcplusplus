function blockx = block (wA, wAEA, wB, wBEB, factor)
  l = length(wA);
  bs = floor(l/factor);
  blockx = zeros(1,bs);
  bi = 1;
  for first = [1:factor:length(wA)]
    last = min(first+factor-1,length(wA));
    r = [first:last];
    if (bi<=bs)
      blockx(bi) = sum(wAEA(r))/sum(wA(r)) - sum(wBEB(r))/ ...
          sum(wB(r));
    end;
    bi = bi+1;
  end
 
  
  
