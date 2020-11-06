function val = max_val(UU)
  
  val = UU(1);
  for I = 2:length(UU)
    val = max(UU(I), val);
  endfor
  
endfunction
