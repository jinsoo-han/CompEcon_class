  # This code builds consumer types;
  nbGrid=3;
  Grid=zeros(nbGrid); for t =1:nbGrid; Grid[t]=t/(nbGrid+1); end;
  vEta=quantile(Normal(), Grid);
  vWeight=pdf(Normal(),vEta)./sum(pdf(Normal(),vEta));