#Demand (share) function
  function demand(vDelta_t::Array,vParam::Array,t);
    mMu=zeros(size(mZ[:,:,t],1),nbGrid);
    for k = 1:size(vParam,1); mMu+=vParam[k]*(mZ[:,k,t]*transpose(mEta[:,k])); end;
    eV=exp.(vDelta_t.+mMu);
    mS=eV./(1+sum(eV,1));
    return vShat_t=mS*vWeight;
  end;

# End of Code
