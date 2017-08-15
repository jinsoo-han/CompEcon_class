#BLP function
  #0. export
  #export ivreg, demand, inverse, jacobian, moment, gmm_obj, gmm, gmm_blp

  #1. ivreg
  function ivreg(mY::Array,mX::Array,mIV::Array,mWeight::Array);
    if size((mX'mIV)*mWeight*(mIV'mX))[1]==1;
      vLParam=([1.0]./((mX'mIV)*mWeight*(mIV'mX)))[1]*(mX'mIV)*mWeight*(mIV'mY)[1];
    else;
      vLParam=inv((mX'mIV)*mWeight*(mIV'mX))*(mX'mIV)*mWeight*(mIV'mY);
    end;
    return vLParam;
  end;

  #2. demand (share)
  function demand(vDelta_t::Array,vParam::Array,t);
    mMu=zeros(size(mZ[id[t],1],1),nbGrid);
    for i = 1:size(vParam,1); mMu+=vParam[i]*(mZ[id[t],i].*transpose(mEta[:,i])); end;
    eV=exp.(vDelta_t.+mMu);
    mS=eV./(1+sum(eV,1));
    vShat_t=mS*vWeight;
    #Cross Partial
      #mWeight=eye(nbGrid).*vWeight;
      #mD=diagm((mS.*(1-mS))*vWeight)-(mS*mWeight*mS'-diagm(diag(mS*mWeight*mS')));
      #return vShat_t, mD;
    return vShat_t;
  end;

  #3. inverse
  function inverse(vDelta::Array,vParam::Array,t);
    #set up
    it=0; maxit=10000; eps0=1e-12; f=ones(size(id[t][id[t].==true],1));
    #transform parameters
    vParam=exp.(vParam)
    #initial delta
    vDelta_t=vDelta[id[t]];
    #contraction mapping
    if sum(isnan.(vDelta))==0 && sum(isnan.(vParam))==0;
      while norm(f)>=eps0 && it<=maxit ;
        vShat_t=demand(vDelta_t,vParam,t);
        f=log.(vShare[id[t]])-log.(vShat_t);
        vDelta_t=vDelta_t+f;
        it+=1;
      end;
      if it==maxit ;
        vDelta_t=fill(NaN,size(id[t][id[t].==true],1),1);
      end;
    else;
      vDelta_t=fill(NaN,size(id[t][id[t].==true],1),1);
    end;
    #println("Market=", t), println("Norm=", norm(f)), println("Final Iter=", it);
    return vDelta_t;
  end;

  #4. gmm_obj
  function gmm_obj(vDelta::Array,vParam::Array,A::Array);
    #Inverse / Contraction mapping
      #local
      #temp=map(t->inverse(vDelta,vParam,t),[i for i=1:T]);
      #parallel
      temp=pmap(t->inverse(vDelta,vParam,t),[i for i=1:T]);
      vDelta=zeros(n); for t =1:T; vDelta[id[t],:]=temp[t]; end;
    # Quality Decomp
    if sum(isnan.(vDelta))==0;
      vLParam[:]=ivreg(vDelta,mX,mIV,A);
      if size(vLParam)[1]==1; vXi=vDelta-mX.*vLParam;
      else; vXi=vDelta-mX*vLParam;
      end;
    # GMM Function
      mG=vXi.*mIV;
      g=sum(mG,1);
      f = g*A*g'/100;
    else;
      f = NaN
    end;
    return f = f[1];
  end;

# End of Code
