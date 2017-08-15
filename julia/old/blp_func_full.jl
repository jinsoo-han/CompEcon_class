#BLP function module
#module BLP_module
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

  #2. demand / share
  function demand(vDelta_t::Array,vParam::Array,t);
    eV=exp.(vDelta_t.+vParam[1]*aType[id[t],:]);
    mS=eV./(1+sum(eV,1));
    vShat_t=mS*vWeight;
    #Cross Partial
    mWeight=eye(nbGrid).*vWeight;
    mD=diagm((mS.*(1-mS))*vWeight)-(mS*mWeight*mS'-diagm(diag(mS*mWeight*mS')));
    return vShat_t, mD;
  end;

  #3. inverse / contraction mapping
  function inverse(vDelta::Array,vParam::Array,t);
    it=0; maxit=10000; eps0=1e-14; eps1=0.1;
    vDelta_t=vDelta[id[t]];
    vLogShare=log.(vShare[id[t]]);
    f=ones(size(id[t][id[t].==true],1));
    if sum(isnan.(vDelta))==0 && sum(isnan.(vParam))==0;
      while norm(f)>=eps0 && it<=maxit ;
        if norm(f)>=eps1 ;
          vShat_t=demand(vDelta_t,vParam,t)[1];
          f=log.(vShat_t)-vLogShare;
          vDelta_t=vDelta_t-f;
        else;
          vShat_t=demand(vDelta_t,vParam,t)[1];
          mD=demand(vDelta_t,vParam,t)[2];
          f=log.(vShat_t)-vLogShare;
          mD./=vShat_t;
          if sum(isnan.(mD))==0 && rank(mD)==size(id[t][id[t].==true],1);
            mInvD=-inv(mD);
          else;
            mInvD=-eye(size(id[t][id[t].==true],1));
          end;
          vDelta_t=vDelta_t+mInvD*f;
        end;
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

  #4. jacobian
  function jacobian(vDelta,vParam);
    mJacobian=zeros(n);
    if sum(int(isnan(vParam)))==0 && sum(int(isnan(vDelta)))==0;
      for t = 1:T;
        # Demand Func
        eV=exp(vDelta[id[t]].+vParam[1]*aType[id[t],:]);
        mS=eV./(1+sum(eV,1));
        vShat_t=mS*vWeight;
        # Cross-partial Matrix
        mWeight=eye(nbGrid).*vWeight;
        mD=diagm((mS.*(1-mS))*vWeight)-(mS*mWeight*mS'-diagm(diag(mS*mWeight*mS')));
        dShare_dSigma=zeros(size(id[t][id[t].==true],1),1);
        if rank(mD)==size(mD,1);
          mInvD=inv(mD)
          mShare_type=mS.*aType[id[t],:];
          dShare_dSigma[:]=(mShare_type-mS.*sum(mShare_type,1))*vWeight*vParam[1];
          mJacobian[id[t]]=-mInvD*dShare_dSigma;
        else;
          mJacobian=0;
          break;
        end;
      end;
    else;
      mJacobian=0;
    end;
    return mJacobian;
  end;

  #5. moment
  function moment(vDelta::Array,vParam::Array,A::Array);
    # To be replaced to: function moment(adFunc::Array,vDelta::Array,vParam::Array);
    # Inverse Mapping
    #temp=map(t->inverse(vDelta,vParam,t),[i for i=1:T]);
    temp=pmap(t->inverse(vDelta,vParam,t),[i for i=1:T]);
    vDelta=zeros(n); for t =1:T; vDelta[id[t],:]=temp[t]; end;
    if sum(isnan.(vDelta))==0;
    # Quality Decomp
    vLParam[:]=ivreg(vDelta,mX,mIV,A);
      if size(vLParam)[1]==1;
      vXi=vDelta-mX.*vLParam;
      else;
      vXi=vDelta-mX*vLParam;
      end;
    # GMM Function
    mG=vXi.*mIV;
    g=sum(mG,1);
    else;
    g=zeros(1,size(mIV)[2]);
    fill!(g,NaN);
    end;
    return g;
  end;

  #6. gmm_obj
  function gmm_obj(vDelta::Array,vParam::Array,A::Array);
    vG=moment(vDelta,vParam,A);
    if sum(isnan.(vG))==0;
    f = vG*A*vG';
    else
    f = NaN
    end;
    return f = f[1];
  end;

  #6.1. gmm
  function gmm(aDelta::Array,aParam::Array,aLParam::Array,aSE::Array,adFunc::Array,aPvalue::Array,aDF::Array,aJac::Array,aB::Array,aConv);
    # Parameters to be stored
      vDelta0=aDelta;
      vParam=aParam; vLParam=aLParam;
      vSE=aSE; chi2=adFunc; pvalue=aPvalue; df=aDF;
      vB=aB;
      mJacobian=aJac;
      iConv="";

    # Convergence Setup
      step=1; max_step=2;
      it=1; maxit=30;
      eps1=1e-5;

    # Optimal Weighting matrix
      A=inv((vXi.*mIV)'*(vXi.*mIV));

    # GNR Algo maxit times
      while (norm(vB[end])>eps1 && it<=maxit);
        #vDelta0
        temp=map(t->inverse(vDelta0,vParam,t),[1:T]);
        temp1=zeros(n); for t =1:T; temp1[id[t],:]=temp[t]; end;
        if sum(int(isnan(temp1)))!=0 || sum(int(isinf(temp1)))!=0;
          #vParam[:]+=vB[end];
          break;
        end;
        vDelta0[:]=temp1
        #mJacobian
        mJacobian=jacobian(vDelta0,vParam);
        if mJacobian==0;
          break;
        end;
        #IVRegression
        mW=[mX mJacobian];
        vB=ivreg(vDelta0,mW,mIV,A);
        #Newton Step
        vParam[:]-=vB[end];
        vLParam[:]=vB[1:size(mX,2)];
        #Next Iteration
        println("it: ",it, " b:", vB[end])
        it+=1;
      end;

    #= No convergence case
      if norm(vB[end])>eps1 || it>=maxit;
        if it>=maxit;
          iConv="MAX_NOCONV";
        else;
          # Nedler Mead Algo
          println("Nelder Mead Trial");
          @time res=optimize(vParam->gmm_obj(vDelta0,vParam,A), [vParam, vParam], method=:nelder_mead, iterations=1000, xtol=1.0e-5, show_trace=true);
          vParam=[res.minimum[1]];
          # Calculate GNR B
          mJacobian=jacobian(vDelta0,vParam);
          if mJacobian==0;
            iConv="MAX_NOCONV";
          else;
            mW=[mX mJacobian];
            vB=ivreg(vDelta0,mW,mIV,A);
          end;
        end;
      end;
    =#;

    # Report Convergence
      if iConv!="MAX_NOCONV";
        if norm(vB[end])<=eps1 ;
          iConv="MAX_CONV";
        else;
          iConv="MAX_NOCONV";
        end;
      end;

    # Convergence
      if iConv=="MAX_CONV";
        df=size(mIV,2)-size(mX,2)-minimum([1 size(vParam,1)]);
        if df==0; pvalue=1
        else; pvalue=1-cdf(Chisq(df),chi2); end;
      else;
        chi2=[NaN]; pvalue=[NaN]; df=[NaN];
      end;

    # Standard errors
      if iConv=="MAX_CONV" && mJacobian!=0;
        if size(vLParam,1)==1;
        mG=(vDelta0-mX.*vLParam).*mIV;
        m2SLS=(1.0./(mX'mIV)*A*(mIV'mX))*(mX'mIV)*A*mIV';
        else;
        mG=(vDelta0-mX*vLParam).*mIV;
        m2SLS=inv((mX'mIV)*A*(mIV'mX))*(mX'mIV)*A*mIV';
        end;
        mMx=eye(n)-mX*m2SLS;
        dG=((mMx*mJacobian)'*mIV)';
        mVnl=inv(dG'A*dG);
        vSE_nl=sqrt(diag(mVnl))';
      else;
        vSE_nl=NaN;
      end;

    # Linear parameter estimates and variances
      if iConv=="MAX_CONV";
        if size(vLParam,1)==1;
          mVl=1.0./[mX'mIV*A*mIV'mX];
          vSE_l=sqrt(mVl);
        else;
          mVl=inv(mX'mIV*A*mIV'mX);
          vSE_l=sqrt(diag(mVl));
      end;
      else;
        #vLParam=fill(NaN,size(vLParam));
        vSE_l=fill(NaN,size(vLParam));
      end;

    # Store Values
      aParam[:]=vParam;
      aLParam[:]=vLParam;
      aSE[:]=[vSE_l, vSE_nl];
      adFunc[:]=chi2;
      aPvalue[:]=pvalue;
      aDF[:]=df;
      aJac[:]=mJacobian;
      aB[:]=vB;
      aConv[:]=Int(iConv=="MAX_CONV");

    #End of Function
  end;

#6.2. gmm
function gmm_blp(aDelta::Array,aParam::Array,aLParam::Array,aSE::Array,adFunc::Array,aPvalue::Array,aDF::Array,aJac::Array,aB::Array,aConv);
    # Parameters to be stored
      vDelta0=aDelta;
      vParam=aParam; vLParam=aLParam;
      vSE=aSE; chi2=adFunc; pvalue=aPvalue; df=aDF;
      vB=aB;
      mJacobian=aJac;
      iConv="";

    # Convergence Setup
      step=1; max_step=2;
      it=1; maxit=30;
      eps1=1e-5;

    # Optimal Weighting matrix
      A=inv((vXi.*mIV)'*(vXi.*mIV));

    # Simplex
      @time res=optimize(x->gmm_obj(vDelta0,[x],A), 0.0,5.0  ,method=BFGS(), iterations=1000, show_trace=true);
      vParam=[res.minimizer[1]];
      if res.x_converged==true || res.f_converged==true ;
        iConv="MAX_CONV"
      end;

    # Report Convergence
      if iConv!="MAX_NOCONV";
        if norm(vB[end])<=eps1 ;
          iConv="MAX_CONV";
        else;
          iConv="MAX_NOCONV";
        end;
      end;

    # Convergence
      if iConv=="MAX_CONV";
        df=size(mIV,2)-size(mX,2)-minimum([1 size(vParam,1)]);
        if df==0; pvalue=1
        else; pvalue=1-cdf(Chisq(df),chi2); end;
      else;
        chi2=[NaN]; pvalue=[NaN]; df=[NaN];
      end;

    # Standard errors
      if iConv=="MAX_CONV" && mJacobian!=0;
        if size(vLParam,1)==1;
        mG=(vDelta0-mX.*vLParam).*mIV;
        m2SLS=(1.0./(mX'mIV)*A*(mIV'mX))*(mX'mIV)*A*mIV';
        else;
        mG=(vDelta0-mX*vLParam).*mIV;
        m2SLS=inv((mX'mIV)*A*(mIV'mX))*(mX'mIV)*A*mIV';
        end;
        mMx=eye(n)-mX*m2SLS;
        dG=((mMx*mJacobian)'*mIV)';
        mVnl=inv(dG'A*dG);
        vSE_nl=sqrt(diag(mVnl))';
      else;
        vSE_nl=NaN;
      end;

    # Linear parameter estimates and variances
      if iConv=="MAX_CONV";
        if size(vLParam,1)==1;
          mVl=1.0./[mX'mIV*A*mIV'mX];
          vSE_l=sqrt(mVl);
        else;
          mVl=inv(mX'mIV*A*mIV'mX);
          vSE_l=sqrt(diag(mVl));
      end;
      else;
        #vLParam=fill(NaN,size(vLParam));
        vSE_l=fill(NaN,size(vLParam));
      end;

    # Store Values
      aParam[:]=vParam;
      aLParam[:]=vLParam;
      aSE[:]=[vSE_l, vSE_nl];
      adFunc[:]=chi2;
      aPvalue[:]=pvalue;
      aDF[:]=df;
      aJac[:]=mJacobian;
      aB[:]=vB;
      aConv[:]=Int(iConv=="MAX_CONV");

    #End of Function
  end;

# End of module
#end
# End of Code
