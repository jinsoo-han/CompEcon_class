ivreg(const mY,const mVariables,const mInstruments,const mWeight)  
{
  decl mInvXZX=invert((mVariables'mInstruments)*mWeight*(mInstruments'mVariables));
  decl vIV;
  if(mInvXZX==0) vIV=constant(.NaN,columns(mVariables),1);
  else vIV=mInvXZX*(mVariables'mInstruments)*mWeight*(mInstruments'mY);  
  return vIV;
}
demand(const aShare,const vDelta,const t,const vParam)
{
  decl i;
  decl id=aProductID[t];
  decl mMu=new matrix[rows(id)][nbGrid];
  for(i=0;i<rows(vParam);i++) mMu+=vParam[i]*(mZ[id][i].*mEta[][i]');
  decl eV=exp(vDelta[id]+mMu);
  decl mS=eV./(1+sumc(eV));
  decl vShat=mS*vWeight;
  decl mWeight=unit(nbGrid).*vWeight;
  aShare[0]=vShat;
  return 1;
}

inverse(const aDelta, const vP)
{
  decl vShat,vDelta=vDelta0;
  decl f,eps=10^(-12);
  decl eps1=0.01;
  decl id,t;
  /* Transform parameter */
  decl vParam=exp(vP);
  decl dParam=exp(vP);
  decl it,maxit=10000;
  decl vIT=new matrix[T][1];
  parallel for(t=0;t<T;t++)
    {
      id=aProductID[t];

      vIT[t]=0;
      do{
	demand(&vShat,vDelta,t,vParam);
	f=log(vShare[id])-log(vShat);
	vDelta[id]=vDelta[id]+f;
	vIT[t]+=1;
      }while(norm(f)>eps && vIT[t]<maxit);
      if(norm(f)>eps) vDelta[id]=constant(.NaN,id);
    }
  aDelta[0]=vDelta;
  return 1;
}
/* GMM objective function */
gmm_obj(const vP, const adFunc, const avScore, const amHessian)
{
  /* Invert demand */
  inverse(&vDelta0,vP);
  /* Quality decomposition */
  decl vLParam=ivreg(vDelta0,mX,mIV,A);  
  decl vXi=vDelta0-mX*vLParam;
  /* GMM objective function */
  decl mG=vXi.*mIV;
  decl g=sumc(mG);
  if(isnan(vDelta0)==1) adFunc[0]=.NaN;
  else adFunc[0]=-g*A*g'/100;
  return 1;
}

