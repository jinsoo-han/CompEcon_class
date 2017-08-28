#include<oxstd.h>
#import <maximize>
#import <solvenle>
#include <oxdraw.h>
#include <oxfloat.h>
#include <oxprob.h>

/* Declare common variables */
decl n,T,aProductID;
decl vShare,vShare0,mZ,mX,mIV;
decl vDelta0,A;
decl mEta,nbGrid,vWeight;

#import<blp_func>

main()
{
  /* Sample size */
  T=20;
  decl J=50;
  n=J*T;
  decl nbRC=2;
  println("////////////////////");
  println("Specification");
  println("nb RC = ",nbRC);
  println("////////////////////");
  
  /* Number of variables and random-coefficients */
  decl nbX1=1;
  decl K=nbRC+nbX1;
  decl step=0.25;
  decl sig_xi=1;
  decl sig_x=1;

  /* Linear preference parameters */
  decl iBeta0=-5;
  decl vLparam_true=(iBeta0|constant(1,K,1));
  println("Linear parameters: ",vLparam_true);

  /* Non-linear preference parameters */
  decl sig=2;
  decl vNLparam_true=log(constant(sig^2,nbRC,1));
  println("true nl param: ",vNLparam_true);
  
  /* Integration grid */
  nbGrid=100;
  mEta=rann(nbGrid,nbRC);
  vWeight=constant(1/nbGrid,nbGrid,1);

  /* MC Param */
  decl vSigma_x=constant(sig_x,K,1);
  decl mSigma=new matrix[K][K];
  mSigma=setdiagonal(mSigma,vSigma_x.^2);
  decl mCholeski_Sigma=choleski(mSigma);
  mX=ones(n,1)~(rann(n,K)*mCholeski_Sigma');
  mZ=mX[][1+K-nbRC:];
  decl vXi_true=rann(n,1)*sig_xi;
  decl vDelta_true=mX*vLparam_true+vXi_true;
  vDelta0=vDelta_true;

  /* Panel Structure */
  decl i;
  decl mPanelID=new matrix[n][2];
  aProductID=new array[T];
  for(i=1;i<=T;i++)
    {
      aProductID[i-1]=range((i-1)*J,i*J-1)';
      mPanelID[range((i-1)*J,i*J-1)'][]=range(1,J)'~constant(i,J,1);
    }
  
  decl mJacobian=1;
  decl vShare_t;
  vShare=new matrix[n][1];
  decl vShare0=new matrix[n][1];
  for(i=0;i<T;i++)
    {
      demand(&vShare_t,vDelta_true,i,exp(vNLparam_true));
      vShare[aProductID[i]]=vShare_t;
      vShare0[aProductID[i]]=constant(1-sumc(vShare_t),vShare_t);
    }
  println("distribution of share 0: ",quantilec(vShare0,<0,1/4,1/2,3/4,1>)');
  decl mInvJacobian=-invert(mJacobian);
  

  decl time0=timer();
  vDelta0=log(vShare);
  inverse(&vDelta0,vNLparam_true);
  println("norm between true and inverted delta: ",norm(vDelta0-vDelta_true));
  println("cpu time (inverse): ",(timer()-time0)/100);

  /* Instruments */
  decl id,k,mD;
  decl mDiffIV=new matrix[n][columns(mZ)];
  for(i=0;i<T;i++)
    {
      id=aProductID[i];
      for(k=0;k<columns(mZ);k++)
	{
	  mD=setdiagonal(mZ[id][k]-mZ[id][k]',-mZ[id][k]');
	  mDiffIV[id][k]=sumr(mD.^2);
	}
    }
  println("mean/sd shares: ",meanc(vShare)~sqrt(varc(vShare)));
  println("mean/sd delta: ",meanc(vDelta_true)~sqrt(varc(vDelta_true)));
  /* Save data */
  savemat("blp_example_param.dta",vNLparam_true);
  savemat("blp_example_data.dta",mPanelID~vShare~vShare0~vDelta_true~vXi_true~mX~mDiffIV,
	  {"productid","marketid","share","share0","delta_true","xi_true","intercept","x1","x2","x3","z1","z2","z3"});
  savemat("blp_example_grid.dta",mEta~vWeight,{"eta1","eta2","weight"});
  savemat("blp_example_data.csv",mPanelID~vShare~vShare0~vDelta_true~vXi_true~mX~mDiffIV,
	  {"productid","marketid","share","share0","delta_true","xi_true","intercept","x1","x2","x3","z1","z2","z3"});
  savemat("blp_example_grid.csv",mEta~vWeight,{"eta1","eta2","weight"});
  
}