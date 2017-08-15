# Code: Data Build;
  # This code builds a dataset;

# 0. Setup;
  # Turn this section on if the code needs to be run stand alone;
  #using Distributions;

# 1. Number of markets and products;
  T=1000;
  vJt=zeros(T);
  i=1;
  while i<T+1;
    vJt[i]=50;
    i+=1;
  end;
  n=Int(sum(vJt));

# 2. IDs;
  # PanelID;
  mPanelID=ones(Int(vJt[1]));
  for (index,value) in enumerate(vJt[2:end]);
    mPanelID=[mPanelID; [index+1 for i=1:value]];
  end;
  # ID;
  id=[mPanelID.==t for t in unique(mPanelID)];
  # ProductID
  vProductID=[i for i=1:vJt[1]];
  for t=2:T; vProductID=[vProductID; [i for i=1:vJt[t]]]; end;

# 3. Parameter values;
  #=
    ijt is suppresed
    Model: Util = vDelta + sigma_eta*vEta*vX + eps
                = ('beta0' + 'beta1'*vX + vXi) + sigma_eta*vEta*vX + eps
                where vX ~ N(0,'sigma_x')
                      eta ~ N(0,1), 'Sim' simulations
                      vXi ~ N(0,'sigma_xi')
                      eps ~ N(0,'sigma_eps')
  =#;
  sigma_eps=1.0;
    #Parameters to be estimated;
  beta0=-5.0/sigma_eps; beta1=1.0/sigma_eps; sigma_eta=2.0/sigma_eps;
    #Parameters
  sigma_x=1.0; sigma_xi=1;

# 4. Simulate data;
  # Characteristics: X;
    vX=rand(Normal(0,sigma_x),n);
  # vXi;
    vXi=rand(Normal(0,sigma_xi),n);
  # vDelta;
    vDelta=beta0+beta1*vX+vXi;
  # aType;
    aType=vX*transpose(vEta);

# 5. Market shares
  # Share
  temp=map(t->demand(vDelta[id[t]],[sigma_eta],t)[1],[i for i=1:T]);
  vShare=zeros(n); for t =1:T; vShare[id[t],:]=temp[t]; end; temp=0;
  vShare_outside=zeros(n); for t=1:T; vShare_outside[id[t]]=1-sum(vShare[id[t]]); end;
  vDelta0=log.(vShare)-log.(vShare_outside);

# 6. IV Construction
  # BLPIV;
  mBLP_IV=zeros(n);
  for t=1:T; nb=Int(vJt[t]);
  mBLP_IV[id[t]]=sum(vX[id[t]])-vX[id[t]];
  end;

  # DiffIV_dist = sqrt(sum((X_j-x_j').^2,2))
  mDiff_IV_dist=zeros(n);
  for t=1:T; mDiff_IV_dist[id[t]]=sum((vX[id[t]].-transpose(vX[id[t]])).^2,2); end;

  #= DiffIV_
  sd_x=sqrt(var(vX));
  mDiff_IV_dummy_halfsd=zeros(n);
  mDiff_IV_dummy_onesd=zeros(n);
  for t=1:T; nb=vJt[t];
  mDiff_IV_dummy_halfsd[id[t]]=float(sum(abs.(vX[id[t]].-transpose(vX[id[t]])).<0.5*sd_x)-eye(nb),2);
  mDiff_IV_dummy_onesd[id[t]]=float(sum(abs(vX[id[t]].-transpose(vX[id[t]])).<sd_x)-eye(nb),2);
  end;
  =#;



# End of Code
