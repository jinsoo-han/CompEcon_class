# Main Code: Demand Calculation;

# 0. Introduction;
	# Number of Processors;
		#addprocs();
	# Load standard packages;
		@everywhere using DataFrames;
	# Load user-written packages;
		@everywhere include("$(homedir())/Dropbox/CompEcon_class/julia/exercise/fn_demand.jl");

# 1. Prep data;
	# Load characteristics data;
	@everywhere global mData=readtable("$(homedir())/Dropbox/CompEcon_class/julia/blp/blp_data/blp_example_data.csv");
	# Check data with eye;
	#head(mData); #tail(mData);
	# Define market parameters;
	@everywhere global T=size(unique(mData[:marketid]),1);
	@everywhere global J=size(unique(mData[:productid]),1);
	@everywhere global n=J*T;
	# Define covariate matrix with the dimesion J*size(mX,2)*T;
	@everywhere global mX=zeros(J,4,T);
	@everywhere for t=1:T; mX[:,:,t]=Array(mData[mData[:marketid].==t,[:intercept, :x1, :x2, :x3]]); end;
	@everywhere global mZ=zeros(J,2,T);
	@everywhere for t=1:T; mZ[:,:,t]=Array(mData[mData[:marketid].==t,[:x2, :x3]]); end;
	@everywhere global xi=zeros(J,1,T);
	@everywhere for t=1:T; xi[:,:,t]=Array(mData[mData[:marketid].==t,[:xi_true]]); end;
	# Load integration grid/weight;
	@everywhere global mGrid=readtable("$(homedir())/Dropbox/CompEcon_class/julia/blp/blp_data/blp_example_grid.csv");
	@everywhere global mEta=Array(mGrid[:,[:eta1,:eta2]]);
	@everywhere global vWeight=Array(mGrid[:weight]);
	@everywhere global nbGrid=size(vWeight,1);
	# Store true share and out-side share;
	@everywhere global vShare=Array(mData[:,:share]);
	@everywhere global vShare0=Array(mData[:,:share0]);
	# Define marketid;
	marketid=Array(mData[:,:marketid]);

# 2. Evaluate the demand f'n;
	# Define parameters;
	@everywhere vNLparam=[1.4, 1.4];
	@everywhere vLparam=[-5.0, 1.0, 1.0, 1.0];
	@everywhere vNLparam2=[1.3863, 1.3863];
	println("vNLparam : ",vNLparam, "vLparam : ",vLparam);

	# Pre-calculate vDelta;
	@everywhere vDelta=zeros(J,1,T);
	@everywhere for t=1:T; vDelta[:,:,t]=mX[:,:,t]*vLparam+xi[:,:,t]; end;
	#temp=zeros(n); for t =1:T; temp[marketid.==t,:]=vDelta[:,:,t]; end;
	#norm(Array(mData[:,:delta_true])-temp);

	# Calculate demand;
	temp=map(t->demand(vDelta[:,:,t],vNLparam,t),[i for i=1:T]);
	vShare_est=zeros(n); for t =1:T; vShare_est[marketid.==t,:]=temp[t]; end;


# 3. Parallelization
	# Add processors;
	addprocs();

	# measure the time difference b/t single-processor vs. n-processor operation
	@time temp=map(t->demand(vDelta[:,:,t],vNLparam,t),[i for i=1:T]);
	vShare_est=zeros(n); for t =1:T; vShare_est[marketid.==t,:]=temp[t]; end;
	@time temp_par=pmap(t->demand(vDelta[:,:,t],vNLparam,t),[i for i=1:T]);
	vShare_est_par=zeros(n); for t =1:T; vShare_est_par[marketid.==t,:]=temp_par[t]; end;
	vShare_est==vShare_est_par

# End of the Code;
