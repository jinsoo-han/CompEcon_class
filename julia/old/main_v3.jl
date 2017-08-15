# Main Code: BLP Estimation;

# 0. Introduction;
	# Julia version check;
		#VERSION
	# Number of Processors;
		addprocs();
	# Directories;
		#UserPath = "C:/Users/jinshan";
		@everywhere UserPath = "/Users/jinsoohan/";
		@everywhere Path = "/Dropbox/BLP/";
		@everywhere Path_Data = string(UserPath,Path,"data/");
		cd(string(UserPath,Path));
		pwd()
	# Load standard packages;
		@everywhere using DataFrames, Optim, Distributions;
	# Load user-written packages;
		@everywhere include(string(UserPath,Path,"/BLP_module.jl"));
		#using BLP_module
	#include(string(UserPath,Path,"/buildquadrature_v1.jl"));
		@everywhere include(string(UserPath,Path,"/buildconsumertype_v1.jl"));

# 1. Prep Data
	# Construct DataSet
	@everywhere srand(1);
	@everywhere include(string(UserPath,Path,"/builddata_v1.jl"));
	println([mean(vDelta), var(vDelta)]);
	println("Avg Outside Share: ",mean([1-sum(vShare[id[i]]) for i=1:T]));
	println("Tot Nb of Products: ",n);
	@everywhere mX=[ones(n) vX];
	@everywhere mIV=[mX mDiff_IV_dist];

# 2. Estimation: Distance Diff_IV;
  #declare GMM Parameters;
		@everywhere vParam=[sigma_eta]; vLParam=zeros(size(mX,2));
		vSE=zeros(size([vLParam, vParam]));
		chi2=zeros(1); pvalue=zeros(1); df=zeros(1);
		mJacobian=zeros(n); vB=ones(size(mX,2)+1); iConv=zeros(1);
	#Estimation;
		println("Begin Estimation")
		@time res=optimize(x->gmm_obj(vDelta0,[x],A), 0.0,5.0, iterations=1000, show_trace=true);
		vParam=[res.minimizer[1]];

		#gmm_blp(vDelta0,vParam,vLParam,vSE,chi2,pvalue,df,mJacobian,vB,iConv);
	#Store Parameters;
    vParam_Diff=vParam;
    if sum(int(isnan(vParam)))!=0; vLParam_Diff=NaN; else; vLParam_Diff=vLParam; end;
    vSE_Diff=vSE;
    chi2_Diff=chi2; pvalue_Diff=pvalue; df_Diff=df;
    mJacobian_Diff=mJacobian; vB_Diff=vB[end]; iConv_Diff=iConv;

# 3. Summary Report;
	#=
	println("Diff IV: est/se of linear parameter: ",[vLParam_Diff vSE_Diff[1] (vLParam_Diff./vSE_Diff[1])]);
	println("Diff IV: est/se of non-linear parameters: ",[vParam_Diff vSE_Diff[2] (vParam_Diff./vSE_Diff[2])]);
	println("Diff IV: J-test: ",[chi2_Diff pvalue_Diff]);
	println("BLP IV: est/se of linear parameter: ",[vLParam_BLP vSE_BLP[1] (vLParam_BLP./vSE_BLP[1])]);
	println("BLP IV: est/se of non-linear parameters: ",[vParam_BLP vSE_BLP[2] (vParam_BLP./vSE_BLP[2])]);
	println("BLP IV: J-test: ",[chi2_BLP pvalue_BLP]);
	println("BLP99 IV: est/se of linear parameter: ",[vLParam_BLP99 vSE_BLP99[1] (vLParam_BLP99./vSE_BLP99[1])]);
	println("BLP99 IV: est/se of non-linear parameters: ",[vParam_BLP99 vSE_BLP99[2] (vParam_BLP99./vSE_BLP99[2])]);
	println("BLP99 IV: J-test: ",[chi2_BLP99 pvalue_BLP99]);
	println("Model Stat1: Mean/VAR of True Delta: ", [round(mean(vDelta),4) round(var(vDelta),4)])
	println("Model Stat2: Correlation b/t X and True Delta: ", beta1)
	println("Model Stat3: Transportation Cost: ", lambda)
	println("Diff IV: est/se/sig of non-linear parameters: ",[round(vParam_Diff,4) round(vSE_Diff[2],4) round((vParam_Diff./vSE_Diff[2]),4)]);
	println("BLP IV: est/se/sig of non-linear parameters: ",[round(vParam_BLP,4) round(vSE_BLP[2],4) round((vParam_BLP./vSE_BLP[2]),4)]);
	println("BLP99 IV: est/se/sig of non-linear parameters: ",[round(vParam_BLP99,4) round(vSE_BLP99[2],4) round((vParam_BLP99./vSE_BLP99[2]),4)]);
	=#

# End of the Code;
