%Dino Sejdinovic, 2013
%D. Sejdinovic, A. Gretton and W. Bergsma.  A KERNEL TEST FOR THREE-VARIABLE INTERACTIONS, 2013.
%runs a total independence test using a total independence statistic for a
%three-variable dataset
function Total = total_independence_test(xx,yy,zz,param)

n=size(xx,1);

if param.use_median
    sigx=median_heur(xx);
    sigy=median_heur(yy);
    sigz=median_heur(zz);
else
    sigx=param.sigx;
    sigy=param.sigy;
    sigz=param.sigz;
end


K=GaussKern(xx,xx,sigx);
L=GaussKern(yy,yy,sigy);
M=GaussKern(zz,zz,sigz);

%---total independence statistic & null
V_stat=TotalIndependenceStatistic(K,L,M);
V_dist=estimate_Total_null(K,L,M,param);
Total.Tot.pvalue  = nnz(V_dist   > V_stat)   / param.num_shuffles;
Total.Tot.rejectNull = Total.Tot.pvalue   < param.alpha;


end


function [totalsum] = TotalIndependenceStatistic(K,L,M)
bar=K.*L.*M;
bar1 = mean(bar(:));
bar2 = -2*mean(mean(K).*mean(L).*mean(M));
bar3 = mean(K(:))*mean(L(:))*mean(M(:));
totalsum = bar1 + bar2 + bar3;
end

function [V_dist] = estimate_Total_null(K,L,M,param)
n=size(K,1);
V_dist=zeros(1,param.num_shuffles);
for jj=1:param.num_shuffles
    
    pY=randperm(n);
    pZ=randperm(n);
    
    %shuffling uncentred matrices 
    Lsh=L(pY,pY);
    Msh=M(pZ,pZ);
    
    V_dist(jj)=TotalIndependenceStatistic(K,Lsh,Msh);
    
end

end


