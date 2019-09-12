%Dino Sejdinovic, 2013
%D. Sejdinovic, A. Gretton and W. Bergsma.  A KERNEL TEST FOR THREE-VARIABLE INTERACTIONS, 2013.

%---standard kernel independence test


function [pvalue, rejectNull] = two_var_test(xx,yy,param)

n=size(xx,1);

if param.use_median
    sigx=median_heur(xx);
    sigy=median_heur(yy);
else
    sigx=param.sigx;
    sigy=param.sigy;
end


K=GaussKern(xx,xx,sigx);
L=GaussKern(yy,yy,sigy);

%---centering
H=eye(n)-repmat(1/n,n,n);
Kc=H*K*H;
Lc=H*L*H;

%---Independence V-statistic
foo=Kc.*Lc;
V_stat=mean(foo(:));


%---permutation test: estimate the null distribution
%---by repeatedly permuting the indices of Y

V_dist=zeros(1,param.num_shuffles);

for jj=1:param.num_shuffles
    
    p1=randperm(n);
    Lshc=Lc(p1,p1);
    
    foo=Kc.*Lshc;
    V_dist(jj)=mean(foo(:));
    
    
end
pvalue    = nnz(V_dist   > V_stat)   / param.num_shuffles;
rejectNull   = pvalue   <  param.alpha;
end

