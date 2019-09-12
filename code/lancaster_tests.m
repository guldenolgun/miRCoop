%Dino Sejdinovic, 2013
%D. Sejdinovic, A. Gretton and W. Bergsma.  A KERNEL TEST FOR THREE-VARIABLE INTERACTIONS, 2013.

%---Tests using a Lancaster statistic

function Lancaster = lancaster_tests(xx,yy,zz,param)

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

%---centering
H=eye(n)-repmat(1/n,n,n);
Kc=H*K*H;
Lc=H*L*H;
Mc=H*M*H;

%---Lancaster statistic
V_stat=LancasterStatistic(Kc,Lc,Mc);

V_distTot = estimate_Lancaster_null(Kc,Lc,Mc,param, 'Tot');
Lancaster.Tot.pvalue=nnz(V_distTot > V_stat) / param.num_shuffles;
Lancaster.Tot.rejectNull   = Lancaster.Tot.pvalue   <  param.alpha;

V_distYZ_X  = estimate_Lancaster_null(Kc,Lc,Mc,param, 'YZ_X');
Lancaster.YZ_X.pvalue = nnz(V_distYZ_X    > V_stat)   / param.num_shuffles;
Lancaster.YZ_X.rejectNull   = Lancaster.YZ_X.pvalue   <  param.alpha;

V_distXZ_Y  = estimate_Lancaster_null(Kc,Lc,Mc,param, 'XZ_Y');
Lancaster.XZ_Y.pvalue = nnz(V_distXZ_Y    > V_stat)   / param.num_shuffles;
Lancaster.XZ_Y.rejectNull   = Lancaster.XZ_Y.pvalue   <  param.alpha;

V_distXY_Z  = estimate_Lancaster_null(Kc,Lc,Mc,param, 'XY_Z');
Lancaster.XY_Z.pvalue = nnz(V_distXY_Z    > V_stat)   / param.num_shuffles;
Lancaster.XY_Z.rejectNull   = Lancaster.XY_Z.pvalue   <  param.alpha;

end

function [V_stat] = LancasterStatistic(Kc,Lc,Mc)
%Hadamard product of three matrices
foo=Kc.*Lc.*Mc;
%average of all matrix elements
V_stat=mean(foo(:));
end


function [V_dist] = estimate_Lancaster_null(Kc,Lc,Mc,param,which_null)

n=size(Kc,1);
V_dist=zeros(1,param.num_shuffles);
for jj=1:param.num_shuffles
    
    p1=randperm(n);
    p2=randperm(n);
    
    %initialize to no shuffling
    Kshc=Kc;
    Lshc=Lc;
    Mshc=Mc;

    if strcmp(which_null,'Tot')
        Lshc=Lc(p1,p1);
        Mshc=Mc(p2,p2);
    elseif strcmp(which_null,'YZ_X')
        Kshc=Kc(p1,p1);
    elseif strcmp(which_null,'XZ_Y')
        Lshc=Lc(p1,p1);
    elseif strcmp(which_null,'XY_Z')
        Mshc=Mc(p1,p1);
    else
        fprintf('ERROR! Wrong input to Lancaster null estimation!\n')
        return
    end

  
V_dist(jj)=LancasterStatistic(Kshc,Lshc,Mshc);
    
end

end



