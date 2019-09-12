%Dino Sejdinovic, 2013
%D. Sejdinovic, A. Gretton and W. Bergsma.  A KERNEL TEST FOR THREE-VARIABLE INTERACTIONS, 2013.
%runs three-variable interaction tests based on the file located in data/

%input format:
%--filename: name of the file (must be in 'data/' subfolder)
%-------data in the file needs to consist of the following three matrices:
%-------xx (size n x dx)
%-------yy (size n x dy)
%-------zz (size n x dz)
%-------n is the number of samples, dx, dy, dz are respective dimensions
%--save_results [OPTIONAL]: boolean specifying whether the results should
%---------------------------be saved; DEFAULT: false

%output format
%--------param: the structure of all parameters used in the test
%--------pvalues: the structure of all pvalues obtained in the tests,
%--------rejectionrate: the rates at which null is rejected in the tests
%-----------------------(averages over number of trials)
%--------pvalues and rejectionrate may contain some or all of substructures
%-----------------------TwoVar, Lancaster, Total, CI
%--------access individual pvalue vectors with, e.g., pvalues.TwoVar.XY_Z

%CI code is from [Zhang et al, 2011] ---- code available at:
%--------http://people.tuebingen.mpg.de/kzhang/KCI-test.zip
%we use the function indtest_new which can be found in '../KCI_test/'

%function [param,pvalues,rejectionrate]=experiment_from_file(filename,save_results)
function [param,pvalues,rejectionrate]=experiment_from_file(xx,yy,zz,filename,save_results)

if nargin==1
    save_results=false;
end

%init_seed

%load(['data/' filename]);
param.num_shuffles=2000;    %number of shuffles in the permutation test
param.use_median=true;      %using median heuristic for kernel bandwidth
param.alpha=0.05;           %confidence level
param.num_trials=1;       %number of trials to be run
param.noise_std=0;          %standard deviation of noise added
%param.subsample_size=796;    %subsample size (must be <n)
param.subsample_size=size(xx,1);
%specify which tests should be run below:
param.run_lancaster=true;
param.run_twovar=true;
param.run_total=false;
param.run_ci=false;

%if CI tests to be run, add KCI code to the path
if param.run_ci
    addpath(genpath('../KCI_test'));
end

%error raised if subsample size is larger than the sample size in the file 
% if param.subsample_size>size(xx,1);
%     fprintf '\nError! Subsample size larger than the sample size!\n';
%     return;
% end

%standardize if required
% muxx=mean(xx); stdxx=std(xx); xx=(xx-muxx)./stdxx;
% muyy=mean(yy); stdyy=std(yy); yy=(yy-muyy)./stdyy;
% muzz=mean(zz); stdzz=std(zz); zz=(zz-muzz)./stdzz;


for jj=1:param.num_trials
    fprintf(strcat('Trial', num2str(jj),'\n'));
        
    %generate a new subsample
    %tau=randperm(size(xx,1));
    %xxn=xx(tau(1:param.subsample_size),:);
    %yyn=yy(tau(1:param.subsample_size),:);
    %zzn=zz(tau(1:param.subsample_size),:);
    
    %add noise to check robustness if required
    %if param.noise_std>0
     %   xxn=xxn+param.noise_std*randn(size(xxn));
      %  yyn=yyn+param.noise_std*randn(size(yyn));
      %  zzn=zzn+param.noise_std*randn(size(zzn));
    %end

    if param.run_twovar
    results(jj).TwoVar   =  all_two_var_tests(xx,yy,zz,param);
    end
  
    if param.run_lancaster
    results(jj).Lancaster = lancaster_tests(xx,yy,zz,param);
    end
    
    if param.run_total
    results(jj).Total = total_independence_test(xx,yy,zz,param);
    end
    
    if param.run_ci
    %calling the UAI2011 CI test for "X independent of Y given Z"
    results(jj).CI.XYgivenZ.pvalue = indtest_new(xx,yy,zz,[]);
    end
   
    
end
[pvalues,rejectionrate] = reshape_results(results,param);

if save_results
filename=['FromFile' filename]; 
save(['results/' filename],'param','pvalues','rejectionrate');
end