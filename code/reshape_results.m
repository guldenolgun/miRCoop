%Dino Sejdinovic, 2013
%D. Sejdinovic, A. Gretton and W. Bergsma.  A KERNEL TEST FOR THREE-VARIABLE INTERACTIONS, 2013.
%helper function---reshapes results
function [pvalues,rejectionrate] = reshape_results(results,param)
%load(['results/' filename]);
    
foo=fieldnames(results);
for ii=1:param.num_trials
    for jj=1:size(foo,1)
        bar=fieldnames( results(ii).(foo{jj}) );
        for ll=1:size(bar,1)
            pvalues.(foo{jj}).(bar{ll})(ii)=results(ii).(foo{jj}).(bar{ll}).pvalue;
        end
    end
    
end

clear foo
clear bar

foo=fieldnames(pvalues);
for jj=1:size(foo,1)
    bar=fieldnames( pvalues.(foo{jj}) );
    for ll=1:size(bar,1)
        rejectionrate.(foo{jj}).(bar{ll})=mean( pvalues.(foo{jj}).(bar{ll}) < param.alpha );
    end
end



%Lancaster: Holm-Bonferroni correction for the 'Factorization' hypothesis
if param.run_lancaster
Lancaster_rej = true(param.num_trials,1);
for ii=1:param.num_trials
    foo=[pvalues.Lancaster.YZ_X(ii) pvalues.Lancaster.XZ_Y(ii) pvalues.Lancaster.XY_Z(ii)];
    foosort=sort(foo);   
    for jj=1:3
        if foosort(jj)>param.alpha/(3-jj+1)
            Lancaster_rej(ii)=false;
        end
    end    
end
rejectionrate.Lancaster.Factor = mean(Lancaster_rej);
clear foo
clear bar
end

%TwoVar: Holm-Bonferroni correction for the 'Factorization' hypothesis
if param.run_twovar
TwoVar_rej    = true(param.num_trials,1);
for ii=1:param.num_trials
    bar=[pvalues.TwoVar.YZ_X(ii) pvalues.TwoVar.XZ_Y(ii) pvalues.TwoVar.XY_Z(ii)];
    barsort=sort(bar);
    for jj=1:3        
        if barsort(jj)>param.alpha/(3-jj+1)
            TwoVar_rej(ii)=false;
        end
    end    
end
rejectionrate.TwoVar.Factor    = mean(TwoVar_rej);
clear foo
clear bar
end




end
