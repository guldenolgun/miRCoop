% Gulden Olgun, 2019
% 
% Perform two-variable and three variable tests
%
% input:
%   mirna : miRNA expression data
%   mrna : mRNA expression data
%   potTriplet : Potential triplets that passed the expression filter (Step 2).
%   mrnaGenes : mRNA genes. It has be same order as in mRNA expression data.
%   mirnaGenes : miRNA genes. It has be same order as in miRNA expression data.
%
% output: 
%   all_reject : Reject values of the tests
%   all_pvalue : P-values for the tests
function [all_reject all_pvalue] = step3(mirna, mrna, potTriplet, mrnaGenes, mirnaGenes)
all_reject = [];
all_pvalue = [];
for i = 1:size(potTriplet,1)
        [param,pvalues,rejectionrate]=experiment_from_file(mirna(:,potTriplet(i,1)),mirna(:,potTriplet(i,2)),mrna(:,potTriplet(i,3)),'a',false);
        [tmp_reject, tmp_pvalue] = creator(potTriplet(i,1),potTriplet(i,2),potTriplet(i,3),rejectionrate, pvalues);
        all_reject = [all_reject; tmp_reject];
        all_pvalue = [all_pvalue; tmp_pvalue];
end

end