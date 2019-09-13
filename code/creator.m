% Gulden Olgun, 2019
% 
% Merge interaction results
% input:
%   i : First miRNA index in the miRNA expression data and gene name data
%   j : Second miRNA index in the miRNA expression data and gene name data
%   k : mRNA index in the mRNA expression data and gene name data
% output: 
%   tmp_reject : Reject values of the tests
%   tmp_pvalue : P-values for the tests
function [tmp_reject, tmp_pvalue] = creator(i,j,k,rejectionrate, pvalues)
    field = fieldnames(rejectionrate);
    tmp_reject=[i,j,k];

    for ii = 1:size(field,1)
        temp = getfield(rejectionrate,field{ii});
        fields = fieldnames(temp);

        for kk = 1:size(fields,1)
             tmp_reject = [tmp_reject getfield(temp,fields{kk})];
        end
    end

    tmp_pvalue = [i,j,k];
    field = fieldnames(pvalues);

    for ii = 1:size(field,1)
        temp = getfield(pvalues,field{ii});
        fields = fieldnames(temp);

        for kk = 1:size(fields,1)
            tmp_pvalue = [tmp_pvalue getfield(temp,fields{kk})];
        end
    end
end