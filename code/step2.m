% Gulden Olgun, 2019
% Expression filter for the potential miRNA-miRNA-mRNA interactions.
%
% input:
%   mirna : miRNA expression data
%   mrna : mRNA expression data
%   ta : miRNA-mRNA target data
%   mrnaGenes : mRNA genes. It has be same order as in mRNA expression data.
%   mirnaGenes : miRNA genes. It has be same order as in miRNA expression data.
%
% output: 
%   att : Potential miRNA-miRNA-mRNA and index of the genes after the
%   expression filter
function [att] = step2(mirna, mrna, ta, mrnaGenes, mirnaGenes)
att=[];
for i = 1:size(mrnaGenes,1)
    i
    ext_mi =  find(strcmp(ta(:,2),mrnaGenes(i)));
    if(ext_mi>0)
        t_mir = ta(ext_mi,1);
        m = [];
        for j = 1 : size(t_mir,1)
            m = [m;find(strcmp(mirnaGenes,t_mir(j)))];
        end
        
        for l = 1 : size(m,1)-1
            for k = l + 1 : size(m,1)
                m2 = quantile(mirna(:,m(l)),0.25);
                m1 = quantile(mirna(:,m(k)),0.25);
                m3 = quantile(mirna(:,m(l)),0.75);
                m4 = quantile(mirna(:,m(k)),0.75);
                sample3 = find(mirna(:,m(l)) < m2 & mirna(:,m(k)) < m1);
                sample2 = find(mirna(:,m(l)) > m3 & mirna(:,m(k)) > m4);
                %sampl = find((mirna(sample2,m(l)).* mirna(sample2,m(k)))>(mirna(sample2,m(l))+ mirna(sample2,m(k))));
                %sample2 = sample2(sampl,:);
                sampledu = find(mirna(:,m(l)) < m2 & mirna(:,m(k)) > m4);
                sampleud = find(mirna(:,m(l)) > m3 & mirna(:,m(k)) < m1);
                
                % if(median(mrna(sample2,i)) < median(mrna(sample3,i)) & median(mrna(sample2,i)) < median(mrna(sampleud,i)) & median(mrna(sample2,i)) < median(mrna(sampledu,i)) & median(mrna(sample3,i)) > median(mrna(sampleud,i)) & median(mrna(sample3,i)) > median(mrna(sampledu,i)))
                % if(median(mrna(sample2,i)) < median(mrna(sample3,i)) & median(mrna(sample2,i)) < median(mrna(sampleud,i)) & median(mrna(sample2,i)) < median(mrna(sampledu,i)))
                if(median(mrna(sample2,i)) < median(mrna(sample3,i)))
                    [pval h2] = ranksum(mrna(sample2,i), mrna(sample3,i), 'tail' , 'left');
                    if(pval <0.05 && m(l)~= m(k))
                        att = [att; m(l),m(k),i ];
                    end
                end
                % att = [att;m(l),m(k),i];
            end
        end
        
    end
end
end