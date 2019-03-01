%f_score_sort_foward(pair,gidP,expE,gidE,labelE,method)
load("cons_pairs_c1.mat")

pair = p1;
gidP = unique(p1);
clear p1

filename = "../../../TPM_TCGA.txt";
delimiterIn = '\t';
headerlinesIn = 1;
T = importdata(filename,delimiterIn,headerlinesIn);
expE = T.data;
clear filename headerlinesIn delimiterIn

rownames = T.textdata(2:length(T.textdata), :);
rownames = str2double(rownames);
gidE = rownames;
clear rownames T

filename = "../../../cluster_table_TCGA.txt";
delimiterIn = '\t';
L = importdata(filename,delimiterIn);
for ii = 1:length(L)
    if L(ii)== 3
       L(ii) = 1;
    else
       L(ii) = 0;
    end
end
clear filename delimiterIn

[t1, tf1, f1ff] = (pair,gidP,expE,gidE,L,1);

clear expE f1ff gidE gidP ii L pair

save("res_m_1_TCGA")