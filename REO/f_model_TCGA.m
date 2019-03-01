%f_score_sort_foward(pair,gidP,expE,gidE,labelE,method)
cd /home/xnestea/Master/Data/REO_to_K/Identify_stable_gene_pair/Binomial_test
load("TCGA_7.mat")

pair = C1;
gidP = unique(C1);
clear C1 C2 C3

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
    if L(ii)== 1
       L(ii) = 1;
    else
       L(ii) = 0;
    end
end
clear filename delimiterIn

[p1, f1, f1ff] = f_score_sort_foward(pair,gidP,expE,gidE,L,1);

clear expE f1ff gidE gidP ii L pair

%%% 2cluster

load("TCGA_7.mat")

pair = C2;
gidP = unique(C2);
clear C1 C2 C3

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
    if L(ii)== 2
       L(ii) = 1;
    else
       L(ii) = 0;
    end
end
clear filename delimiterIn

[p2, f2, f1ff] = f_score_sort_foward(pair,gidP,expE,gidE,L,1);

clear expE f1ff gidE gidP ii L pair

%%% 3cluster

load("TCGA_7.mat")

pair = C3;
gidP = unique(C3);
clear C1 C2 C3

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

[p3, f3, f1ff] = f_score_sort_foward(pair,gidP,expE,gidE,L,1);

clear expE f1ff gidE gidP ii L pair

save("TCGA_pairs_f.mat")