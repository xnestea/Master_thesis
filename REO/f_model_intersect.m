%f_score_sort_foward(pair,gidP,expE,gidE,labelE,method)
%f_score_sort_foward(pair,gidP,expE,gidE,labelE,method)
load("consensus_rev.mat")

pair = C1_c;
gidP = unique(C1_c);
clear C1_c C2_c C3_c

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

%C 2
load("consensus_rev.mat")

pair = C2_c;
gidP = unique(C2_c);
clear C1_c C2_c C3_c

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

load("consensus_rev.mat")

pair = C3_c;
gidP = unique(C3_c);
clear C1_c C2_c C3_c

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

save("cons_pairs_f")