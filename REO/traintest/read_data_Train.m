filename = "../../../TPM_Merged_train_1.txt";
delimiterIn = '\t';
headerlinesIn = 1;
T = importdata(filename,delimiterIn,headerlinesIn);
rownames = T.textdata(2:length(T.textdata), :);
rownames = str2double(rownames);
str = string(T.textdata(1,1));
colnames = strsplit(str, '"');
colnames(cellfun('length',colnames)<3) = [];
T1 = T;
colnames1 = colnames;

filename = "../../../TPM_Merged_train_2.txt";
delimiterIn = '\t';
headerlinesIn = 1;
T = importdata(filename,delimiterIn,headerlinesIn);
rownames = T.textdata(2:length(T.textdata), :);
rownames = str2double(rownames);
str = string(T.textdata(1,1));
colnames = strsplit(str, '"');
colnames(cellfun('length',colnames)<3) = [];
T2 = T;
colnames2 = colnames;

filename = "../../../TPM_Merged_train_3.txt";
delimiterIn = '\t';
headerlinesIn = 1;
T = importdata(filename,delimiterIn,headerlinesIn);
rownames = T.textdata(2:length(T.textdata), :);
rownames = str2double(rownames);
str = string(T.textdata(1,1));
colnames = strsplit(str, '"');
colnames(cellfun('length',colnames)<3) = [];
T3 = T;
colnames3 = colnames;

clear colnames delimiterIn filename str headerlinesIn T

save("Train_data_REO")