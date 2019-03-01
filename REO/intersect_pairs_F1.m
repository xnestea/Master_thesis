load("JP_GP_3.mat")

JP1 = C1;
JP2 = C2;
JP3 = C3;

clear C1 C2 C3

load("TCGA_GP_7.mat");
C1_c = intersect(C1, JP3, "rows");
C2_c = intersect(C2, JP1, "rows");
C3_c = intersect(C3, JP2, "rows");

clear C1 C2 C3 ans JP1 JP2 JP3 

%save("consensus_rev")