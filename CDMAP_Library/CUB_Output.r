#Codon Usage Bias Output Script
#purpose: this script generates a text output for each codon and the distribution of triplets
#for each codon

matrix_obj <- input_matrix

col_names <- c("Codon", "Count")

Phe_mat <- matrix(0L, nrow = 0, ncol =2)
colnames(Phe_mat) <- col_names
TTT <- c("TTT", matrix_obj[1,1])
TTC <- c("TTC", matrix_obj[3,1])
Phe_mat <- rbind(Phe_mat, TTT)
Phe_mat <- rbind(Phe_mat, TTC)
frame <- data.frame(Phe_mat)
Phe_mat <- Phe_mat[order(frame$Count, decreasing = FALSE),]
rownames(Phe_mat) <- NULL

Leu_mat <- matrix(0L, nrow = 0, ncol =2)
colnames(Leu_mat) <- col_names
TTA <- c("TTA", matrix_obj[4,1])
TTG <- c("TTG", matrix_obj[2,1])
CTT <- c("CTT", matrix_obj[9,1])
CTC <- c("CTC", matrix_obj[11,1])
CTA <- c("CTA", matrix_obj[12,1])
CTG <- c("CTG", matrix_obj[10,1])
Leu_mat <- rbind(Leu_mat, TTA)
Leu_mat <- rbind(Leu_mat, TTG)
Leu_mat <- rbind(Leu_mat, CTT)
Leu_mat <- rbind(Leu_mat, CTC)
Leu_mat <- rbind(Leu_mat, CTA)
Leu_mat <- rbind(Leu_mat, CTG)
frame <- data.frame(Leu_mat)
Leu_mat <- Leu_mat[order(frame$Count, decreasing = FALSE),]
rownames(Leu_mat) <- NULL

Ile_mat <- matrix(0L, nrow = 0, ncol =2)
colnames(Ile_mat) <- col_names
ATT <- c("ATT", matrix_obj[13,1])
ATC <- c("ATC", matrix_obj[15,1])
ATA <- c("ATA", matrix_obj[16,1])
Ile_mat <- rbind(Ile_mat, ATT)
Ile_mat <- rbind(Ile_mat, ATC)
Ile_mat <- rbind(Ile_mat, ATA)
frame <- data.frame(Ile_mat)
Ile_mat <- Ile_mat[order(frame$Count, decreasing = FALSE),]
rownames(Ile_mat) <- NULL

Met_mat <- matrix(0L, nrow = 0, ncol =2)
colnames(Met_mat) <- col_names
ATG <- c("TTT", matrix_obj[14,1])
Met_mat <- rbind(Met_mat, ATG)
frame <- data.frame(Met_mat)
Met_mat <- Met_mat[order(frame$Count, decreasing = FALSE),]
rownames(Met_mat) <- NULL

Val_mat <- matrix(0L, nrow = 0, ncol =2)
colnames(Val_mat) <- col_names
GTT <- c("GTT", matrix_obj[5,1])
GTC <- c("GTC", matrix_obj[7,1])
GTA <- c("GTA", matrix_obj[8,1])
GTG <- c("GTG", matrix_obj[6,1])
Val_mat <- rbind(Val_mat, GTT)
Val_mat <- rbind(Val_mat, GTC)
Val_mat <- rbind(Val_mat, GTA)
Val_mat <- rbind(Val_mat, GTG)
frame <- data.frame(Val_mat)
Val_mat <- Val_mat[order(frame$Count, decreasing = FALSE),]
rownames(Val_mat) <- NULL

Ser_mat <- matrix(0L, nrow = 0, ncol =2)
colnames(Ser_mat) <- col_names
TCT <- c("TCT", matrix_obj[1,3])
TCC <- c("TCC", matrix_obj[3,3])
TCA <- c("TCA", matrix_obj[4,3])
TCG <- c("TCG", matrix_obj[2,3])
AGT <- c("AGT", matrix_obj[13,2])
AGC <- c("AGC", matrix_obj[15,2])
Ser_mat <- rbind(Ser_mat, TCT)
Ser_mat <- rbind(Ser_mat, TCC)
Ser_mat <- rbind(Ser_mat, TCA)
Ser_mat <- rbind(Ser_mat, TCG)
Ser_mat <- rbind(Ser_mat, AGT)
Ser_mat <- rbind(Ser_mat, AGC)
frame <- data.frame(Ser_mat)
Ser_mat <- Ser_mat[order(frame$Count, decreasing = FALSE),]
rownames(Ser_mat) <- NULL

Pro_mat <- matrix(0L, nrow = 0, ncol =2)
colnames(Pro_mat) <- col_names
CCT <- c("CCT", matrix_obj[9,3])
CCC <- c("CCC", matrix_obj[11,3])
CCA <- c("CCA", matrix_obj[12,3])
CCG <- c("CCG", matrix_obj[10,3])
Pro_mat <- rbind(Pro_mat, CCT)
Pro_mat <- rbind(Pro_mat, CCC)
Pro_mat <- rbind(Pro_mat, CCA)
Pro_mat <- rbind(Pro_mat, CCG)
frame <- data.frame(Pro_mat)
Pro_mat <- Pro_mat[order(frame$Count, decreasing = FALSE),]
rownames(Pro_mat) <- NULL

Thr_mat <- matrix(0L, nrow = 0, ncol =2)
colnames(Thr_mat) <- col_names
ACT <- c("ACT", matrix_obj[13,3])
ACC <- c("ACC", matrix_obj[15,3])
ACA <- c("ACA", matrix_obj[16,3])
ACG <- c("ACG", matrix_obj[14,3])
Thr_mat <- rbind(Thr_mat, ACT)
Thr_mat <- rbind(Thr_mat, ACC)
Thr_mat <- rbind(Thr_mat, ACA)
Thr_mat <- rbind(Thr_mat, ACG)
frame <- data.frame(Thr_mat)
Thr_mat <- Thr_mat[order(frame$Count, decreasing = FALSE),]
rownames(Thr_mat) <- NULL

Ala_mat <- matrix(0L, nrow = 0, ncol =2)
colnames(Ala_mat) <- col_names
GCT <- c("GCT", matrix_obj[5,3])
GCC <- c("GCC", matrix_obj[7,3])
GCA <- c("GCA", matrix_obj[8,3])
GCG <- c("GCG", matrix_obj[6,3])
Ala_mat <- rbind(Ala_mat, GCT)
Ala_mat <- rbind(Ala_mat, GCC)
Ala_mat <- rbind(Ala_mat, GCA)
Ala_mat <- rbind(Ala_mat, GCG)
frame <- data.frame(Ala_mat)
Ala_mat <- Ala_mat[order(frame$Count, decreasing = FALSE),]
rownames(Ala_mat) <- NULL

Tyr_mat <- matrix(0L, nrow = 0, ncol =2)
colnames(Tyr_mat) <- col_names
TAT <- c("TAT", matrix_obj[1,4])
TAC <- c("TAC", matrix_obj[3,4])
Tyr_mat <- rbind(Tyr_mat, TAT)
Tyr_mat <- rbind(Tyr_mat, TAC)
frame <- data.frame(Tyr_mat)
Tyr_mat <- Tyr_mat[order(frame$Count, decreasing = FALSE),]
rownames(Tyr_mat) <- NULL

His_mat <- matrix(0L, nrow = 0, ncol =2)
colnames(His_mat) <- col_names
CAT <- c("CAT", matrix_obj[9,4])
CAC <- c("CAC", matrix_obj[11,4])
His_mat <- rbind(His_mat, CAT)
His_mat <- rbind(His_mat, CAC)
frame <- data.frame(His_mat)
His_mat <- His_mat[order(frame$Count, decreasing = FALSE),]
rownames(His_mat) <- NULL

Gln_mat <- matrix(0L, nrow = 0, ncol =2)
colnames(Gln_mat) <- col_names
CAA <- c("CAA", matrix_obj[12,4])
CAG <- c("CAG", matrix_obj[10,4])
Gln_mat <- rbind(His_mat, CAA)
Gln_mat <- rbind(His_mat, CAG)
frame <- data.frame(Gln_mat)
Gln_mat <- Gln_mat[order(frame$Count, decreasing = FALSE),]
rownames(Gln_mat) <- NULL

STOP_mat <- matrix(0L, nrow = 0, ncol =2)
colnames(STOP_mat) <- col_names
TAA <- c("TAA", matrix_obj[4,4])
TAG <- c("TAG", matrix_obj[2,4])
TGA <- c("TGA", matrix_obj[4,2])
STOP_mat <- rbind(STOP_mat, TAA)
STOP_mat <- rbind(STOP_mat, TAG)
STOP_mat <- rbind(STOP_mat, TGA)
frame <- data.frame(STOP_mat)
STOP_mat <- STOP_mat[order(frame$Count, decreasing = FALSE),]
rownames(STOP_mat) <- NULL

Asn_mat <- matrix(0L, nrow = 0, ncol =2)
colnames(Asn_mat) <- col_names
AAT <- c("AAT", matrix_obj[13,4])
AAC <- c("AAC", matrix_obj[15,4])
Asn_mat <- rbind(Asn_mat, AAT)
Asn_mat <- rbind(Asn_mat, AAC)
frame <- data.frame(Asn_mat)
Asn_mat <- Asn_mat[order(frame$Count, decreasing = FALSE),]
rownames(Asn_mat) <- NULL

Lys_mat <- matrix(0L, nrow = 0, ncol =2)
colnames(Lys_mat) <- col_names
AAA <- c("AAA", matrix_obj[16,4])
AAG <- c("AAG", matrix_obj[14,4])
Lys_mat <- rbind(Lys_mat, AAA)
Lys_mat <- rbind(Lys_mat, AAG)
frame <- data.frame(Lys_mat)
Lys_mat <- Asn_mat[order(frame$Count, decreasing = FALSE),]
rownames(Lys_mat) <- NULL

Asp_mat <- matrix(0L, nrow = 0, ncol =2)
colnames(Asp_mat) <- col_names
GAT <- c("GAT", matrix_obj[5,4])
GAC <- c("GAC", matrix_obj[7,4])
Asp_mat <- rbind(Asp_mat, GAT)
Asp_mat <- rbind(Asp_mat, GAC)
frame <- data.frame(Asp_mat)
Asp_mat <- Asp_mat[order(frame$Count, decreasing = FALSE),]
rownames(Asp_mat) <- NULL

Glu_mat <- matrix(0L, nrow = 0, ncol =2)
colnames(Glu_mat) <- col_names
GAA <- c("GAA", matrix_obj[8,4])
GAG <- c("GAG", matrix_obj[6,4])
Glu_mat <- rbind(Glu_mat, GAA)
Glu_mat <- rbind(Glu_mat, GAG)
frame <- data.frame(Glu_mat)
Glu_mat <- Glu_mat[order(frame$Count, decreasing = FALSE),]
rownames(Glu_mat) <- NULL

Cys_mat <- matrix(0L, nrow = 0, ncol =2)
colnames(Cys_mat) <- col_names
TGT <- c("TGT", matrix_obj[1,2])
TGC <- c("TGC", matrix_obj[3,2])
Cys_mat <- rbind(Cys_mat, TGT)
Cys_mat <- rbind(Cys_mat, TGC)
frame <- data.frame(Glu_mat)
Cys_mat <- Cys_mat[order(frame$Count, decreasing = FALSE),]
rownames(Cys_mat) <- NULL

Trp_mat <- matrix(0L, nrow = 0, ncol =2)
colnames(Trp_mat) <- col_names
TGG <- c("TGG", matrix_obj[2,2])
Trp_mat <- rbind(Trp_mat, TGG)
frame <- data.frame(Trp_mat)
Trp_mat <- Trp_mat[order(frame$Count, decreasing = FALSE),]
rownames(Trp_mat) <- NULL

Arg_mat <- matrix(0L, nrow = 0, ncol =2)
colnames(Arg_mat) <- col_names
CGT <- c("CGT", matrix_obj[9,2])
CGC <- c("CGC", matrix_obj[11,2])
CGA <- c("CGA", matrix_obj[12,2])
CGG <- c("CGG", matrix_obj[10,2])
AGA <- c("AGA", matrix_obj[16,2])
AGG <- c("AGG", matrix_obj[14,2])
Arg_mat <- rbind(Arg_mat, CGT)
Arg_mat <- rbind(Arg_mat, CGC)
Arg_mat <- rbind(Arg_mat, CGA)
Arg_mat <- rbind(Arg_mat, CGG)
Arg_mat <- rbind(Arg_mat, AGA)
Arg_mat <- rbind(Arg_mat, AGG)
frame <- data.frame(Arg_mat)
Leu_mat <- Arg_mat[order(frame$Count, decreasing = FALSE),]
rownames(Arg_mat) <- NULL

Gly_mat <- matrix(0L, nrow = 0, ncol =2)
colnames(Gly_mat) <- col_names
GGT <- c("GGT", matrix_obj[5,2])
GGC <- c("GGC", matrix_obj[7,2])
GGA <- c("GGA", matrix_obj[6,2])
Gly_mat <- rbind(Gly_mat, GGT)
Gly_mat <- rbind(Gly_mat, GGC)
Gly_mat <- rbind(Gly_mat, GGA)
frame <- data.frame(Gly_mat)
STOP_mat <- Gly_mat[order(frame$Count, decreasing = FALSE),]
rownames(Gly_mat) <- NULL

setwd(output_path)
#starting the output report
#===================
CUB_file_name <- output_name
sink(file = CUB_file_name, type = c("output", "message"), split = TRUE)
print("Phenylalanine - Phe - F codon usage bias: ")
print("\n")
print(Phe_mat)

print("Leucine - Leu - L codon usage bias: ")
print("\n")
print(Leu_mat)

print("Isoleucine - Ile - I codon usage bias: ")
print("\n")
print(Ile_mat)

print("Methionine - Met - M codon usage bias: ")
print("\n")
print(Met_mat)

print("Valine - Val - V codon usage bias: ")
print("\n")
print(Val_mat)

print("Serine - Ser - S codon usage bias: ")
print("\n")
print(Ser_mat)

print("Proline - Pro - P codon usage bias: ")
print("\n")
print(Pro_mat)

print("Theronine - Thr - T codon usage bias: ")
print("\n")
print(Thr_mat)

print("Alanine - Ala - A codon usage bias: ")
print("\n")
print(Ala_mat)

print("Tyrosine - Tyr - Y codon usage bias: ")
print("\n")
print(Tyr_mat)

print("STOP- STOP - STOP codon usage bias: ")
print("\n")
print(STOP_mat)

print("Histodine - His - H codon usage bias: ")
print("\n")
print(His_mat)

print("Glutamine - Gln - Q codon usage bias: ")
print("\n")
print(Gln_mat)

print("Asparagine - Asn - N codon usage bias: ")
print("\n")
print(Asn_mat)

print("Lysine - Lys - K codon usage bias: ")
print("\n")
print(Lys_mat)

print("Aspartic Acid - Asp - D codon usage bias: ")
print("\n")
print(Asp_mat)

print("Glutamic Acid - Glu - E codon usage bias: ")
print("\n")
print(Glu_mat)

print("Cystine - Cys - C codon usage bias: ")
print("\n")
print(Cys_mat)

print("Tryptophan - Trp - W codon usage bias: ")
print("\n")
print(Trp_mat)

print("Arginine - Arg - R codon usage bias: ")
print("\n")
print(Arg_mat)

print("Glycine - Gly - G codon usage bias: ")
print("\n")
print(Gly_mat)
sink()