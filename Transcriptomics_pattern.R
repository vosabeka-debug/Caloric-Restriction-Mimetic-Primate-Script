library(data.table)
library(Mfuzz)
library(dplyr)
library(ggpubr)
library(PupillometryR)
library(ggthemes)
library(plyr)
library(ggplot2)
library(pheatmap)
options(scipen = 10)

#########################
dat <- read.csv("YMOCR_mean_FPKM.csv")
rownames(dat) <- dat$name
dat <- dat[,-1]
dat <- dat[which(rowSums(dat) > 0),]
gene <- as.matrix(dat)
colnames(gene) <- c("Y-AL","M-AL","O-AL","O-CR")
#构建对象
mfuzz_class <- new('ExpressionSet',exprs = gene)
exprs=Biobase::exprs
#预处理缺失值或者异常值
mfuzz_class <- filter.NA(mfuzz_class, thres = 0.25)
mfuzz_class <- fill.NA(mfuzz_class, mode = 'mean')
mfuzz_class <- filter.std(mfuzz_class, min.std = 0)
#标准化数据
mfuzz_class <- standardise(mfuzz_class)

# for (num in c(2:100)) {
#   set.seed(123)
#   number <- num
#   cluster_num <- number
#   mfuzz_cluster <- mfuzz(mfuzz_class, c = cluster_num, m = mestimate(mfuzz_class))
#   setwd("G:\\result\\monkey_RNA_seq\\P2P2/tran_pattern\\")
#   pdf(paste("Aging_72","_",number,"_4_CR.pdf",sep = ""),width=12, height=8)
#   p <- mfuzz.plot2(mfuzz_class, cl = mfuzz_cluster, mfrow = c(4, 5), time.labels = colnames(gene),x11 = F)
#   print(p)
#   dev.off()
# }
set.seed(123)
number <- 16
cluster_num <- number
mfuzz_cluster <- mfuzz(mfuzz_class, c = cluster_num, m = mestimate(mfuzz_class))
setwd("G:\\result\\monkey_RNA_seq\\P2P2/tran_pattern\\")
pdf(paste("Aging_72","_",number,"_4_CR.pdf",sep = ""),width=12, height=8)
p <- mfuzz.plot2(mfuzz_class, cl = mfuzz_cluster, mfrow = c(4, 5), time.labels = colnames(gene),x11 = F)
print(p)
dev.off()
cluster_size <- mfuzz_cluster$size
names(cluster_size) <- 1:cluster_num
cluster_size
protein_cluster <- mfuzz_cluster$cluster
protein_standard <- mfuzz_class@assayData$exprs
protein_standard_cluster <- cbind(protein_standard[names(protein_cluster), ], protein_cluster)
write.csv(protein_standard_cluster, paste("Aging_72","_",number,"_mfuzz_standard_cluster_aging.csv",sep = ""))
dt <- protein_standard_cluster
dt$type1 <- ifelse(dt$protein_cluster==5|dt$protein_cluster==9|dt$protein_cluster==14,"up_down",
                           ifelse(dt$protein_cluster==10,"down_up",
                                  ifelse(dt$protein_cluster==13,"up_up",
                                         ifelse(dt$protein_cluster==2|dt$protein_cluster==4,"down_not",
                                                ifelse(dt$protein_cluster==11,"not_up",
                                                       ifelse(dt$protein_cluster==15,"up_not","others"))))))
dt$type2 <- ifelse(dt$type1=="up_down"|dt$type1=="down_up","rescue",
                ifelse(dt$type1=="up_up","pro-aging",
                       ifelse(dt$type1=="up_not"|dt$type1=="down_not","aging-unique",
                              ifelse(dt$type1=="not_up","CR-unique","others")
                              )))
write.csv(dt,paste("Aging_72","_",number,"_female_male_CR_all_anno.csv",sep = ""),row.names =F)
