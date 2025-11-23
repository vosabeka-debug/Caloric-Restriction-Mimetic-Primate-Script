#Monkey snRNAseq pattern calculation -----by LCY
library(Seurat)
library(dplyr)
library(ggplot2)
library(stringr)
library(DoubletFinder)
library(future)
library(RColorBrewer)
library(patchwork)
library(ggalluvial)
library(ggplot2)
library(ggridges)
library(SoupX)
library(cowplot)
library(reshape)
library(reshape2)
library(ggpubr)
library(ggrepel)
library(monocle)
library(pheatmap)
library(reshape)
library(reshape2)
library(hdWGCNA)
library(glmnet)
library(ggpmisc)
set.seed(5)
#----------------------------------------------------------------------------------------------#
plan("multicore", workers = 25)
options(future.globals.maxSize = 400000 * 1024^2)#100000MB~=100G
#-------------------------------------#
#ts
ts <- "liver"
analysis_dir <- 'liuchengyu/'
qc_dir <- paste0(analysis_dir, 'qc/')

DEG_dir <- paste0(analysis_dir, 'DEG/')

srt <- readRDS(paste0(analysis_dir,'liver_celltype.rds'))
# srt refers to the seurat object obtained after integrating samples and performing cell-type annotation.
# The metadata file must include the following fields: sample name (sample), cell type (cell_abbr), age (age_number), age group (age), and sex(gender).
srt@meta.data$bcd <- rownames(srt@meta.data)
o_donors <- unique(srt@meta.data[srt@meta.data$age == "O",]$bcd)
train_id <- sample(o_donors, length(o_donors) / 2)
test_id  <- setdiff(o_donors, train_id)

srt$O_split <- srt$age
srt$O_split <- as.vector(srt$O_split)
srt@meta.data[srt$bcd %in% train_id,]$O_split <- "O_train"
srt@meta.data[srt$bcd %in% test_id, ]$O_split <- "O_test"


combination <- srt
#01_pat_ana------------------------------------------
wd <- paste0(analysis_dir)

work_wd <- paste(wd,"Age_dependent/",sep = "")
wd_Cluster_smoothline <- paste(work_wd,"Cluster_smoothline/",sep = "")
png_work_wd <- paste(work_wd,"png/",sep = "")
png_wd_Cluster_smoothline <- paste(wd_Cluster_smoothline,"png/",sep = "")

dir.create(work_wd)
dir.create(wd_Cluster_smoothline)
dir.create(png_wd_Cluster_smoothline)
dir.create(png_work_wd)
age_ls <-c('Y','M','O_train','CR')
##young---------------------------------------------------
age=1
cell_list <- list.files(DEG_dir)
cells <- c()
for (cl in cell_list) {
  mmm <- strsplit(cl,paste0('_'))[[1]][3]
  
  cells <- c(cells,mmm)
}
cells <- unique(sapply(cells, function(vec){strsplit(vec,'\\.')[[1]][1]} ))
args = intersect(unique(cells),unique(combination@meta.data$cell_abbr))

for (cell in args) {
  try({
    #cell <- 'Hep-Z2'
    object <- subset(combination,subset=((cell_abbr%in%cell) & (O_split%in%age_ls) ))#cell_abbr
    if (length(table(object@meta.data$O_split))>=4) {#all sample exists
      cell_num_s <- table(object$O_split)
      if ( (cell_num_s[[1]]>=4)&(cell_num_s[[2]]>=4)&(cell_num_s[[3]]>=4)&(cell_num_s[[4]]>=4) ){ #每个样本的细胞数大于4
        print(paste0("-------------- ",ts,'_',cell,"_all samples exist ------------------"))
        Idents(object) <- "date"
        DefaultAssay(object) <- 'RNA'
        obj_Y <- object
        #obj_O <- subset(object,O_split=='O')
        ##preparation of DEGs------------------------------------------------------
        if (length(grep(pattern = paste0('.*','\\+'), x = cell, value = TRUE))==1  ) {
          grepname <- paste0(strsplit(cell,paste0('\\+') )[[1]][1],'\\+',strsplit(cell,paste0('\\+') )[[1]][2])
        }else{grepname <- cell}
        deg_list <- grep(pattern = paste0(grepname), x = list.files(DEG_dir), value = TRUE)
        deg_all <- c()
        for (i in deg_list) {
          deg <- read.csv(paste0(DEG_dir,i))
          deg_all <- rbind(deg_all,deg)
        }
        deg_all <- deg_all[deg_all$p_val_adj<=0.05,]
        deg_all$drc <- ifelse(deg_all$avg_log2FC>0,'U','D')
        deg_all <- deg_all[(deg_all$drc=='U'&deg_all$pct.1>=0.3)|(deg_all$drc=='D'&deg_all$pct.2>=0.3),]
        
        age_select_gene <- unique(as.vector(deg_all$gene))
        print(paste0("--------------",ts,'_',cell, '_'," total DEG number is ",length(age_select_gene),"--------------"))
        if (length(age_select_gene)>=20) {
          
          object_matrix <- GetAssayData(object = obj_Y,assay="RNA",layer='data')[age_select_gene,]
          CellID <- c()
          Sample_ID <- c()
          time_list_Y <- age_ls[1:3]
          for (m in time_list_Y) {
            sub_sample <- subset(object@meta.data,O_split%in%m) %>% rownames(.) %>% as.vector()
            set.seed(176)
            sample_id <- sample(x=sub_sample,size=200,replace=T)
            Sample_ID <- c(Sample_ID,sample_id)
            cellid <- paste0(cell,"_YO_",m,"_",c(1:200))
            CellID <- c(CellID,cellid)
          }
          anno_cell_sample <- data.frame(Sample_ID,CellID)
          
          all_sample_matrix <- c()
          all_cell_num <- length(time_list_Y)*200
          for (i in 1:all_cell_num) {
            if (i==1) {
              all_sample_matrix <- object_matrix[,which(colnames(object_matrix)%in%Sample_ID[i])] %>% as.matrix()
            }else{
              sample_matrix <- object_matrix[,which(colnames(object_matrix)%in%Sample_ID[i])] %>% as.matrix()
              all_sample_matrix <- cbind(all_sample_matrix,sample_matrix)
            }
          }
          
          colnames(all_sample_matrix) <- CellID
          count_matrix_Y <- all_sample_matrix
          
          print(paste0("--------------",ts,'_',cell, " sample_matrix is ready","--------------"))
          
          meta <- object@meta.data[,c("orig.ident","O_split","cell_abbr","celltype")]
          
          meta$Sample_ID <- rownames(meta)##vessel_YD0_ATAGAGACAACGCCCA-1
          for (i in 1:all_cell_num) {
            if (i==1) {
              all_sample_meta <- meta[which(rownames(meta)%in%Sample_ID[i]),] %>% as.data.frame()
            }else{
              sample_meta <- meta[which(rownames(meta)%in%Sample_ID[i]),] %>% as.data.frame()
              all_sample_meta <- rbind(all_sample_meta,sample_meta)
            }
          }
          
          all_sample_meta <- merge(all_sample_meta, anno_cell_sample, by="Sample_ID") %>% .[!duplicated(.),]
          rownames(all_sample_meta) <- as.vector(all_sample_meta$CellID)
          all_sample_meta$CellID <- factor(all_sample_meta$CellID, levels = CellID)
          all_sample_meta <- all_sample_meta[order(all_sample_meta$CellID,decreasing = F),]
          
          
          pd_Y <- all_sample_meta #pd means ordered meta.data
          pd_Y$State <- 1
          pd_Y$branch <- 'branch_1'
          print(paste0("--------------",ts,'_',cell ," monocle-obj-Y created --------------"))
          
          ##creat CR_monocle object-----------------------------------------------
          obj_CR <- object
          object_matrix <- GetAssayData(object = obj_CR)[age_select_gene,]
          CellID <- c()
          Sample_ID <- c()
          time_list_CR <- age_ls[c(1,2,4)]
          for (m in time_list_CR) {
            sub_sample <- subset(object@meta.data,O_split%in%m) %>% rownames(.) %>% as.vector()
            set.seed(176)
            sample_id <- sample(x=sub_sample,size=200,replace=T)
            Sample_ID <- c(Sample_ID,sample_id)
            cellid <- paste0(cell,"_YCR_",m,"_",c(1:200))
            CellID <- c(CellID,cellid)
          }
          anno_cell_sample <- data.frame(Sample_ID,CellID)
          
          all_sample_matrix <- c()
          all_cell_num <- length(time_list_CR)*200
          for (i in 1:all_cell_num) {
            if (i==1) {
              all_sample_matrix <- object_matrix[,which(colnames(object_matrix)%in%Sample_ID[i])] %>% as.matrix()
            }else{
              sample_matrix <- object_matrix[,which(colnames(object_matrix)%in%Sample_ID[i])] %>% as.matrix()
              all_sample_matrix <- cbind(all_sample_matrix,sample_matrix)
            }
          }
          
          colnames(all_sample_matrix) <- CellID
          count_matrix_CR <- all_sample_matrix
          
          print(paste0("--------------",ts,'_',cell, " sample_matrix is ready","--------------"))
          
          meta <- object@meta.data[,c("orig.ident","O_split","cell_abbr","celltype")]
          
          meta$Sample_ID <- rownames(meta)##vessel_YD0_ATAGAGACAACGCCCA-1
          for (i in 1:all_cell_num) {
            if (i==1) {
              all_sample_meta <- meta[which(rownames(meta)%in%Sample_ID[i]),] %>% as.data.frame()
            }else{
              sample_meta <- meta[which(rownames(meta)%in%Sample_ID[i]),] %>% as.data.frame()
              all_sample_meta <- rbind(all_sample_meta,sample_meta)
            }
          }
          
          all_sample_meta <- merge(all_sample_meta, anno_cell_sample, by="Sample_ID") %>% .[!duplicated(.),]
          rownames(all_sample_meta) <- as.vector(all_sample_meta$CellID)
          all_sample_meta$CellID <- factor(all_sample_meta$CellID, levels = CellID)
          all_sample_meta <- all_sample_meta[order(all_sample_meta$CellID,decreasing = F),]
          
          pd_CR <- all_sample_meta
          pd_CR$State=2
          
          pd_CR$branch <- 'branch_2'
          
          print(paste0("--------------",ts,'_',cell ," monocle-obj-O created --------------"))
          
          
          pd <- rbind(pd_Y,pd_CR)
          count_matrix <- cbind(count_matrix_Y,count_matrix_CR)
          fd <- as.data.frame(rownames(count_matrix_CR))
          colnames(fd) <- "gene_short_name"
          rownames(fd) <- fd$gene_short_name
          
          #creat monocle object
          my_cds <- newCellDataSet(count_matrix,#as(count_matrix, "sparseMatrix")
                                   phenoData = new("AnnotatedDataFrame", data = pd),
                                   featureData = new("AnnotatedDataFrame", data = fd))
          
          print(paste0("-------------- ",ts,'_',cell,"-my_cds are loaded ","------------------"))
          
          
          DelayedArray:::set_verbose_block_processing(TRUE)
          DelayedArray:::set_verbose_block_processing(TRUE)
          ##normalization-------------------------------------------
          my_cds <- estimateSizeFactors(my_cds)
          #norm_method = 'log'
          #my_cds <- estimateDispersions(my_cds)
          norm_method <- "vstExprs"
          ##############################################---
          #judgement-------------------------------------------------------#
          judgement = tryCatch({
            #correct pipeline
            my_cds <- estimateDispersions(my_cds)
            
            print(paste0('----------------------',cell,'_estimateDispersions finished','----------------------'))
          },  error=function(e) e
          )
          if(inherits(judgement, 'error')) {
            print(paste0(cell,'_estimateSizeFactors_alternative'))
            norm_method <- "log"
          }
          ######################################---
          agetime <- seq(0,100,length.out = 600)
          pData(my_cds)$Agetime[1:600] <- agetime
          pData(my_cds)$Agetime[601:1200] <- agetime
          pData(my_cds)$Agetime <- as.numeric(pData(my_cds)$Agetime)
          print("--------------RDS my_cds_subset are saved ------------------")
          
          
          
          agetime <- seq(0,100,length.out = 600)
          
          DelayedArray:::set_verbose_block_processing(TRUE)
          DelayedArray:::set_verbose_block_processing(TRUE)
          
          cluster_rows = TRUE
          hclust_method = "ward.D2"
          hmcols = NULL
          add_annotation_row = NULL
          add_annotation_col = NULL
          show_rownames = FALSE
          use_gene_short_name = TRUE
          scale_max = 3
          scale_min = -3
          trend_formula = "~sm.ns(Agetime, df=3) * branch"
          return_heatmap = FALSE
          cores = 20
          
          newdataA <- data.frame(Pseudotime = seq(0, 100, length.out = 200),
                                 Branch = as.factor(unique(as.character(pData(my_cds)$branch))[1]
                                 )
          )
          
          newdataB <- data.frame(Pseudotime = seq(0, 100, length.out = 200),
                                 Branch = as.factor(unique(as.character(pData(my_cds)$branch))[2]
                                 )
          )
          judgement = tryCatch({
            #correct pipeline
            BranchAB_exprs <- genSmoothCurves(my_cds[, ], cores = 21,
                                              trend_formula = trend_formula, relative_expr = T,
                                              new_data = rbind(newdataA,newdataB))
            print(paste0('----------------------',cell,'_genSmoothCurves finished','----------------------'))
          },  error=function(e) e
          )
          if(inherits(judgement, 'error')) {
            print(paste0(cell,'_genSmoothCurves_error'))
            next}
          ############################################################################--
          
          BranchA_exprs <- BranchAB_exprs[, 1:600]
          BranchB_exprs <- BranchAB_exprs[, 601:1200]
          
          branch_states <- c('1','2')
          common_ancestor_cells <- row.names(pData(my_cds)[pData(my_cds)$State ==
                                                             setdiff(pData(my_cds)$State, branch_states), ])
          
          BranchP_num <- 200
          BranchA_num <- 1200
          BranchB_num <- BranchA_num
          
          if (norm_method == "vstExprs") {
            BranchA_exprs <- vstExprs(my_cds, expr_matrix = BranchA_exprs)
            BranchB_exprs <- vstExprs(my_cds, expr_matrix = BranchB_exprs)
            
          }else if (norm_method == "log") {
            BranchA_exprs <- log10(BranchA_exprs + 1)
            BranchB_exprs <- log10(BranchB_exprs + 1)
          }
          
          col_gap_ind <- 601
          heatmap_matrix <- cbind(BranchA_exprs,BranchB_exprs)
          
          heatmap_matrix = heatmap_matrix[!apply(heatmap_matrix, 1,
                                                 sd) == 0, ]
          
          age_select_gene <- as.data.frame(rownames(heatmap_matrix))
          colnames(age_select_gene) <- "gene"
          write.csv(age_select_gene, paste(work_wd,cell,"_age_DEG_genes_dyn.csv",sep = ""))
          num_clusters = ifelse(length(age_select_gene$gene)<=50,length(age_select_gene$gene),50)
          
          
          heatmap_matrix = Matrix::t(scale(Matrix::t(heatmap_matrix),
                                           center = TRUE))
          heatmap_matrix = heatmap_matrix[is.na(row.names(heatmap_matrix)) ==
                                            FALSE, ]
          heatmap_matrix[is.nan(heatmap_matrix)] = 0
          heatmap_matrix[heatmap_matrix > scale_max] = scale_max
          heatmap_matrix[heatmap_matrix < scale_min] = scale_min
          heatmap_matrix_ori <- heatmap_matrix
          heatmap_matrix <- heatmap_matrix[is.finite(heatmap_matrix[,
                                                                    1]) & is.finite(heatmap_matrix[, col_gap_ind]), ]
          row_dist <- as.dist((1 - cor(Matrix::t(heatmap_matrix)))/2)
          row_dist[is.na(row_dist)] <- 1
          exp_rng <- range(heatmap_matrix)
          bks <- seq(exp_rng[1] - 0.1, exp_rng[2] + 0.1, by = 0.1)
          if (is.null(hmcols)) {
            hmcols <- c(colorRampPalette(c('#1AB7FC', '#F1F1F1'))(30),
                        colorRampPalette(c('#F1F1F1', '#FF3F12'))(30))
            
          }
          
          
          ph <- pheatmap(heatmap_matrix, useRaster = T, cluster_cols = FALSE,
                         cluster_rows = TRUE, show_rownames = F, show_colnames = F,
                         clustering_distance_rows = row_dist, clustering_method = hclust_method,
                         cutree_rows = num_clusters, silent = TRUE, filename = NA,
                         breaks = bks, color = hmcols)
          annotation_row <- data.frame(Cluster = factor(cutree(ph$tree_row,
                                                               num_clusters)))
          annotation_row <- data.frame(group=factor(cutree(ph$tree_row,num_clusters)))
          colnames(heatmap_matrix) <- c(1:ncol(heatmap_matrix))
          
          
          
          annotation_col <- data.frame(row.names = c(1:ncol(heatmap_matrix)),
                                       `Stage` = c(rep(paste0('S',c(1,2,3)),each=200)),
                                       `Age`=c(rep(c('age','CR'),each=600))
          )
          colnames(annotation_col) <- c('Stage','Age')
          annotation_colors=list(`group`=colorRampPalette(brewer.pal(num_clusters,"Dark2"))(num_clusters),
                                 `Stage`=c('#d8d3cd','#94d6d4','#05717c'),
                                 `Age`=c('#3aa9d3','#e8861b')
          )
          names(annotation_colors$`group`)=c(1:num_clusters)
          names(annotation_colors$`Stage`)=c('S1','S2','S3')
          names(annotation_colors$`Age`)=c('age','CR')
          
          names(annotation_colors$`group`)=c(1:num_clusters)
          
          if (!is.null(add_annotation_col)) {
            annotation_col <- cbind(annotation_col, add_annotation_col[fData(cds[row.names(annotation_col),
            ])$gene_short_name, 1])
          }
          ph_res <- pheatmap(heatmap_matrix[, ],
                             useRaster = T,
                             cluster_cols = FALSE,
                             cluster_rows = cluster_rows,
                             show_rownames = F,
                             show_colnames = F,
                             clustering_distance_rows = row_dist,
                             clustering_method = hclust_method,
                             cutree_rows = num_clusters,
                             annotation_row = annotation_row,
                             annotation_col = annotation_col,
                             annotation_colors = annotation_colors,
                             annotation_names_row = T,
                             annotation_names_col = F,
                             treeheight_row = 0,
                             breaks = bks,
                             fontsize = 2,
                             color = hmcols,
                             border_color = NA,
                             silent = TRUE,
                             filename = paste(work_wd,cell,"_age_DEGs_heatmap_mapping","_MPY.pdf",sep = ""),
                             width=5,
                             height=7)
          hmcols <- c(colorRampPalette(c('#101C3C', '#085EA5'))(7),
                      colorRampPalette(c('#085EA5', '#F1F1F1'))(23),
                      colorRampPalette(c('#F1F1F1', '#CE5A36'))(23),
                      colorRampPalette(c('#CE5A36', '#912406'))(7)
          )
          print(paste0("--------------",ts,'_',cell," clustered heatmap is saved ------------------") )
          
          clusters <- cutree(ph_res$tree_row, k = num_clusters)
          clustering <- data.frame(clusters)
          clustering[,1] <- as.numeric(clustering[,1])
          colnames(clustering) <- "Gene_Clusters"
          clustering$Gene <- rownames(clustering)
          clustering <- clustering[order(clustering$Gene_Clusters,decreasing = F),]
          table(clustering$Gene_Clusters)
          
          write.csv(clustering,paste(work_wd,cell,"_branched_htmp_genes_mapping_dyn",".csv",sep = ""),row.names = F)
          
          write.csv(heatmap_matrix,paste(work_wd,cell,"_branched_htmp_matrix_mapping_dyn",".csv",sep = ""))
          
          row_gap <- c()
          gap_lth_t <- 0
          for (i in 1:num_clusters) {
            tmp <- subset(clustering,clustering$Gene_Clusters==i)
            gap_lth <- length(tmp$Gene_Clusters)
            gap_lth_t <- gap_lth_t+gap_lth
            row_gap <-c(row_gap,gap_lth_t)
            
          }
          
          heatmap_matrix2 <- heatmap_matrix
          heatmap_matrix2 <- heatmap_matrix2[rownames(clustering),,drop=FALSE]
          
          
          png(filename=paste(work_wd,cell,"_age_DEGs_heatmap_mapping","_clt_odr_dyn.png",sep = ""), width=3000, height=4200, res=300)
          print(pheatmap(heatmap_matrix2[, ],
                         useRaster = T,
                         cluster_cols = FALSE,
                         cluster_rows = F,
                         show_rownames = T,
                         show_colnames = F,
                         clustering_distance_rows = row_dist,
                         clustering_method = hclust_method,
                         cutree_rows = num_clusters,
                         annotation_row = annotation_row,
                         annotation_col = annotation_col,
                         annotation_colors = annotation_colors,
                         annotation_names_row = F,
                         annotation_names_col = F,
                         treeheight_row = 0,
                         breaks = bks,
                         fontsize = 8,
                         number_color='red',
                         fontsize_row=6,
                         gaps_row = row_gap,
                         color = hmcols,
                         border_color = NA,
                         silent = TRUE)
          )
          dev.off()
          
          
          ph_res_heatmap_rownames <- clustering
          for (i in 1:num_clusters) {
            #i=2
            newgene<-rownames(subset(ph_res_heatmap_rownames,Gene_Clusters==i))
            newdata <- heatmap_matrix[newgene,]
            if (length(newgene)>1) {
              newdata_Y <- newdata[,1:600]
              newdata_O <- newdata[,601:1200]
              #---------------#
              tmp_data_Y<-t(newdata_Y)
              
              myFUN<-function(x){#done-----
                df = data.frame(count = x, cell= rownames(tmp_data_Y))
                return(df)
              }
              line_data_Y = do.call(rbind,apply(tmp_data_Y,2,myFUN))
              line_data_Y$group = rep(colnames(tmp_data_Y),each=length(rownames(tmp_data_Y)))
              line_data_Y$cell = factor(line_data_Y$cell,levels = 1:600)
              line_data_Y$smooth = rep('smooth',length(rownames(line_data_Y)))
              line_data_Y$Age = rep('Y',length(rownames(line_data_Y)))
              line_data_Y$cell = as.numeric(line_data_Y$cell)
              #---------------#
              colnames(newdata_O) <- 1:600
              tmp_data_O<-t(newdata_O)
              
              myFUN<-function(x){#done-----
                df = data.frame(count = x, cell= rownames(tmp_data_O))
                return(df)
              }
              line_data_O = do.call(rbind,apply(tmp_data_O,2,myFUN))
              line_data_O$group = rep(colnames(tmp_data_O),each=length(rownames(tmp_data_O)))
              line_data_O$cell = factor(line_data_O$cell,levels = 1:600)
              line_data_O$smooth = rep('smooth',length(rownames(line_data_O)))
              line_data_O$Age = rep('O',length(rownames(line_data_O)))
              line_data_O$cell = as.numeric(line_data_O$cell)
            }else{
              newdata_Y <- newdata[1:600]
              newdata_O <- newdata[601:1200]
              tmp_data_Y <- as.data.frame(newdata_Y)
              myFUN<-function(x){#done-----
                df = data.frame(count = x, cell= rownames(tmp_data_Y))
                return(df)
              }
              line_data_Y = data.frame(count=tmp_data_Y[,1],cell=rownames(tmp_data_Y))
              line_data_Y$group = rep(colnames(tmp_data_Y),each=length(rownames(tmp_data_Y)))
              line_data_Y$cell = factor(line_data_Y$cell,levels = 1:600)
              line_data_Y$smooth = rep('smooth',length(rownames(line_data_Y)))
              line_data_Y$Age = rep('YMO',length(rownames(line_data_Y)))
              line_data_Y$cell = as.numeric(line_data_Y$cell)
              #-------------------------#
              names(newdata_O) <- 1:600
              tmp_data_O <- as.data.frame(newdata_O)
              myFUN<-function(x){#done-----
                df = data.frame(count = x, cell= rownames(tmp_data_O))
                return(df)
              }
              line_data_O = data.frame(count=tmp_data_O[,1],cell=rownames(tmp_data_O))
              line_data_O$group = rep(colnames(tmp_data_O),each=length(rownames(tmp_data_O)))
              line_data_O$cell = factor(line_data_O$cell,levels = 1:600)
              line_data_O$smooth = rep('smooth',length(rownames(line_data_O)))
              line_data_O$Age = rep('YMCR',length(rownames(line_data_O)))
              line_data_O$cell = as.numeric(line_data_O$cell)
              ####
            }#if (length(newgene)>1){}else
            
            line_data <- rbind(line_data_Y,line_data_O)
            line_data$symbol <- paste0(line_data$group,line_data$Age,sep='_')
            gene_num <- length(newgene)
            title <- paste0("Cluster ",i,"\n","n = ",length(newgene))
            ylable <- expression(atop("Z-score","Expression"))
            
            g<-ggplot(data = line_data,aes(x = cell,y = count, group = symbol))
            plot1 <- g +
              stat_smooth(aes(colour=Age),method = lm,
                          formula = y~poly(x,15),se=F,
                          #fill = "",
                          #color=NA,
                          size = 0.25,linetype=2, alpha = 0.05)+
              geom_smooth(aes(group=Age,colour=Age),
                          formula = y~poly(x,15),method = lm,
                          size = 1.8,alpha = 1)+
              scale_colour_manual(values =c('#EE3B3B','#6495ED'))+
              labs(x="Stage",y=ylable,title = title)+
              coord_cartesian(ylim=c(-3, 3))+
              scale_y_continuous(breaks=seq(-3, 3, 1))+ 
              scale_x_continuous(limits=c(1,600), breaks=seq(100,600, 200),label = c("stage1","stage2","stage3"))+ #淇敼鍧愭爣杞存樉绀哄埢????breaks=seq(璧峰鍊硷紝缁堟鍊硷紝闂撮???))  label  璁剧疆鏍囩鍚嶏紝change the axis labels
              geom_vline(xintercept=c(1,200,400,600),lty=4,col=c('#919191','#919191','#919191','#919191'),lwd=0.6)+
              theme_bw() +
              theme(panel.grid = element_blank()) +
              theme(plot.title = element_text(hjust = 0.5))+
              theme(axis.text.y = element_text(size=14, colour = "black"))+
              theme(axis.text.x = element_text(size=14, colour = "black"))
            ggsave(paste(wd_Cluster_smoothline,cell,'_mapping','Cluster_',i,"_smoothline_dyn.pdf",sep = ""), plot = plot1, width = 10, height = 7)
            ggsave(paste(png_wd_Cluster_smoothline,cell,'_mapping','Cluster_',i,"_smoothline_dyn.png",sep = ""), plot = plot1, width = 10, height = 7)
            print(paste0("--------------Cluster[",i,"]_smoothlines are saved ------------------"))
          }#for (i in 1:num_clusters)
          
        }else{
          print(paste0("-------------- ",ts,'_',cell," DEGs not enough to analysis --------------"))
        }#if (length(age_select_gene)>=21){} else{}
      }else{
        print(paste0("-------------- ",ts,'_',cell," sample not enough to analysis --------------"))
      }#if (length(table(object@meta.data$sample))==4) else{}
    }else{
      print(paste0("-------------- ",ts,'_',cell," sample missing --------------"))
    }#if 8 sample exists{} else{}
  })
  
}#for (cell in args)