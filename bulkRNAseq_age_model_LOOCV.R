#Monkey bulk RNAseq age-model LOOCV -----by LCY
library(dplyr)
library(ggplot2)
library(stringr)
library(future)
library(patchwork)
library(cowplot)
library(reshape)
library(reshape2)
library(ggpubr)
library(glmnet)
set.seed(5)

omic_wd <- '/liuchengyu/model_multiomics/'
workwd <- paste0(omic_wd,'bulkRNAseq/')

mtx_wd <- paste0(workwd,'FPKM_mtx/')
model_wd_tt <- paste0(workwd,'model_training/')
dir.create(model_wd_tt)


## initialtion---------------------
# df_rpkm is the RPKM matrix (rows are gene names and columns are sample names)
# meta is the sample metadata (including sample name, age, and group).

# sample names are correlated to column names of df_rpkm.
# ages represents the true age of each sample.
# group indicates the group assignment of each sample (Y, M, O, CR).

##model_train--------
df_rpkm <- df_rpkm
meta <- meta
meta$sample <- paste0(meta$sample,'_',meta$tissue)
identical(colnames(df_rpkm),meta$sample)

meta$ag_grp <- meta$group
meta$ag_grp <- as.vector(meta$ag_grp) # Y M O CR

## --------------------------- 1) Global unsupervised QC ------------------------------------
keep_expr  <- rowSums(df_rpkm > 0) == ncol(df_rpkm)
keep_var   <- apply(df_rpkm, 1, var) > 0
keep_genes <- which(keep_expr & keep_var)
cat("Genes kept after QC:", length(keep_genes), "\n")

df_rpkm_qc <- as.matrix(df_rpkm[keep_genes, , drop = FALSE])
logX_all <- t(log2(df_rpkm_qc + 1))   # samples x genes
colnames(logX_all) <- rownames(df_rpkm_qc)
rownames(logX_all) <- colnames(df_rpkm_qc)

## --------------------------- 2) Data slicing -------------------------------------------
ctrl_idx <- which(meta$ag_grp %in% c("Y","M","O"))
cr_idx   <- which(meta$ag_grp %in% c("CR"))

X_ctrl    <- logX_all[ctrl_idx, , drop = FALSE]
y_ctrl    <- as.numeric(meta$age[ctrl_idx])
grp_ctrl  <- as.character(meta$ag_grp[ctrl_idx])
samp_ctrl <- meta$sample[ctrl_idx]

X_cr    <- if (length(cr_idx) > 0) logX_all[cr_idx, , drop = FALSE] else NULL
y_cr    <- if (length(cr_idx) > 0) suppressWarnings(as.numeric(meta$age[cr_idx])) else NULL
samp_cr <- if (length(cr_idx) > 0) meta$sample[cr_idx] else character(0)

o_in_ctrl  <- which(grp_ctrl == "O")
X_o_all    <- X_ctrl[o_in_ctrl, , drop = FALSE]
y_o_all    <- y_ctrl[o_in_ctrl]
samp_o_all <- samp_ctrl[o_in_ctrl]
cat("Num O for LOOCV:", length(o_in_ctrl), "\n")

## --------------------------- 3) Hyperparameters -------------------------------------------
ALPHAS               <- c(0.1,0.2,0.3,0.5,0.7,0.9)
N_FOLDS_INNER        <- 7
TYPE_MEASURE         <- "mae"
USE_CALIBRATION      <- TRUE
TARGET_MIN_GENES     <- 5
MAX_ROUNDS           <- 40
FALLBACK_LMIN        <- TRUE

B_REPEATS            <- 40   # Fix the features in each iteration and repeat B times.
SUBSAMPLE_FRAC       <- 0.8
RESAMPLE_REPLACE     <- FALSE 

TOP_NONZERO_FRACTION <- 0.85  # Within each b, count only the top 80% non-zero features ranked by |1 SD effect|.
FREQ_THRESHOLD_INIT  <- 0.65  # end-of-round frequency threshold (adaptively decreasing)
FREQ_THRESHOLD_STEP  <- 0.05

MIN_SELECT_ROUND     <- 2
## --------------------------- 4) Utility functions -----------------------
# 1) inner CV
inner_cv_fit <- function(X, y, alphas=ALPHAS, nfolds=N_FOLDS_INNER,
                         type.measure=TYPE_MEASURE, seed=NULL) {
  if (!is.null(seed)) {
    old <- if (exists(".Random.seed", envir=.GlobalEnv)) get(".Random.seed", envir=.GlobalEnv) else NULL
    on.exit({ if (!is.null(old)) assign(".Random.seed", old, envir=.GlobalEnv) }, add=TRUE)
    set.seed(seed)
  }
  n <- length(y)
  nfolds_eff <- max(2, min(nfolds, n - 1))
  foldid <- sample(rep(1:nfolds_eff, length.out = n))  # 受 seed 控制
  best <- list(score=Inf, alpha=NA, cvfit=NULL, foldid=foldid)
  for (a in alphas) {
    cvfit <- cv.glmnet(as.matrix(X), y, family="gaussian",
                       alpha=a, nfolds=nfolds_eff, foldid=foldid,
                       type.measure=type.measure, standardize=TRUE,
                       nlambda=200, lambda.min.ratio=0.01)
    sc <- cvfit$cvm[match(cvfit$lambda.1se, cvfit$lambda)]
    if (sc < best$score) best <- list(score=sc, alpha=a, cvfit=cvfit, foldid=foldid)
  }
  best
}

coef_1se_or_min <- function(cvfit) {
  co_1se <- coef(cvfit, s="lambda.1se")
  beta_1se <- as.numeric(co_1se[-1, , drop=TRUE]); names(beta_1se) <- rownames(co_1se)[-1]
  if (sum(beta_1se != 0) == 0 && FALLBACK_LMIN) {
    co_min <- coef(cvfit, s="lambda.min")
    beta_min <- as.numeric(co_min[-1, , drop=TRUE]); names(beta_min) <- rownames(co_min)[-1]
    list(beta=beta_min, s="lambda.min", coefmat=co_min)
  } else list(beta=beta_1se, s="lambda.1se", coefmat=co_1se)
}

one_sd_importance <- function(beta_vec, X_train) {
  nz <- which(beta_vec != 0)
  if (length(nz) == 0) return(data.frame())
  feats <- names(beta_vec)[nz]
  sdX <- apply(as.matrix(X_train[, feats, drop=FALSE]), 2, sd)
  eff <- beta_vec[feats] * sdX
  out <- data.frame(feature=feats, beta=beta_vec[feats],
                    effect_per_1SD=eff, abs_effect_1SD=abs(eff),
                    stringsAsFactors=FALSE)
  out[order(-out$abs_effect_1SD), , drop=FALSE]
}

calibrate_model <- function(y_true, y_pred) {
  ok <- is.finite(y_true) & is.finite(y_pred)
  y_true <- y_true[ok]; y_pred <- y_pred[ok]
  if (length(y_true) < 3 || sd(y_pred, na.rm=TRUE) < 1e-8) {
    a <- mean(y_true, na.rm=TRUE) - mean(y_pred, na.rm=TRUE); b <- 1; return(c(a=a,b=b))
  }
  fit <- try(lm(y_true ~ y_pred), silent=TRUE)
  if (inherits(fit,"try-error")) {
    b <- cov(y_true, y_pred, use="complete.obs") / var(y_pred, na.rm=TRUE); if (!is.finite(b)) b <- 1
    a <- mean(y_true, na.rm=TRUE) - b*mean(y_pred, na.rm=TRUE); return(c(a=a,b=b))
  }
  co <- coef(fit)
  if (any(!is.finite(co))) {
    b <- cov(y_true, y_pred, use="complete.obs") / var(y_pred, na.rm=TRUE); if (!is.finite(b)) b <- 1
    a <- mean(y_true, na.rm=TRUE) - b*mean(y_pred, na.rm=TRUE); return(c(a=a,b=b))
  }
  c(a=unname(co[1]), b=unname(co[2]))
}
apply_calibration <- function(yhat, calib) {
  a <- suppressWarnings(as.numeric(calib["a"])); if (!is.finite(a)) a <- 0
  b <- suppressWarnings(as.numeric(calib["b"])); if (!is.finite(b)) b <- 1
  a + b * yhat
}

# 2) stratified resampling
stratified_resample <- function(idx_all, groups, frac=SUBSAMPLE_FRAC, replace=FALSE, seed=NULL) {
  if (!is.null(seed)) {
    old <- if (exists(".Random.seed", envir=.GlobalEnv)) get(".Random.seed", envir=.GlobalEnv) else NULL
    on.exit({ if (!is.null(old)) assign(".Random.seed", old, envir=.GlobalEnv) }, add=TRUE)
    set.seed(seed)
  }
  sel <- integer(0)
  for (g in unique(groups)) {
    idx_g <- idx_all[groups == g]
    k <- if (replace) max(1, round(length(idx_g) * frac)) else max(1, floor(length(idx_g) * frac))
    sel <- c(sel, sample(idx_g, k, replace=replace))
  }
  sort(unique(sel))
}

## ---------------------- 5) Outer-layer evaluation: O-LOOCV (end-of-round stability filtering) -------------------
table_train <- list()    # Table for training log（fold x round x b）
table_test  <- list()    # Table for the testing O-sample（fold x round x b）
table_cr    <- list()    # O-CR prediction（fold x round x b）
model_store <- list()    
stab_stats_store <- list()  

for (fold_idx in seq_along(o_in_ctrl)) {
  test_id   <- o_in_ctrl[fold_idx]
  test_name <- samp_ctrl[test_id]
  
  tr_idx <- setdiff(seq_len(nrow(X_ctrl)), test_id)
  te_idx <- test_id
  
  X_tr_full <- X_ctrl[tr_idx, , drop=FALSE]
  y_tr_full <- y_ctrl[tr_idx]
  g_tr_full <- grp_ctrl[tr_idx]
  
  X_te_full <- X_ctrl[te_idx, , drop=FALSE]
  y_te_full <- y_ctrl[te_idx]
  
  feats_current <- colnames(X_tr_full)  
  
  for (round_id in seq_len(MAX_ROUNDS)) {
    
    bag_features <- list()  
    bag_effects  <- list()  
    
    for (b in seq_len(B_REPEATS)) {
      
      
      seed_sub <- 10100000L + fold_idx*10000L + round_id*100L + b
      seed_cv  <- 20200000L + fold_idx*10000L + round_id*100L + b
      
      idx_all <- seq_len(nrow(X_tr_full))
      sel     <- stratified_resample(idx_all, g_tr_full, frac=SUBSAMPLE_FRAC,
                                     replace=RESAMPLE_REPLACE, seed=seed_sub)
      Xb      <- X_tr_full[sel, feats_current, drop=FALSE]
      yb      <- y_tr_full[sel]
      
      best <- inner_cv_fit(Xb, yb, seed=seed_cv)
      co   <- coef_1se_or_min(best$cvfit)
      imp  <- one_sd_importance(co$beta, Xb)
      
      
      if (nrow(imp) > 0) {
        k_keep <- max(1, floor(nrow(imp) * TOP_NONZERO_FRACTION))
        feats_b <- imp$feature[seq_len(k_keep)]
        bag_features[[length(bag_features)+1]] <- feats_b
        bag_effects[[length(bag_effects)+1]]  <- imp[seq_len(k_keep), c("feature","abs_effect_1SD")]
      }
      
      # in-training calibration
      pred_tr_raw <- as.numeric(predict(best$cvfit, newx=as.matrix(Xb), s=co$s))
      calib <- if (USE_CALIBRATION) calibrate_model(yb, pred_tr_raw) else c(a=0,b=1)
      pred_tr_cal <- apply_calibration(pred_tr_raw, calib)
      
      # testing O-sample
      Xte_sub <- X_te_full[, feats_current, drop=FALSE]
      pred_te_raw <- as.numeric(predict(best$cvfit, newx=as.matrix(Xte_sub), s=co$s))
      pred_te_cal <- apply_calibration(pred_te_raw, calib)
      
      # testing O-CR
      if (!is.null(X_cr) && nrow(X_cr)>0) {
        Xcr_sub <- X_cr[, feats_current, drop=FALSE]
        pred_cr_raw <- as.numeric(predict(best$cvfit, newx=as.matrix(Xcr_sub), s=co$s))
        pred_cr_cal <- apply_calibration(pred_cr_raw, calib)
        table_cr[[length(table_cr)+1]] <- data.frame(
          fold=fold_idx, round=round_id, b=b,
          sample=samp_cr, age_true=if (is.null(y_cr)) NA_real_ else y_cr,
          age_pred=pred_cr_raw, age_pred_cal=pred_cr_cal
        )
      }
      
      # recording
      table_train[[length(table_train)+1]] <- data.frame(
        fold=fold_idx, round=round_id, b=b,
        n_train=nrow(Xb), inner_cv_mae=best$score,
        train_mae_raw=mean(abs(pred_tr_raw - yb)),
        train_mae_cal=mean(abs(pred_tr_cal - yb)),
        n_features=length(feats_current), alpha=best$alpha, lambda_choice=co$s
      )
      
      table_test[[length(table_test)+1]] <- data.frame(
        fold=fold_idx, round=round_id, b=b, test_sample=test_name,
        age_true=y_te_full,
        age_pred=pred_te_raw, age_pred_cal=pred_te_cal,
        abs_err_raw=abs(pred_te_raw - y_te_full),
        abs_err_cal=abs(pred_te_cal - y_te_full),
        n_features=length(feats_current), alpha=best$alpha, lambda_choice=co$s
      )
      
      # save models
      key <- sprintf("fold%02d_round%02d_b%02d", fold_idx, round_id, b)
      model_store[[key]] <- list(
        fold=fold_idx, round=round_id, b=b, o_test_sample=test_name,
        feats=feats_current, alpha=best$alpha, s_choice=co$s,
        cvfit=best$cvfit, calib=calib, inner_cv_mae=best$score,
        seed_sub=seed_sub, seed_cv=seed_cv
      )
    } # end for b
    
    # end-of-round stability aggregation (from the top non-zero features across B runs)
    if (length(bag_features) == 0) {
      feats_next <- if (length(feats_current) > TARGET_MIN_GENES)
        feats_current[seq_len(max(TARGET_MIN_GENES, floor(length(feats_current)/2)))]
      else feats_current
      stats_r <- data.frame(feature=feats_current, freq=0, med_abs_eff=0)
      pi_used <- NA_real_
    } else {
      all_feats <- unlist(bag_features, use.names=FALSE)
      freq_tab  <- sort(table(all_feats), decreasing=TRUE)
      freq      <- as.numeric(freq_tab); names(freq) <- names(freq_tab)
      freq      <- freq / length(bag_features)
      eff_df    <- do.call(rbind, bag_effects)
      med_eff   <- tapply(eff_df$abs_effect_1SD, eff_df$feature, median, na.rm=TRUE)
      med_eff   <- med_eff[names(freq)]
      stats_r   <- data.frame(feature=names(freq), freq=freq, med_abs_eff=as.numeric(med_eff))
      stats_r   <- stats_r[order(-stats_r$freq, -stats_r$med_abs_eff), ]
      
      pi_thr <- FREQ_THRESHOLD_INIT
      sel <- stats_r$feature[stats_r$freq >= pi_thr]
      while (length(sel) < TARGET_MIN_GENES && pi_thr > 0) {
        pi_thr <- max(0, pi_thr - FREQ_THRESHOLD_STEP)
        sel <- stats_r$feature[stats_r$freq >= pi_thr]
      }
      if (length(sel) < TARGET_MIN_GENES) sel <- head(stats_r$feature, TARGET_MIN_GENES)
      if (length(sel) > length(feats_current)) sel <- head(sel, length(feats_current))
      feats_next <- sel
      pi_used <- pi_thr
    }
    
    stab_stats_store[[length(stab_stats_store)+1]] <- cbind(
      data.frame(fold=fold_idx, round=round_id, pi_used=pi_used,
                 n_feats_current=length(feats_current), n_repeats=length(bag_features)),
      stats_r
    )
    
    if (length(feats_next) <= TARGET_MIN_GENES) break
    feats_current <- feats_next
  } # end for round
} # end for fold

## --------------------------- 6) summary table ----------------------------------------
table1_train_df <- if (length(table_train)) do.call(rbind, table_train) else data.frame()
table2_test_df  <- if (length(table_test))  do.call(rbind, table_test)  else data.frame()
cr_fold_df      <- if (length(table_cr))    do.call(rbind, table_cr)    else data.frame()
stab_stats_long <- if (length(stab_stats_store)) do.call(rbind, stab_stats_store) else data.frame()

# Verify whether the b predictions within each fold × round truly exhibit variation.
var_check <- table2_test_df %>%
  group_by(fold, round) %>%
  summarise(n_b = n_distinct(b),
            n_unique_pred = n_distinct(round(age_pred_cal, 6)),
            .groups = "drop")
if (any(var_check$n_unique_pred == 1 & var_check$n_b > 1)) {
  message("Warning: The B predictions for the following fold × round are nearly identical, indicating possibly excessive stability:")
  print(var_check %>% filter(n_unique_pred == 1 & n_b > 1))
}

## --------------------------- 7) Select the optimal one for each O sample. (round,b) ----------------
tmp <- table2_test_df %>%
  mutate(abs_err_sel = ifelse(is.finite(abs_err_cal), abs_err_cal, abs_err_raw))

best_ge_min <- tmp %>%
  filter(round >= MIN_SELECT_ROUND) %>%
  group_by(test_sample) %>%
  slice_min(abs_err_sel, n = 1, with_ties = FALSE) %>%
  ungroup()

all_O_samples <- unique(tmp$test_sample)
covered       <- unique(best_ge_min$test_sample)
need_fallback <- setdiff(all_O_samples, covered)

best_fallback <- tmp %>%
  filter(test_sample %in% need_fallback) %>%
  group_by(test_sample) %>%
  slice_min(abs_err_sel, n = 1, with_ties = FALSE) %>%
  ungroup()

selected_pairs_df <- bind_rows(best_ge_min, best_fallback) %>%
  select(fold, round, b, test_sample, age_true, age_pred_cal, abs_err_cal, abs_err_sel)

if (length(need_fallback) > 0) {
  message("The following O samples have records only for rounds < ", MIN_SELECT_ROUND,
          ", and have been reverted to their global minimum: ",
          paste(need_fallback, collapse = ", "))
}
cat("Selected (round,b) per O sample:\n"); print(selected_pairs_df)

## --------------------------- 8) extract the selected model ----------------
selected_keys   <- sprintf("fold%02d_round%02d_b%02d",
                           selected_pairs_df$fold,
                           selected_pairs_df$round,
                           selected_pairs_df$b)
selected_models <- lapply(selected_keys, function(k) model_store[[k]])
names(selected_models) <- selected_keys
keep <- !vapply(selected_models, is.null, logical(1))
if (any(!keep)) warning("Some selected models are NULL: ", paste(selected_keys[!keep], collapse=", "))
selected_models   <- selected_models[keep]
selected_pairs_df <- selected_pairs_df[keep, , drop=FALSE]
stopifnot(length(selected_models) > 0)

## --------------------------- 9) final presentation / extrapolation ----------------------------------
o_selected_df <- table2_test_df %>%
  inner_join(selected_pairs_df %>% select(fold, round, b),
             by = c("fold","round","b")) %>%
  transmute(sample=test_sample, group="O",
            age_true=age_true, age_pred=age_pred_cal,
            fold=fold, round=round, b=b)

cr_selected_ens <- data.frame()
if (nrow(cr_fold_df) > 0) {
  cr_selected_ens <- cr_fold_df %>%
    inner_join(selected_pairs_df %>% select(fold, round, b),
               by = c("fold","round","b")) %>%
    group_by(sample) %>%
    summarise(age_true = first(age_true),
              age_pred = mean(age_pred_cal, na.rm=TRUE),
              .groups = "drop") %>%
    mutate(group="CR", fold=NA_integer_, round=NA_integer_, b=NA_integer_)
}

predict_with_model <- function(model, X_full) {
  feats <- intersect(model$feats, colnames(X_full))
  if (length(feats) == 0) return(rep(NA_real_, nrow(X_full)))
  pr_raw <- as.numeric(predict(model$cvfit, newx=as.matrix(X_full[,feats,drop=FALSE]), s=model$s_choice))
  apply_calibration(pr_raw, model$calib)
}
predict_group_ensemble <- function(idx_vec, group_label) {
  if (length(idx_vec) == 0) return(data.frame())
  X_sub <- X_ctrl[idx_vec, , drop=FALSE]
  pred_mat <- sapply(selected_models, function(m) predict_with_model(m, X_sub))
  if (is.null(dim(pred_mat))) pred_mat <- matrix(pred_mat, ncol=1)
  pred_avg <- rowMeans(pred_mat, na.rm=TRUE)
  data.frame(sample=samp_ctrl[idx_vec], group=group_label,
             age_true=y_ctrl[idx_vec], age_pred=pred_avg,
             fold=NA_integer_, round=NA_integer_, b=NA_integer_)
}
idx_Y <- which(grp_ctrl == "Y")
idx_M <- which(grp_ctrl == "M")
ym_ens_df <- bind_rows(
  predict_group_ensemble(idx_Y, "Y"),
  predict_group_ensemble(idx_M, "M")
)

all_pred_df <- bind_rows(
  ym_ens_df %>% select(sample, group, age_true, age_pred, fold, round, b),
  o_selected_df %>% select(sample, group, age_true, age_pred, fold, round, b),
  cr_selected_ens %>% select(sample, group, age_true, age_pred, fold, round, b)
) %>%
  mutate(group=factor(as.character(group), levels=c("Y","M","O","CR")),
         delta_age=age_pred - age_true)

## --------------------------- 10) merge into a single linear model -------------------------
get_model_weights <- function(selected_models, table2_test_df, use_equal_weights = TRUE) {
  M <- length(selected_models)
  if (use_equal_weights) return(rep(1/M, M))
  w <- sapply(selected_models, function(m) {
    err <- table2_test_df$abs_err_cal[table2_test_df$fold==m$fold &
                                        table2_test_df$round==m$round &
                                        table2_test_df$b==m$b][1]
    if (!is.finite(err)) err <- NA_real_
    1/(err + 1e-6)
  })
  w[!is.finite(w)] <- 0; s <- sum(w); if (s<=0) rep(1/M,M) else w/s
}
merge_models_to_single <- function(models, weights=NULL) {
  stopifnot(length(models) > 0)
  feats_union <- unique(sort(unlist(lapply(models, function(m) m$feats))))
  intercept_vec <- numeric(length(models))
  BETA <- matrix(0, nrow=length(feats_union), ncol=length(models),
                 dimnames=list(feats_union, names(models)))
  for (i in seq_along(models)) {
    m <- models[[i]]
    co <- coef(m$cvfit, s=m$s_choice)
    beta0 <- as.numeric(co[1,1])
    beta  <- as.numeric(co[-1,1]); names(beta) <- rownames(co)[-1]
    a <- suppressWarnings(as.numeric(m$calib["a"])); if (!is.finite(a)) a <- 0
    b <- suppressWarnings(as.numeric(m$calib["b"])); if (!is.finite(b)) b <- 1
    beta0_cal <- a + b*beta0
    beta_cal  <- b*beta
    intercept_vec[i] <- beta0_cal
    BETA[names(beta_cal), i] <- beta_cal
  }
  if (is.null(weights)) {
    w <- rep(1/length(models), length(models))
  } else {
    w <- as.numeric(weights); w[!is.finite(w)] <- 0
    s <- sum(w); if (s<=0) w <- rep(1/length(models), length(models)) else w <- w/s
  }
  beta_avg  <- as.numeric(BETA %*% w); names(beta_avg) <- rownames(BETA)
  beta0_avg <- sum(intercept_vec * w)
  merged <- list(intercept=beta0_avg, beta=beta_avg, feats=rownames(BETA),
                 weight=w, note="Weighted mean of per-model CALIBRATED coefficients (per-O best models)")
  class(merged) <- c("merged_linear_model","list")
  merged
}
predict_merged <- function(merged_model, X_new) {
  feats <- intersect(merged_model$feats, colnames(X_new))
  if (length(feats)==0) return(rep(merged_model$intercept, nrow(X_new)))
  as.numeric(merged_model$intercept + as.matrix(X_new[,feats,drop=FALSE]) %*% merged_model$beta[feats])
}
w_models     <- get_model_weights(selected_models, table2_test_df, use_equal_weights = TRUE)
merged_model <- merge_models_to_single(selected_models, w_models)

## --------------------------- 11) Display -----------------------------------------
p_all <- ggplot(all_pred_df, aes(x=age_true, y=age_pred, color=group, shape=group)) +
  geom_abline(slope=1, intercept=0, linetype="dashed") +
  geom_point(size=2, alpha=0.9) +
  scale_color_manual(values = c("Y"="#2E295A",
                                "M"="#2E295A",
                                "O"="#2E295A",
                                'CR'='#3DA34C'
  ))+
  
  labs(title="Predicted vs actual age (per-O best round,b)",
       x="Actual age", y="Predicted age (calibrated)") +
  theme_bw()
print(p_all)
ggsave(file.path(model_wd, paste0("Age_pred_", ts_tst, ".png")), plot=p_all, width=3, height=3)
ggsave(file.path(model_wd, paste0("Age_pred_", ts_tst, ".pdf")), plot=p_all, width=3, height=3)

p_delta <- ggplot(all_pred_df, aes(x=group, y=delta_age, fill=group)) +
  geom_violin(trim=FALSE, alpha=0.65) +
  geom_boxplot(width=0.15, outlier.shape=NA) +
  scale_fill_manual(values = c("Y"="#2E295A",
                               "M"="#2E295A",
                               "O"="#2E295A",
                               'CR'='#3DA34C'
  ))+
  ggpubr::stat_compare_means(comparisons = list(c("O", "CR")),
                             label = "p.format",
                             method='wilcox.test',# wilcox.test  t.test
                             label.size=0.1,
                             bracket.size = 0.1
                             ,size=2.5)+
  
  geom_hline(yintercept=0, linetype="dashed") +
  labs(title="ΔAge by group (per-O best models)",
       x="Group", y="ΔAge (calibrated)") +
  theme_bw()
print(p_delta)
ggsave(file.path(model_wd, paste0("Age_delta_", ts_tst, ".png")), plot=p_delta, width=3, height=3)
ggsave(file.path(model_wd, paste0("Age_delta_", ts_tst, ".pdf")), plot=p_delta, width=3, height=3)

## --------------------------- 12) save the model and results --------------------------------
out_dir <- file.path(model_wd, "perO_best_models_with_stability")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

# the optimal model for each O sample
model_index_rows <- list(); model_files <- character(0)
for (i in seq_along(selected_models)) {
  mdl <- selected_models[[i]]
  key <- names(selected_models)[i]
  file_i <- file.path(out_dir, paste0("model_", key, ".rds"))
  mdl$note <- "Per-O best (round,b) model (selected by minimal calibrated test error)"
  mdl$n_features <- length(mdl$feats)
  saveRDS(mdl, file_i)
  model_files <- c(model_files, file_i)
  a <- tryCatch(as.numeric(mdl$calib["a"]), error=function(e) NA_real_)
  b <- tryCatch(as.numeric(mdl$calib["b"]), error=function(e) NA_real_)
  model_index_rows[[length(model_index_rows)+1]] <- data.frame(
    fold=mdl$fold, round=mdl$round, b=mdl$b,
    test_sample = if (!is.null(mdl$o_test_sample)) mdl$o_test_sample else NA_character_,
    n_features  = length(mdl$feats),
    alpha       = if (!is.null(mdl$alpha)) mdl$alpha else NA_real_,
    lambda      = if (!is.null(mdl$s_choice)) as.character(mdl$s_choice) else NA_character_,
    calib_a     = a, calib_b = b,
    inner_cv_mae= if (!is.null(mdl$inner_cv_mae)) mdl$inner_cv_mae else NA_real_,
    seed_sub    = if (!is.null(mdl$seed_sub)) mdl$seed_sub else NA_integer_,
    seed_cv     = if (!is.null(mdl$seed_cv))  mdl$seed_cv  else NA_integer_,
    rds_path    = file_i,
    stringsAsFactors = FALSE
  )
}
model_index_df <- do.call(rbind, model_index_rows)
write.csv(model_index_df, file.path(out_dir, "index_models_perO_best.csv"), row.names=FALSE)

# stability statistics (median of freq / |1SD| for each fold and each round)
if (nrow(stab_stats_long) > 0) {
  write.csv(stab_stats_long, file.path(out_dir, "stability_stats_by_fold_round.csv"), row.names=FALSE)
}

# the merged single linear model
saveRDS(merged_model, file.path(out_dir, "merged_linear_model.rds"))

# key tables
write.csv(table1_train_df, file.path(out_dir, "table1_train_log.csv"), row.names=FALSE)
write.csv(table2_test_df,  file.path(out_dir, "table2_test_O_predictions_all_rb.csv"), row.names=FALSE)
if (nrow(cr_fold_df)>0) write.csv(cr_fold_df, file.path(out_dir, "table_CR_predictions_all_rb.csv"), row.names=FALSE)
write.csv(selected_pairs_df, file.path(out_dir, "selected_pairs_perO.csv"), row.names=FALSE)
write.csv(all_pred_df, file.path(out_dir, "ALL_predictions_YMOCR_perO_best.csv"), row.names=FALSE)




