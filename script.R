pacman::p_load(readxl, survey, survival, gtsummary, survminer, splines, mice, timeROC, 
               dplyr, ggplot2, ggsci, rms, tidycmprsk, tidyr, nricens, tableone, patchwork)

raw_df <- read_excel("data.xlsx", sheet = "RAW DATA")

# ==============================================================================
# BƯỚC 1: LỌC DATA CƠ BẢN VÀ CHUẨN BỊ ĐIỀN KHUYẾT (IMPUTATION)
# ==============================================================================
cat("\n--- BƯỚC 1: LỌC DATA VÀ ĐIỀN KHUYẾT ---\n")

df_step1 <- raw_df %>% 
  filter(RIDAGEYR >= 20) %>%
  filter(
    !is.na(LBXTR) & !is.na(LBDSGLSI) & !is.na(LBDHDD) & 
      !is.na(LBXSATSI) & !is.na(LBXSASSI) &               
      !is.na(BMXBMI) & !is.na(BMXWAIST) & !is.na(BMXWT) & !is.na(BMXHT) & 
      !is.na(WTSAF2YR) & WTSAF2YR > 0 # Lọc luôn trọng số 0 để tránh warning
  )

df_step1 <- df_step1 %>%
  mutate(cvd_hx = ifelse(MCQ160B == 1 | MCQ160C == 1 | MCQ160E == 1 | MCQ160F == 1, 1, 0))

vars_to_impute <- c("INDFMPIR", "DMDEDUC2", "DMDMARTL", "LBXWBCSI", "LBXSCR", 
                    "LBXSUA", "cvd_hx", "diabetes", "hyperten")

cat_vars <- c("RIAGENDR", "RIDRETH1", "DMDEDUC2", "DMDMARTL", "diabetes", "hyperten", "cvd_hx")
df_step1[cat_vars] <- lapply(df_step1[cat_vars], as.factor)

imputed_data <- mice(df_step1[, vars_to_impute], m = 1, method = 'pmm', seed = 123, printFlag = FALSE)

df_imputed <- df_step1
df_imputed[, vars_to_impute] <- complete(imputed_data, 1)

df_imputed$diabetes <- as.numeric(as.character(df_imputed$diabetes))
df_imputed$hyperten <- as.numeric(as.character(df_imputed$hyperten))

# ==============================================================================
# BƯỚC 2: TÍNH TOÁN CHỈ SỐ VÀ CHẨN ĐOÁN MASLD
# ==============================================================================
cat("\n--- BƯỚC 2: TÍNH TOÁN CHỈ SỐ VÀ CHẨN ĐOÁN MASLD ---\n")

df_final <- df_imputed %>%
  mutate(
    HSI = 8 * (LBXSATSI / LBXSASSI) + BMXBMI + ifelse(RIAGENDR == 2, 2, 0) + ifelse(diabetes == 1, 2, 0),
    steatosis = ifelse(HSI > 36, 1, 0),
    
    cm_risk = ifelse(
      BMXBMI >= 25 | (RIAGENDR == 1 & BMXWAIST >= 94) | (RIAGENDR == 2 & BMXWAIST >= 80) | 
        LBDSGLSI >= 5.55 | diabetes == 1 | hyperten == 1 | LBXTR >= 150 | 
        (RIAGENDR == 1 & LBDHDD < 40) | (RIAGENDR == 2 & LBDHDD < 50), 1, 0),
    
    is_masld = ifelse(steatosis == 1 & cm_risk == 1, 1, 0),
    
    TyHGB = (LBXTR / LBDHDD) + (0.7 * LBDSGLSI) + (0.1 * BMXBMI),
    ABSI = BMXWAIST / ((BMXBMI^(2/3)) * (BMXHT^(1/2))),
    WWI = BMXWAIST / sqrt(BMXWT),
    WHtR = BMXWAIST / BMXHT,
    
    # Giữ nguyên dạng SỐ (0/1) cho phân tích Survival
    status_acm = ifelse(mortstat == 1, 1, 0),
    # Cập nhật format mã CVM bọc lót tất cả các trường hợp để không rớt ca nào
    status_cvm = ifelse(mortstat == 1 & ucod_leading %in% c(1, 5, "1", "5", "01", "05", "001", "005"), 1, 0),
    time_months = permth_int
  ) %>%
  filter(is_masld == 1, (is.na(RIDEXPRG) | RIDEXPRG != 1), !is.na(mortstat), !is.na(time_months), eligstat == 1) %>%
  mutate(
    TyHGB_Q = factor(ntile(TyHGB, 4), labels = c("Q1", "Q2", "Q3", "Q4")),
    # Tạo riêng 2 biến factor để dùng hiển thị gtsummary
    status_acm_factor = factor(status_acm),
    status_cvm_factor = factor(status_cvm),
    # Biến cho rủi ro cạnh tranh (Competing Risk)
    event_cr = case_when(mortstat == 0 ~ 0, status_cvm == 1 ~ 1, mortstat == 1 & status_cvm == 0 ~ 2)
  )

# Xử lý biến học vấn (Gộp nhóm)
df_final <- df_final %>%
  mutate(
    DMDEDUC2 = case_when(
      DMDEDUC2 %in% c(1, 2) ~ "1. Duoi cap 3",
      DMDEDUC2 == 3 ~ "2. Cap 3",
      DMDEDUC2 %in% c(4, 5) ~ "3. Dai hoc tro len",
      TRUE ~ NA_character_
    )
  ) %>%
  drop_na(DMDEDUC2)

## Xử lý biến hôn nhân & chủng tộc
df_final <- df_final %>%
  mutate(
    # Hôn nhân: Gộp thành 3 nhóm lớn
    DMDMARTL = case_when(
      DMDMARTL %in% c(1, 6) ~ "1. Married/Partner", # Đã kết hôn / Sống chung
      DMDMARTL %in% c(2, 3, 4) ~ "2. Widowed/Divorced/Separated", # Góa / Ly dị / Ly thân
      DMDMARTL == 5 ~ "3. Never married", # Chưa từng kết hôn
      TRUE ~ NA_character_
    ),
    # Chủng tộc: Gộp thành 3 nhóm lớn
    RIDRETH1 = case_when(
      RIDRETH1 == 3 ~ "1. Non-Hispanic White",
      RIDRETH1 == 4 ~ "2. Non-Hispanic Black",
      TRUE ~ "3. Other Races"
    )
  ) %>%
  drop_na(DMDMARTL, RIDRETH1)

# 3. Khai báo lại factor cho covariates TRƯỚC KHI tạo design mẫu
cat_vars_final <- c("RIAGENDR", "RIDRETH1", "DMDEDUC2", "DMDMARTL", "diabetes", "hyperten", "cvd_hx")
df_final <- df_final %>% mutate(across(all_of(cat_vars_final), as.factor))

# 4. TẠO DESIGN MẪU (Lúc này data đã sạch bong, không còn 1 mống NA)
masld_design <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTSAF2YR, nest = TRUE, data = df_final)

# ==============================================================================
# BƯỚC 3: THỐNG KÊ MÔ TẢ (TABLE 1)
# ==============================================================================
cat("\n--- BƯỚC 3: XUẤT TABLE 1 BẰNG GTSUMMARY ---\n")

vars_table1 <- c("RIDAGEYR", "RIAGENDR", "RIDRETH1", "DMDEDUC2", "DMDMARTL", "INDFMPIR",
                 "BMXBMI", "BMXWAIST", "TyHGB", "ABSI", "WWI", "WHtR", 
                 "LBXTC", "LBXTR", "LBDHDD", "LBDSGLSI", "LBXSATSI", "LBXSASSI",
                 "diabetes", "hyperten", "cvd_hx", "status_acm_factor", "status_cvm_factor")

tab1_gt <- masld_design %>%
  tbl_svysummary(
    by = TyHGB_Q, 
    include = all_of(vars_table1),
    statistic = list(all_continuous() ~ "{mean} ({sd})", all_categorical() ~ "{n_unweighted} ({p}%)"),
    digits = list(all_continuous() ~ 2, all_categorical() ~ c(0, 1)), missing = "no"
  ) %>%
  add_p() %>% add_overall() %>%
  modify_header(stat_0 = "**Overall**\n(N = {N_unweighted})", all_stat_cols() ~ "**{level}**\n(N = {n_unweighted})") %>%
  bold_labels() 

tab1_gt

# ==============================================================================
# BƯỚC 4: CHẠY MÔ HÌNH COX VÀ RÚT ĐIỂM DỰ BÁO (LINEAR PREDICTOR)
# ==============================================================================
cat("\n--- BƯỚC 4: RÚT ĐIỂM DỰ BÁO ---\n")

# Định nghĩa các covariates dạng chuỗi (KHÔNG CẦN dùng update)
covs_str <- "RIDAGEYR + RIAGENDR + RIDRETH1 + INDFMPIR + DMDEDUC2 + DMDMARTL + diabetes + hyperten + cvd_hx"

# CHẠY 4 MÔ HÌNH ACM
cox_acm_m1 <- svycoxph(as.formula(paste("Surv(time_months, status_acm) ~ TyHGB +", covs_str)), design = masld_design)
cox_acm_m2 <- svycoxph(as.formula(paste("Surv(time_months, status_acm) ~ TyHGB + ABSI +", covs_str)), design = masld_design)
cox_acm_m3 <- svycoxph(as.formula(paste("Surv(time_months, status_acm) ~ TyHGB + WWI +", covs_str)), design = masld_design)
cox_acm_m4 <- svycoxph(as.formula(paste("Surv(time_months, status_acm) ~ TyHGB + WHtR +", covs_str)), design = masld_design)

# CHẠY 4 MÔ HÌNH CVM
cox_cvm_m1 <- svycoxph(as.formula(paste("Surv(time_months, status_cvm) ~ TyHGB +", covs_str)), design = masld_design)
cox_cvm_m2 <- svycoxph(as.formula(paste("Surv(time_months, status_cvm) ~ TyHGB + ABSI +", covs_str)), design = masld_design)
cox_cvm_m3 <- svycoxph(as.formula(paste("Surv(time_months, status_cvm) ~ TyHGB + WWI +", covs_str)), design = masld_design)
cox_cvm_m4 <- svycoxph(as.formula(paste("Surv(time_months, status_cvm) ~ TyHGB + WHtR +", covs_str)), design = masld_design)

# Rút điểm dự báo
df_final <- df_final %>%
  mutate(
    lp_acm_m1 = predict(cox_acm_m1, type="lp"), Q_acm_m1 = factor(ntile(lp_acm_m1, 4), labels=paste0("Q", 1:4)),
    lp_acm_m2 = predict(cox_acm_m2, type="lp"), Q_acm_m2 = factor(ntile(lp_acm_m2, 4), labels=paste0("Q", 1:4)),
    lp_acm_m3 = predict(cox_acm_m3, type="lp"), Q_acm_m3 = factor(ntile(lp_acm_m3, 4), labels=paste0("Q", 1:4)),
    lp_acm_m4 = predict(cox_acm_m4, type="lp"), Q_acm_m4 = factor(ntile(lp_acm_m4, 4), labels=paste0("Q", 1:4)),
    
    lp_cvm_m1 = predict(cox_cvm_m1, type="lp"), Q_cvm_m1 = factor(ntile(lp_cvm_m1, 4), labels=paste0("Q", 1:4)),
    lp_cvm_m2 = predict(cox_cvm_m2, type="lp"), Q_cvm_m2 = factor(ntile(lp_cvm_m2, 4), labels=paste0("Q", 1:4)),
    lp_cvm_m3 = predict(cox_cvm_m3, type="lp"), Q_cvm_m3 = factor(ntile(lp_cvm_m3, 4), labels=paste0("Q", 1:4)),
    lp_cvm_m4 = predict(cox_cvm_m4, type="lp"), Q_cvm_m4 = factor(ntile(lp_cvm_m4, 4), labels=paste0("Q", 1:4))
  )

masld_design <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTSAF2YR, nest = TRUE, data = df_final)

# ==============================================================================
# BƯỚC 5: VẼ KAPLAN-MEIER GRID (2 HÀNG x 4 CỘT)
# ==============================================================================
cat("\n--- BƯỚC 5: VẼ KAPLAN-MEIER GHÉP KHỐI CHO ACM & CVM ---\n")

create_km <- function(group_var, event_var, model_name) {
  # 1. Chạy hàm lấy P-value từ thiết kế mẫu
  form_svy <- as.formula(paste("Surv(time_months,", event_var, ") ~", group_var))
  res_logrank <- svylogrank(form_svy, design = masld_design)
  
  # 2. Bóc tách P-value bọc thép (Xử lý việc survey package giấu P khi có >2 nhóm)
  if (!is.null(res_logrank$p.value)) {
    p_val <- res_logrank$p.value
  } else if (is.list(res_logrank) && !is.null(res_logrank[[2]]$p)) {
    p_val <- res_logrank[[2]]$p
  } else if (is.list(res_logrank) && !is.null(res_logrank[[1]]$p.value)) {
    p_val <- res_logrank[[1]]$p.value
  } else {
    p_val <- NA
  }
  
  p_val <- as.numeric(p_val)
  
  # 3. Format dòng chữ P-value in lên biểu đồ
  if (length(p_val) == 0 || is.na(p_val)) {
    p_str <- "Weighted P: N/A"
  } else if (p_val < 0.001) {
    p_str <- "Weighted P < 0.001"
  } else {
    p_str <- paste0("Weighted P = ", round(p_val, 3))
  }
  
  # 4. Fit mô hình unweighted để lấy số cho bảng Risk Table
  fit <- eval(parse(text = paste0("survfit(Surv(time_months, ", event_var, ") ~ ", group_var, ", data = df_final)")))
  
  # 5. Bảng màu chuẩn theo yêu cầu sếp
  my_palette <- c("#4daf4a", "#d95f02", "#7570b3", "#984ea3") 
  
  # 6. Vẽ biểu đồ với các tinh chỉnh thẩm mỹ
  p <- ggsurvplot(
    fit, 
    data = df_final, 
    palette = my_palette, 
    linewidth = 1,
    title = "",                               # Ẩn tên biểu đồ trên cùng
    legend.title = model_name,                # Tên mô hình (Sẽ được bôi đậm ở dưới)
    legend.labs = c("Q1", "Q2", "Q3", "Q4"),  # Chỉ hiện Q1-Q4 cho sạch sẽ
    xlab = "Months", 
    ylab = "Survival Probability",
    pval = p_str,                             # Bắn dòng P-value đã format lên hình
    pval.size = 4,                            
    conf.int = TRUE,                          
    conf.int.style = "ribbon",                
    conf.int.alpha = 0.15,                    
    censor = FALSE,                           # Tắt sạch các dấu gạch chéo
    risk.table = TRUE, 
    risk.table.height = 0.28, 
    fontsize = 3.5, 
    ggtheme = theme_classic()
  )
  
  # 7. Tinh chỉnh font chữ (In đậm tên mô hình, bỏ chữ trục Y của bảng Risk)
  p$plot <- p$plot + theme(
    legend.title = element_text(face = "bold", size = 11) 
  )
  p$table <- p$table + theme(axis.title.y = element_blank())
  
  return(p)
}

p_acm1 <- create_km("Q_acm_m1", "status_acm", "ACM: TyHGB")
p_acm2 <- create_km("Q_acm_m2", "status_acm", "ACM: TyHGB+ABSI")
p_acm3 <- create_km("Q_acm_m3", "status_acm", "ACM: TyHGB+WWI")
p_acm4 <- create_km("Q_acm_m4", "status_acm", "ACM: TyHGB+WHtR")

p_cvm1 <- create_km("Q_cvm_m1", "status_cvm", "CVM: TyHGB")
p_cvm2 <- create_km("Q_cvm_m2", "status_cvm", "CVM: TyHGB+ABSI")
p_cvm3 <- create_km("Q_cvm_m3", "status_cvm", "CVM: TyHGB+WWI")
p_cvm4 <- create_km("Q_cvm_m4", "status_cvm", "CVM: TyHGB+WHtR")

arrange_ggsurvplots(list(p_acm1, p_acm2, p_acm3, p_acm4, p_cvm1, p_cvm2, p_cvm3, p_cvm4), print = TRUE, ncol = 4, nrow = 2)

# ==============================================================================
# BƯỚC 6: VẼ ĐƯỜNG CONG ROC GHÉP KHỐI (2 HÀNG x 3 CỘT)
# ==============================================================================
cat("\n--- BƯỚC 6: VẼ ROC CURVES GHÉP KHỐI (ACM vs CVM tại 3 mốc TG) ---\n")

# Cho ACM
roc_a1 <- timeROC(T=df_final$time_months, delta=df_final$status_acm, marker=df_final$lp_acm_m1, cause=1, weighting="marginal", times=c(36,60,120))
roc_a2 <- timeROC(T=df_final$time_months, delta=df_final$status_acm, marker=df_final$lp_acm_m2, cause=1, weighting="marginal", times=c(36,60,120))
roc_a3 <- timeROC(T=df_final$time_months, delta=df_final$status_acm, marker=df_final$lp_acm_m3, cause=1, weighting="marginal", times=c(36,60,120))
roc_a4 <- timeROC(T=df_final$time_months, delta=df_final$status_acm, marker=df_final$lp_acm_m4, cause=1, weighting="marginal", times=c(36,60,120))

# Cho CVM (Xử lý Competing Risk với event_cr)
roc_c1 <- timeROC(T=df_final$time_months, delta=df_final$event_cr, marker=df_final$lp_cvm_m1, cause=1, weighting="marginal", times=c(36,60,120))
roc_c2 <- timeROC(T=df_final$time_months, delta=df_final$event_cr, marker=df_final$lp_cvm_m2, cause=1, weighting="marginal", times=c(36,60,120))
roc_c3 <- timeROC(T=df_final$time_months, delta=df_final$event_cr, marker=df_final$lp_cvm_m3, cause=1, weighting="marginal", times=c(36,60,120))
roc_c4 <- timeROC(T=df_final$time_months, delta=df_final$event_cr, marker=df_final$lp_cvm_m4, cause=1, weighting="marginal", times=c(36,60,120))

plot_roc <- function(r1, r2, r3, r4, time_idx, title) {
  df_roc <- data.frame(
    FPR = c(r1$FP[,time_idx], r2$FP[,time_idx], r3$FP[,time_idx], r4$FP[,time_idx]),
    TPR = c(r1$TP[,time_idx], r2$TP[,time_idx], r3$TP[,time_idx], r4$TP[,time_idx]),
    Model = c(rep(paste0("TyHGB (AUC=", round(r1$AUC[time_idx], 3), ")"), length(r1$FP[,time_idx])),
              rep(paste0("TyHGB+ABSI (AUC=", round(r2$AUC[time_idx], 3), ")"), length(r2$FP[,time_idx])),
              rep(paste0("TyHGB+WWI (AUC=", round(r3$AUC[time_idx], 3), ")"), length(r3$FP[,time_idx])),
              rep(paste0("TyHGB+WHtR (AUC=", round(r4$AUC[time_idx], 3), ")"), length(r4$FP[,time_idx])))
  )
  ggplot(df_roc, aes(x = FPR, y = TPR, color = Model)) +
    geom_line(linewidth = 0.8) + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    scale_color_nejm() + theme_bw() +
    theme(legend.position = c(0.65, 0.2), legend.title = element_blank(),
          legend.background = element_rect(fill = alpha("white", 0.5), color = NA),
          legend.text = element_text(size = 8), title = element_text(size = 11, face = "bold")) +
    labs(title = title, x = "1 - Specificity", y = "Sensitivity")
}

p_roc_a1 <- plot_roc(roc_a1, roc_a2, roc_a3, roc_a4, 1, "ACM - 3 Years")
p_roc_a2 <- plot_roc(roc_a1, roc_a2, roc_a3, roc_a4, 2, "ACM - 5 Years")
p_roc_a3 <- plot_roc(roc_a1, roc_a2, roc_a3, roc_a4, 3, "ACM - 10 Years")

p_roc_c1 <- plot_roc(roc_c1, roc_c2, roc_c3, roc_c4, 1, "CVM - 3 Years")
p_roc_c2 <- plot_roc(roc_c1, roc_c2, roc_c3, roc_c4, 2, "CVM - 5 Years")
p_roc_c3 <- plot_roc(roc_c1, roc_c2, roc_c3, roc_c4, 3, "CVM - 10 Years")

final_roc_plot <- (p_roc_a1 | p_roc_a2 | p_roc_a3) / (p_roc_c1 | p_roc_c2 | p_roc_c3)
print(final_roc_plot)

# ==============================================================================
# BƯỚC 7, 8, 9 (GIỮ NGUYÊN NHƯ FILE TRƯỚC LÀ CHẠY ĐƯỢC 100%)
# ==============================================================================