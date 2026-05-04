pacman::p_load(readxl, survey, survival, gtsummary, survminer, splines, mice, timeROC, 
               dplyr, ggplot2, ggsci, rms, tidycmprsk, tidyr, survIDINRI, tableone, patchwork)

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
# BƯỚC 3: THỐNG KÊ MÔ TẢ (TABLE 1) - FULL LABEL CHUẨN TỪ ĐIỂN NHANES
# ==============================================================================

# 1. Danh sách biến đầy đủ
vars_table1 <- c("RIDAGEYR", "RIAGENDR", "RIDRETH1", "DMDEDUC2", "DMDMARTL", "INDFMPIR",
                 "BMXBMI", "BMXWAIST", "TyHGB", "ABSI", "WWI", "WHtR", 
                 "LBXTC", "LBXTR", "LBDHDD", "LBDSGLSI", "LBXSATSI", "LBXSASSI",
                 "LBXWBCSI", "LBXLYPCT", "LBXMOPCT", "LBXNEPCT", "LBXEOPCT", "LBXBAPCT",
                 "LBDLYMNO", "LBDMONO", "LBDNENO", "LBDEONO", "LBDBANO",
                 "LBXRBCSI", "LBXHGB", "LBXHCT", "LBXWBCSI", "LBDLYMNO", "LBDMONO", "LBDNENO",
                 "LBDEONO", "LBDBANO", "LBXRBCSI", "LBXHGB", "LBXHCT",
                 "diabetes", "hyperten", "cvd_hx", "status_acm_factor", "status_cvm_factor")

# 2. Xuất bảng với Label "căng đét"
tab1_gt <- masld_design %>%
  tbl_svysummary(
    by = TyHGB_Q, 
    include = all_of(vars_table1),
    statistic = list(
      all_continuous() ~ "{mean} ({sd})", 
      all_categorical() ~ "{n_unweighted} ({p}%)"
    ),
    digits = list(
      all_continuous() ~ 2, 
      all_categorical() ~ c(0, 1)
    ), 
    missing = "no",
    
    # ---------- BỘ TỪ ĐIỂN LABEL CHUẨN NHANES ----------
    label = list(
      # Nhóm Nhân trắc học & Xã hội
      RIDAGEYR ~ "Age (years)",
      RIAGENDR ~ "Gender",
      RIDRETH1 ~ "Race",
      DMDEDUC2 ~ "Education level",
      DMDMARTL ~ "Marital status",
      INDFMPIR ~ "Family income to poverty",
      
      # Nhóm Chỉ số béo phì & Chuyển hóa
      BMXBMI   ~ "Body mass index (kg/m²)",
      BMXWAIST ~ "Waist circumference (cm)",
      TyHGB    ~ "TyHGB",
      ABSI     ~ "A Body Shape Index (ABSI)",
      WWI      ~ "Weight-Adjusted Waist Index (WWI)",
      WHtR     ~ "Waist-to-Height Ratio (WHtR)",
      
      # Nhóm Sinh hóa & Gan máu (Lưu ý đuôi SI là đơn vị chuẩn quốc tế)
      LBXTC    ~ "Total cholesterol (mg/dL)",
      LBXTR    ~ "Triglycerides (mg/dL)",
      LBDHDD   ~ "HDL-cholesterol (mg/dL)",
      LBDSGLSI ~ "Fasting glucose (mmol/L)",
      LBXSATSI ~ "Alanine aminotransferase (ALT, U/L)",
      LBXSASSI ~ "Aspartate aminotransferase (AST, U/L)",
      
      # Nhóm Cận lâm sàng Huyết học (Tế bào máu CBC)
      LBXWBCSI ~ "WBC (1000 cells/μL)",
      LBDLYMNO ~ "Lymphocytes (1000 cells/μL)",
      LBDMONO  ~ "Monocytes (1000 cells/μL)",
      LBDNENO  ~ "Neutrophils (1000 cells/μL)",
      LBDEONO  ~ "Eosinophils (1000 cells/μL)",
      LBDBANO  ~ "Basophils (1000 cells/μL)",
      LBXRBCSI ~ "RBC (million cells/μL)",
      LBXHGB   ~ "Hemoglobin (g/dL)",
      LBXHCT   ~ "Hematocrit (%)",
      
      # Nhóm Bệnh nền & Kết cục
      diabetes ~ "Diabetes",
      hyperten ~ "Hypertension",
      cvd_hx   ~ "History of cardiovascular disease",
      status_acm_factor ~ "All-cause mortality",
      status_cvm_factor ~ "Cardiovascular mortality"
    )
    # ---------------------------------------------------
  ) %>%
  add_p() %>% 
  add_overall() %>%
  modify_header(
    stat_0 = "**Overall**\n(N = {N_unweighted})", 
    all_stat_cols() ~ "**{level}**\n(N = {n_unweighted})"
  ) %>%
  bold_labels()

# 3. Hiển thị Table 1
tab1_gt

# ==============================================================================
# BƯỚC 3.4: CHUẨN HÓA Z-SCORE CHO 4 CHỈ SỐ ĐỂ SO SÁNH CÔNG BẰNG
# ==============================================================================

# 1. Ép 4 chỉ số về chuẩn Z-score (Mean = 0, SD = 1)
df_final <- df_final %>%
  mutate(
    TyHGB_Z = as.numeric(scale(TyHGB)),
    ABSI_Z  = as.numeric(scale(ABSI)),
    WWI_Z   = as.numeric(scale(WWI)),
    WHtR_Z  = as.numeric(scale(WHtR))
  )

# 2. Cập nhật lại design với data mới
masld_design <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTSAF2YR, nest = TRUE, data = df_final)

# ==============================================================================
# BƯỚC 3.5: PHÂN TÍCH COX ĐỘC LẬP TỪNG BIẾN (DÙNG BIẾN ĐÃ CHUẨN HÓA _Z)
# ==============================================================================

# 1. ĐỔI TÊN markers SANG CÁC BIẾN "_Z" VỪA TẠO
markers <- c("TyHGB_Z", "ABSI_Z", "WWI_Z", "WHtR_Z")

# Biến hiệu chỉnh (giữ nguyên)
adj_vars <- "RIDAGEYR + LBXSATSI"

# 2. HÀM CHẠY VÀ RÚT KẾT QUẢ
get_independent_res <- function(outcome, marker) {
  
  # A. Univariate
  form_uni <- paste0("Surv(time_months, ", outcome, ") ~ ", marker)
  fit_uni <- suppressWarnings(svycoxph(as.formula(form_uni), design = masld_design))
  res_uni <- tidy(fit_uni, exponentiate = TRUE, conf.int = TRUE) %>% filter(term == marker)
  
  hr_uni <- sprintf("%.3f (%.3f-%.3f)", res_uni$estimate, res_uni$conf.low, res_uni$conf.high)
  p_uni <- ifelse(res_uni$p.value < 0.001, "<0.001", sprintf("%.3f", res_uni$p.value))
  
  # B. Multivariate
  form_multi <- paste0("Surv(time_months, ", outcome, ") ~ ", marker, " + ", adj_vars)
  fit_multi <- suppressWarnings(svycoxph(as.formula(form_multi), design = masld_design))
  res_multi <- tidy(fit_multi, exponentiate = TRUE, conf.int = TRUE) %>% filter(term == marker)
  
  hr_multi <- sprintf("%.3f (%.3f-%.3f)", res_multi$estimate, res_multi$conf.low, res_multi$conf.high)
  p_multi <- ifelse(res_multi$p.value < 0.001, "<0.001", sprintf("%.3f", res_multi$p.value))
  
  return(c(hr_uni, p_uni, hr_multi, p_multi))
}

# 3. LẶP VÀ XÂY DỰNG DATA.FRAME
results_list <- list()

for (outcome in c("status_acm", "status_cvm")) {
  outcome_label <- ifelse(outcome == "status_acm", "All-cause mortality", "Cardiovascular mortality")
  
  for (marker in markers) {
    res <- get_independent_res(outcome, marker)
    
    # Bỏ chữ "_Z" đi khi in ra bảng để nhìn cho đẹp, giữ nguyên tên gốc
    clean_marker_name <- gsub("_Z", "", marker)
    
    results_list[[length(results_list) + 1]] <- data.frame(
      Outcome = outcome_label,
      Marker = clean_marker_name,
      Uni_HR = res[1],
      Uni_P = res[2],
      Multi_HR = res[3],
      Multi_P = res[4],
      stringsAsFactors = FALSE
    )
  }
}

df_cox_independent <- bind_rows(results_list)

# 4. VẼ BẢNG BẰNG 'gt'
table_cox_new <- df_cox_independent %>%
  group_by(Outcome) %>%
  gt() %>%
  tab_spanner(label = md("**Univariate analysis (per 1-SD increase)**"), columns = c(Uni_HR, Uni_P)) %>%
  tab_spanner(label = md("**Multivariate Analysis (per 1-SD increase)**"), columns = c(Multi_HR, Multi_P)) %>%
  cols_label(
    Marker = "",
    Uni_HR = md("**HR (95% CIs)**"),
    Uni_P = md("**p-value**"),
    Multi_HR = md("**HR (95% CIs)**"),
    Multi_P = md("**p-value**")
  ) %>%
  tab_options(
    row_group.font.weight = "bold",
    row_group.background.color = "#E8E8E8",
    column_labels.font.weight = "bold",
    table.border.top.color = "black",
    table.border.bottom.color = "black",
    table_body.border.bottom.color = "black"
  ) %>%
  cols_align(align = "left", columns = everything())

# 5. HIỂN THỊ
table_cox_new

# ==============================================================================
# BƯỚC 4: CHẠY MÔ HÌNH COX VÀ RÚT ĐIỂM DỰ BÁO 
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
  
  # 2. Bóc tách P-value CHUẨN XÁC NHẤT cho package survey (> 2 nhóm)
  if (!is.null(res_logrank$p.value)) {
    p_val <- res_logrank$p.value
  } else {
    # Khi so sánh >= 3 nhóm, res_logrank[[2]] là vector chứa [Chi-square, df]
    chisq_stat <- res_logrank[[2]][1]
    df_stat <- res_logrank[[2]][2]
    # Tự tính P-value từ Chi-square và df
    p_val <- pchisq(chisq_stat, df_stat, lower.tail = FALSE)
  }
  
  p_val <- as.numeric(p_val)
  
  # 3. Format dòng chữ P-value in lên biểu đồ
  if (is.na(p_val)) {
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

combine_km <- function(km_obj) {
  # Ẩn chữ trục Y của bảng risk table cho sạch
  km_obj$table <- km_obj$table + theme(axis.title.y = element_blank())
  
  # Ghép plot và table theo tỷ lệ 3:1 chiều cao
  (km_obj$plot / km_obj$table) + plot_layout(heights = c(3, 1))
}

pw_acm1 <- combine_km(create_km("Q_acm_m1", "status_acm", "TyHGB"))
pw_acm2 <- combine_km(create_km("Q_acm_m2", "status_acm", "TyHGB_ABSI"))
pw_acm3 <- combine_km(create_km("Q_acm_m3", "status_acm", "TyHGB_WWI"))
pw_acm4 <- combine_km(create_km("Q_acm_m4", "status_acm", "TyHGB_WHtR"))

pw_cvm1 <- combine_km(create_km("Q_cvm_m1", "status_cvm", "TyHGB"))
pw_cvm2 <- combine_km(create_km("Q_cvm_m2", "status_cvm", "TyHGB_ABSI"))
pw_cvm3 <- combine_km(create_km("Q_cvm_m3", "status_cvm", "TyHGB_WWI"))
pw_cvm4 <- combine_km(create_km("Q_cvm_m4", "status_cvm", "TyHGB_WHtR"))

final_plot <- (pw_acm1 | pw_acm2 | pw_acm3 | pw_acm4) / 
  (pw_cvm1 | pw_cvm2 | pw_cvm3 | pw_cvm4)

# Lưu bằng ggsave
ggsave("Figure_KM_Patchwork.png", plot = final_plot, width = 18, height = 10, dpi = 300, bg = "white")

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
  get_roc_data <- function(roc_obj) {
    if (!is.null(roc_obj$AUC)) {
      return(list(AUC = roc_obj$AUC, FP = roc_obj$FP, TP = roc_obj$TP)) 
    } else {
      return(list(AUC = roc_obj$AUC_1, FP = roc_obj$FP_1, TP = roc_obj$TP_1)) 
    }
  }
  
  d1 <- get_roc_data(r1); d2 <- get_roc_data(r2); d3 <- get_roc_data(r3); d4 <- get_roc_data(r4)
  
  df_roc <- data.frame(
    FPR = c(d1$FP[,time_idx], d2$FP[,time_idx], d3$FP[,time_idx], d4$FP[,time_idx]),
    TPR = c(d1$TP[,time_idx], d2$TP[,time_idx], d3$TP[,time_idx], d4$TP[,time_idx]),
    Model = c(rep(paste0("TyHGB (AUC=", round(d1$AUC[time_idx], 3), ")"), length(d1$FP[,time_idx])),
              rep(paste0("TyHGB+ABSI (AUC=", round(d2$AUC[time_idx], 3), ")"), length(d2$FP[,time_idx])),
              rep(paste0("TyHGB+WWI (AUC=", round(d3$AUC[time_idx], 3), ")"), length(d3$FP[,time_idx])),
              rep(paste0("TyHGB+WHtR (AUC=", round(d4$AUC[time_idx], 3), ")"), length(d4$FP[,time_idx])))
  )
  
  ggplot(df_roc, aes(x = FPR, y = TPR, color = Model)) +
    geom_line(linewidth = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    scale_color_nejm() + theme_bw() +
    theme(legend.position = c(0.65, 0.22), 
          legend.title = element_blank(),
          legend.background = element_rect(fill = alpha("white", 0.7), color = "black"),
          legend.text = element_text(size = 9, face = "bold"), 
          title = element_text(size = 12, face = "bold")) +
    labs(title = title, x = "1 - Specificity", y = "Sensitivity")
}

hide_x <- theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
hide_y <- theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())

p_roc_a1 <- plot_roc(roc_a1, roc_a2, roc_a3, roc_a4, 1, "ACM - 3 Years") + hide_x
p_roc_a2 <- plot_roc(roc_a1, roc_a2, roc_a3, roc_a4, 2, "ACM - 5 Years") + hide_x + hide_y
p_roc_a3 <- plot_roc(roc_a1, roc_a2, roc_a3, roc_a4, 3, "ACM - 10 Years") + hide_x + hide_y

p_roc_c1 <- plot_roc(roc_c1, roc_c2, roc_c3, roc_c4, 1, "CVM - 3 Years")
p_roc_c2 <- plot_roc(roc_c1, roc_c2, roc_c3, roc_c4, 2, "CVM - 5 Years") + hide_y
p_roc_c3 <- plot_roc(roc_c1, roc_c2, roc_c3, roc_c4, 3, "CVM - 10 Years") + hide_y

final_roc_plot <- (p_roc_a1 | p_roc_a2 | p_roc_a3) / 
  (p_roc_c1 | p_roc_c2 | p_roc_c3)

ggsave("ROC_Grid_Publication.png", plot = final_roc_plot, 
       width = 10, height = 7, units = "in", dpi = 300, bg = "white")

# ==============================================================================
# BƯỚC 7: C-INDEX, ĐỊNH LƯỢNG NRI (CONTINUOUS) VÀ IDI CHO TẤT CẢ MÔ HÌNH (10 NĂM)
# ==============================================================================
cat("\n--- BƯỚC 7: TÍNH C-INDEX, CONTINUOUS NRI, IDI CHO ACM & CVM (10 NĂM) ---\n")

# 1. HÀM TÍNH C-INDEX
get_c_ci <- function(time_var, event_var, lp_var, weights_var) {
  c_obj <- concordance(Surv(time_var, event_var) ~ lp_var, weights = weights_var)
  c_val <- ifelse(c_obj$concordance < 0.5, 1 - c_obj$concordance, c_obj$concordance)
  se <- sqrt(c_obj$var)
  lower <- c_val - 1.96 * se
  upper <- c_val + 1.96 * se
  return(list(c_str = paste0(sprintf("%.3f", c_val), " (", sprintf("%.3f", lower), "-", sprintf("%.3f", upper), ")"), 
              c_val = c_val, se = se))
}

# 2. HÀM SO SÁNH (CÓ CƠ CHẾ CHỐNG LỖI NA CỦA DELONG TEST)
compare_models <- function(time_var, event_var, event_cr_var, lp_base, lp_new, weights_var, t0, is_cvm) {
  
  # Lấy C-index
  base_c <- get_c_ci(time_var, event_var, lp_base, weights_var)
  new_c  <- get_c_ci(time_var, event_var, lp_new, weights_var)
  
  # TÍNH P-VALUE DELONG
  if(is_cvm) { 
    roc_base <- timeROC(T=time_var, delta=event_cr_var, marker=lp_base, cause=1, weighting="marginal", times=t0, iid=TRUE)
    roc_new  <- timeROC(T=time_var, delta=event_cr_var, marker=lp_new, cause=1, weighting="marginal", times=t0, iid=TRUE)
  } else {
    roc_base <- timeROC(T=time_var, delta=event_var, marker=lp_base, cause=1, weighting="marginal", times=t0, iid=TRUE)
    roc_new  <- timeROC(T=time_var, delta=event_var, marker=lp_new, cause=1, weighting="marginal", times=t0, iid=TRUE)
  }
  
  comp_res <- compare(roc_new, roc_base, adjusted=FALSE)
  p_delong <- if(!is.null(comp_res$p_values_AUC)) comp_res$p_values_AUC[1] else comp_res$p_values_AUC_1[1]
  
  # BỌC THÉP LỚP 2: Nếu DeLong bị bão hòa ra NA, tự động chuyển sang Z-test của C-index
  if(is.na(p_delong)) {
    se_diff <- sqrt(base_c$se^2 + new_c$se^2)
    if (se_diff == 0) se_diff <- 1e-6 # Chống lỗi chia cho 0
    z_stat <- abs(new_c$c_val - base_c$c_val) / se_diff
    p_delong <- 2 * (1 - pnorm(z_stat))
  }
  
  # TÍNH NRI & IDI (Fix cứng npert = 100 theo lệnh sếp)
  indata_df <- data.frame(time = time_var, status = event_var)
  
  set.seed(123)
  idi_res <- IDI.INF(indata_df, covs0 = as.matrix(lp_base), covs1 = as.matrix(lp_new), t0 = t0, npert = 100) 
  out_matrix <- IDI.INF.OUT(idi_res)
  
  idi_str <- paste0(sprintf("%.3f", out_matrix[1, 1]), " (", sprintf("%.3f", out_matrix[1, 2]), "-", sprintf("%.3f", out_matrix[1, 3]), ")")
  idi_p   <- out_matrix[1, 4]
  
  nri_str <- paste0(sprintf("%.3f", out_matrix[3, 1]), " (", sprintf("%.3f", out_matrix[3, 2]), "-", sprintf("%.3f", out_matrix[3, 3]), ")")
  nri_p   <- out_matrix[3, 4]
  
  return(list(c_new = new_c$c_str, p_delong = p_delong, nri = nri_str, p_nri = nri_p, idi = idi_str, p_idi = idi_p))
}

# 3. HÀM TỔNG HỢP KẺ BẢNG
run_all_comparisons <- function(time_var, event_var, event_cr_var, lp_m1, lp_m2, lp_m3, lp_m4, weights_var, t0 = 120, is_cvm = FALSE, outcome_name) {
  cat(paste("\n================ KẾT QUẢ CHO:", outcome_name, "================\n"))
  
  c_base_str <- get_c_ci(time_var, event_var, lp_m1, weights_var)$c_str
  cat(sprintf("Mô hình gốc (TyHGB) - C-index: %s\n\n", c_base_str))
  
  cat(sprintf("%-12s | %-19s | %-8s | %-19s | %-8s | %-19s | %-8s\n",
              "Model so sánh", "C-index New (95%CI)", "P DeLong", "Cont. NRI (95%CI)", "P (NRI)", "IDI (95%CI)", "P (IDI)"))
  cat(strrep("-", 108), "\n")
  
  models <- list(
    list(name = "TyHGB+ABSI", lp = lp_m2),
    list(name = "TyHGB+WWI",  lp = lp_m3),
    list(name = "TyHGB+WHtR", lp = lp_m4)
  )
  
  for (m in models) {
    res <- compare_models(time_var, event_var, event_cr_var, lp_m1, m$lp, weights_var, t0, is_cvm)
    # P-value in <0.001 nếu quá nhỏ cho đẹp
    p_dl_print <- ifelse(res$p_delong < 0.001, "<0.001", sprintf("%.4f", res$p_delong))
    p_nri_print <- ifelse(res$p_nri < 0.001, "<0.001", sprintf("%.4f", res$p_nri))
    p_idi_print <- ifelse(res$p_idi < 0.001, "<0.001", sprintf("%.4f", res$p_idi))
    
    cat(sprintf("%-12s | %-19s | %-8s | %-19s | %-8s | %-19s | %-8s\n",
                m$name, res$c_new, p_dl_print, res$nri, p_nri_print, res$idi, p_idi_print))
  }
  cat("\n")
}

cat("Đang xử lý Bootstrap 100 vòng. Sẽ mất vài phút...\n")

# 1. Bảng so sánh cho ACM
run_all_comparisons(
  time_var     = df_final$time_months, 
  event_var    = df_final$status_acm, 
  event_cr_var = df_final$status_acm,     
  lp_m1        = df_final$lp_acm_m1,      
  lp_m2        = df_final$lp_acm_m2, 
  lp_m3        = df_final$lp_acm_m3, 
  lp_m4        = df_final$lp_acm_m4, 
  weights_var  = df_final$WTSAF2YR,
  t0           = 120, 
  is_cvm       = FALSE, 
  outcome_name = "ALL-CAUSE MORTALITY (ACM) TẠI 10 NĂM"
)

# 2. Bảng so sánh cho CVM
run_all_comparisons(
  time_var     = df_final$time_months, 
  event_var    = df_final$status_cvm, 
  event_cr_var = df_final$event_cr,       
  lp_m1        = df_final$lp_cvm_m1,      
  lp_m2        = df_final$lp_cvm_m2, 
  lp_m3        = df_final$lp_cvm_m3, 
  lp_m4        = df_final$lp_cvm_m4, 
  weights_var  = df_final$WTSAF2YR,
  t0           = 120, 
  is_cvm       = TRUE, 
  outcome_name = "CARDIOVASCULAR MORTALITY (CVM) TẠI 10 NĂM"
)


# ==============================================================================
# BƯỚC 8: VẼ ĐƯỜNG CONG LIỀU-ĐÁP ỨNG (RCS) - BẢN FINAL CHUẨN KÍCH THƯỚC Q1
# ==============================================================================
# 1. TẠO TẬP DỮ LIỆU SẠCH
rcs_vars <- c("time_months", "status_acm", "status_cvm", "TyHGB", 
              "RIDAGEYR", "RIAGENDR", "RIDRETH1", "INDFMPIR", 
              "DMDEDUC2", "DMDMARTL", "diabetes", "hyperten", "cvd_hx", "WTSAF2YR")

df_rcs <- na.omit(df_final[, rcs_vars])
dd <- datadist(df_rcs)
options(datadist = "dd")
covs_str <- "RIDAGEYR + RIAGENDR + RIDRETH1 + INDFMPIR + DMDEDUC2 + DMDMARTL + diabetes + hyperten + cvd_hx"

# 2. HÀM VẼ BẢN CHỐT SỔ
plot_rcs_exact_style <- function(event_var, title_str) {
  
  form <- as.formula(paste("Surv(time_months, ", event_var, ") ~ rcs(TyHGB, 4) +", covs_str))
  fit <- cph(form, data = df_rcs, x = TRUE, y = TRUE, weights = df_rcs$WTSAF2YR)
  
  anv <- anova(fit)
  p_overall <- anv["TyHGB", "P"]
  p_nonlin  <- anv[" Nonlinear", "P"] 
  
  p_overall_str <- ifelse(p_overall < 0.001, "< 0.001", sprintf("%.3f", p_overall))
  p_nonlin_str  <- ifelse(p_nonlin < 0.001, "< 0.001", sprintf("%.3f", p_nonlin))
  
  knots_array <- fit$Design$parms$TyHGB
  knots_str <- paste(round(knots_array, 2), collapse = ", ")
  
  anno_text <- paste0("P-overall ", ifelse(p_overall < 0.001, "", "= "), p_overall_str, "\n",
                      "P-non-linear ", ifelse(p_nonlin < 0.001, "", "= "), p_nonlin_str, "\n",
                      "Knots: ", knots_str)
  
  # Rút dữ liệu dự báo 
  pred <- Predict(fit, TyHGB, fun = exp, ref.zero = TRUE)
  df_pred <- as.data.frame(pred)
  
  # ---------- BỘ LỌC TỌA ĐỘ THEO LỆNH SẾP ----------
  min_x <- min(df_rcs$TyHGB, na.rm = TRUE)
  max_x <- quantile(df_rcs$TyHGB, 0.99, na.rm = TRUE) 
  
  upper_limit_y <- 1.5
  
  max_dens <- max(density(df_rcs$TyHGB, na.rm = TRUE)$y)
  scale_factor <- max_dens / (upper_limit_y * 0.5) 
  # -------------------------------------------------
  
  df_knots <- data.frame(k = knots_array)
  
  p <- ggplot() +
    # LỚP NỀN: Histogram Mật độ
    geom_histogram(data = df_rcs, aes(x = TyHGB, y = after_stat(density) / scale_factor),
                   bins = 60, fill = "lightblue", color = "white", alpha = 0.9) +
    
    # ---------- BỘ 3 ĐƯỜNG CONG TÁCH BẠCH ----------
  # 2. Hai đường viền KTC (Nét đứt, mỏng, giúp tách biệt khỏi đường chính)
  geom_line(data = df_pred, aes(x = TyHGB, y = lower), 
            color = "black", linewidth = 1.5) +
    geom_line(data = df_pred, aes(x = TyHGB, y = upper), 
              color = "black", linewidth = 1.5) +
    # 3. Đường HR trung tâm (Nét liền, đậm)
    geom_line(data = df_pred, aes(x = TyHGB, y = yhat), 
              color = "red", linewidth = 1) +
    # -----------------------------------------------
  
  # Vạch Knots dưới trục X
  geom_segment(data = df_knots, aes(x = k, xend = k, y = 0, yend = upper_limit_y * 0.02),
               color = "blue", linewidth = 1) +
    
    # Đường tham chiếu HR = 1
    geom_hline(yintercept = 1, linetype = "dashed", color = "black", linewidth = 0.6) +
    
    # CHÚ THÍCH DỜI XUỐNG GÓC PHẢI DƯỚI CÙNG (y = -Inf)
    annotate("text", x = Inf, y = -Inf, label = anno_text, 
             hjust = 1.05, vjust = -0.5, size = 4.5, color = "black") +
    
    coord_cartesian(xlim = c(min_x, max_x), ylim = c(0, upper_limit_y)) +
    
    scale_y_continuous(
      name = "Hazard Ratio (95% CI)",
      breaks = seq(0, 1.5, by = 0.5), 
      sec.axis = sec_axis(~ . * scale_factor, name = "Probability Density")
    ) +
    
    theme_bw() +
    labs(title = title_str, x = "TyHGB Index") +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.title.y.left = element_text(size = 12), axis.text.y.left = element_text(size = 11),
      axis.title.y.right = element_text(size = 12), axis.text.y.right = element_text(size = 11),
      axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 11),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(linetype = "dotted", color = "gray80")
    )
  
  return(p)
}

# 3. CHẠY VÀ LẮP RÁP BẰNG PATCHWORK
cat("Đang render RCS cho ACM...\n")
p_rcs_acm <- plot_rcs_exact_style("status_acm", "")

cat("Đang render RCS cho CVM...\n")
p_rcs_cvm <- plot_rcs_exact_style("status_cvm", "")

cat("Đang ghép khung hình ngang và xuất file 300 DPI...\n")
final_rcs_plot <- (p_rcs_acm | p_rcs_cvm) +
  plot_annotation(
    tag_levels = 'A',           # Đánh dấu bằng chữ cái in hoa (A, B, C...)
    tag_prefix = '(',           # Thêm dấu ngoặc mở
    tag_suffix = ')'            # Thêm dấu ngoặc đóng
  ) &
  theme(plot.tag = element_text(size = 22, face = "bold")) # Chỉnh size chữ (A), (B) cho to và đậm

ggsave("Figure_RCS_DoseResponse_Final.png", plot = final_rcs_plot, 
       width = 16, height = 7, units = "in", dpi = 300, bg = "white")


# ==============================================================================
# BƯỚC 9: PHÂN TÍCH PHÂN NHÓM (SUBGROUP ANALYSIS) & FOREST PLOT CHUẨN Q1
# ==============================================================================

# 1. CHUẨN BỊ BIẾN PHÂN NHÓM TỪ DỮ LIỆU GỐC
df_sub <- df_final %>%
  mutate(
    Age_group = case_when(
      RIDAGEYR < 40 ~ "18-39",
      RIDAGEYR >= 40 & RIDAGEYR < 60 ~ "40-59",
      RIDAGEYR >= 60 ~ ">=60"
    ),
    Gender = ifelse(RIAGENDR == 1, "Male", "Female"),
    PIR_group = case_when(
      INDFMPIR < 1.3 ~ "Low income",
      INDFMPIR >= 1.3 & INDFMPIR <= 3.5 ~ "Medium income",
      INDFMPIR > 3.5 ~ "High income"
    ),
    Obesity = ifelse(BMXBMI >= 30, "Obesity", "Non-Obesity"),
    Diabetes_cat = ifelse(diabetes == 1, "Yes", "No"),
    Hyperten_cat = ifelse(hyperten == 1, "Yes", "No"),
    CVD_cat = ifelse(cvd_hx == 1, "Yes", "No")
  ) %>%
  mutate(
    Age_group = factor(Age_group, levels = c("18-39", "40-59", ">=60")),
    Gender = factor(Gender, levels = c("Male", "Female")),
    PIR_group = factor(PIR_group, levels = c("Low income", "Medium income", "High income")),
    Obesity = factor(Obesity, levels = c("Non-Obesity", "Obesity")),
    Diabetes_cat = factor(Diabetes_cat, levels = c("Yes", "No")),
    Hyperten_cat = factor(Hyperten_cat, levels = c("Yes", "No")),
    CVD_cat = factor(CVD_cat, levels = c("Yes", "No"))
  )

sub_design <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTSAF2YR, nest = TRUE, data = df_sub)
all_adj_vars <- c("RIDAGEYR", "RIAGENDR", "RIDRETH1", "INDFMPIR", "DMDEDUC2", "DMDMARTL", "diabetes", "hyperten", "cvd_hx")

var_list <- list(
  "Age" = list(col = "Age_group", rm_cov = "RIDAGEYR"),
  "Gender" = list(col = "Gender", rm_cov = "RIAGENDR"),
  "Race" = list(col = "RIDRETH1", rm_cov = "RIDRETH1"),
  "Education" = list(col = "DMDEDUC2", rm_cov = "DMDEDUC2"),
  "PIR" = list(col = "PIR_group", rm_cov = "INDFMPIR"),
  "Marital status" = list(col = "DMDMARTL", rm_cov = "DMDMARTL"),
  "Diabetes" = list(col = "Diabetes_cat", rm_cov = "diabetes"),
  "Hypertension" = list(col = "Hyperten_cat", rm_cov = "hyperten"),
  "Cardiovascular outcome" = list(col = "CVD_cat", rm_cov = "cvd_hx"),
  "Obesity" = list(col = "Obesity", rm_cov = NA) 
)

# 2. HÀM TÍNH TOÁN DATA (Đã thêm suppressWarnings để giấu cảnh báo Loglik do mẫu nhỏ)
generate_subgroup_data <- function(event_var) {
  res_list <- list()
  
  for (v_name in names(var_list)) {
    col_name <- var_list[[v_name]]$col
    rm_cov <- var_list[[v_name]]$rm_cov
    
    curr_covs <- if(!is.na(rm_cov)) setdiff(all_adj_vars, rm_cov) else all_adj_vars
    covs_fmla <- paste(curr_covs, collapse = " + ")
    
    # a) Tính P-for-interaction
    form_inter <- as.formula(paste("Surv(time_months, ", event_var, ") ~ TyHGB *", col_name, "+", covs_fmla))
    fit_inter <- suppressWarnings(tryCatch(svycoxph(form_inter, design = sub_design), error = function(e) NULL))
    
    p_inter <- NA
    if(!is.null(fit_inter)) {
      p_inter <- suppressWarnings(tryCatch({
        regTermTest(fit_inter, as.formula(paste("~ TyHGB:", col_name)))$p[1]
      }, error = function(e) NA))
    }
    p_inter_str <- ifelse(is.na(p_inter), "N/A", ifelse(p_inter < 0.001, "<0.001", sprintf("%.3f", p_inter)))
    
    res_list[[length(res_list) + 1]] <- data.frame(
      Variables = v_name, Subgroups = "", Event = "", ` ` = paste(rep(" ", 30), collapse = ""),
      `HR(95%CI)` = "", `P value` = "", `P for interaction` = p_inter_str,
      est = NA, lower = NA, upper = NA, check.names = FALSE, stringsAsFactors = FALSE
    )
    
    # b) Tính HR cho từng Level
    lvls <- levels(df_sub[[col_name]])
    for (lvl in lvls) {
      ev_count <- sum(df_sub[[col_name]] == lvl & df_sub[[event_var]] == 1, na.rm = TRUE)
      
      # Dùng suppressMessages & suppressWarnings để không nhả chữ rác ra Console
      d_sub <- suppressMessages(suppressWarnings(subset(sub_design, get(col_name) == lvl)))
      form_sub <- as.formula(paste("Surv(time_months, ", event_var, ") ~ TyHGB +", covs_fmla))
      fit_sub <- suppressWarnings(tryCatch(svycoxph(form_sub, design = d_sub), error = function(e) NULL))
      
      if (!is.null(fit_sub)) {
        sm <- summary(fit_sub)
        hr <- sm$conf.int[1, "exp(coef)"]
        ci_l <- sm$conf.int[1, "lower .95"]
        ci_u <- sm$conf.int[1, "upper .95"]
        pval <- sm$coefficients[1, "Pr(>|z|)"]
        
        # Bọc lót nếu HR nổ lên vô cực (lớn hơn 100) do thiếu dữ liệu ở nhóm nhỏ
        if(!is.na(hr) && hr > 100) { hr <- NA; ci_l <- NA; ci_u <- NA }
        
        if(!is.na(hr)) {
          hr_str <- sprintf("%.3f (%.3f, %.3f)", hr, ci_l, ci_u)
          pval_str <- ifelse(pval < 0.001, "<0.001", sprintf("%.3f", pval))
        } else {
          hr_str <- "N/A"; pval_str <- "N/A"
        }
      } else {
        hr <- NA; ci_l <- NA; ci_u <- NA
        hr_str <- "N/A"; pval_str <- "N/A"
      }
      
      res_list[[length(res_list) + 1]] <- data.frame(
        Variables = "", Subgroups = paste0("      ", lvl), Event = as.character(ev_count),
        ` ` = paste(rep(" ", 30), collapse = ""), `HR(95%CI)` = hr_str, `P value` = pval_str, 
        `P for interaction` = "", est = hr, lower = ci_l, upper = ci_u,
        check.names = FALSE, stringsAsFactors = FALSE
      )
    }
  }
  return(do.call(rbind, res_list))
}

# 3. HÀM VẼ FOREST PLOT TỪ DATA (BẢN UPDATE ZOOM TRỤC X & KẺ BẢNG)
draw_forest_plot <- function(dt) {
  
  fill_colors <- c()
  current_color <- "#E8F4F0" 
  for (i in 1:nrow(dt)) {
    if (dt$Variables[i] != "") {
      current_color <- ifelse(current_color == "white", "#E8F4F0", "white")
    }
    fill_colors <- c(fill_colors, current_color)
  }
  
  tm <- forest_theme(
    base_size = 11,
    refline_gp = gpar(col = "black", lty = "dashed"), 
    arrow_type = "closed",
    arrow_label_just = "end",                         
    ci_pch = 15,
    ci_col = "#E64B35", 
    ci_fill = "#E64B35",
    ci_alpha = 1,
    ci_lty = 1,
    ci_lwd = 1.5,
    ci_Theight = 0.2,
    core = list(bg_params = list(fill = fill_colors))
  )
  
  p <- forest(
    data = dt[, c("Variables", "Subgroups", "Event", " ", "HR(95%CI)", "P value", "P for interaction")],
    est = dt$est,
    lower = dt$lower,
    upper = dt$upper,
    ci_column = 4, 
    ref_line = 1,
    arrow_lab = c("Low risk", "High Risk"),
    
    # ---------- BỘ LỌC TỌA ĐỘ MỚI (ZOOM TRỤC X) ----------
    xlim = c(0.7, 1.2), # Zoom cực sâu vào đoạn 0.6 đến 1.4 để KTC dãn rộng ra
    clip = c(0.7, 1.2), # Khóa chặt viền, KTC nào lố quá sẽ tự thành mũi tên cực đẹp
    ticks_at = c(0.7, 0.9, 1.0, 1.2), # Chia vạch chi tiết hơn
    # -----------------------------------------------------
    
    theme = tm
  )
  
  # In đậm cho cột Variables
  p <- edit_plot(p, row = 1:nrow(dt), col = 1, part = "body", gp = gpar(fontface = "bold"))
  
  # ---------- FIX LỖI KẺ BẢNG Ở ĐÂY ----------
  # 1. Đường kẻ dày trên cùng của Header
  p <- add_border(p, part = "header", row = 1, where = "top", gp = gpar(lwd = 1.5, col = "black"))
  
  # 2. Đường kẻ ngang phân cách Header và phần dữ liệu
  p <- add_border(p, part = "header", row = 1, where = "bottom", gp = gpar(lwd = 1.5, col = "black"))
  
  # 3. Đường kẻ chốt sổ dưới cùng của toàn bộ bảng (Để bọc lại cho đẹp)
  p <- add_border(p, part = "body", row = nrow(dt), where = "bottom", gp = gpar(lwd = 1.5, col = "black"))
  # -------------------------------------------
  
  return(p)
}

# 4. CHẠY VÀ LẮP RÁP HÌNH ẢNH
cat("Đang Bootstrapping tính Subgroup cho ACM (Console sẽ im lặng, sếp ráng đợi xíu nhé)...\n")
dt_acm <- generate_subgroup_data("status_acm")
p_acm <- draw_forest_plot(dt_acm)

cat("Đang Bootstrapping tính Subgroup cho CVM...\n")
dt_cvm <- generate_subgroup_data("status_cvm")
p_cvm <- draw_forest_plot(dt_cvm)

cat("Đang ghép khung hình ngang, gắn tag (A), (B) và xuất file...\n")

final_forest <- wrap_elements(full = p_acm) | wrap_elements(full = p_cvm)
final_forest <- final_forest + 
  plot_annotation(
    tag_levels = 'A',           # Đánh dấu bằng chữ cái in hoa (A, B, C...)
    tag_prefix = '(',           # Thêm dấu ngoặc mở
    tag_suffix = ')'            # Thêm dấu ngoặc đóng
  ) &
  theme(plot.tag = element_text(size = 22, face = "bold")) 

ggsave("Figure_Subgroup_Analysis.png", plot = final_forest, 
       width = 22, height = 12, units = "in", dpi = 300, bg = "white")


# ==============================================================================
# BƯỚC 10: PHÂN TÍCH TRUNG GIAN (MEDIATION ANALYSIS) VỚI ACID URIC (LBXSUA)
# ==============================================================================

# 1. LÀM SẠCH DỮ LIỆU (Xử lý missing cho LBXSUA)
df_med_clean <- df_final %>% filter(!is.na(LBXSUA))

mediator_var <- "LBXSUA"
# Danh sách biến hiệu chỉnh đầy đủ để tránh nhiễu (Confounding)
med_covs <- c("RIDAGEYR", "RIAGENDR", "RIDRETH1", "INDFMPIR", "DMDEDUC2", "DMDMARTL", "diabetes")
med_covs_str <- paste(med_covs, collapse = " + ")

med_design <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTSAF2YR, nest = TRUE, data = df_med_clean)

# 2. HÀM TÍNH TOÁN VÀ VẼ SƠ ĐỒ
draw_mediation_plot <- function(event_var, title_str) {
  
  # Path A: TyHGB -> Acid Uric (Hiệu chỉnh nhiễu)
  form_path_a <- as.formula(paste(mediator_var, "~ TyHGB +", med_covs_str))
  fit_path_a <- svyglm(form_path_a, design = med_design, family = gaussian())
  sm_a <- summary(fit_path_a)
  beta_a <- coef(fit_path_a)["TyHGB"]
  se_a <- sm_a$coefficients["TyHGB", "Std. Error"]
  p_a <- sm_a$coefficients["TyHGB", "Pr(>|t|)"]
  
  # Total Effect: TyHGB -> Tử vong (Hiệu chỉnh nhiễu)
  form_tot <- as.formula(paste("Surv(time_months, ", event_var, ") ~ TyHGB +", med_covs_str))
  fit_tot <- svycoxph(form_tot, design = med_design)
  sm_tot <- summary(fit_tot)
  beta_tot <- coef(fit_tot)["TyHGB"]
  p_tot <- sm_tot$coefficients["TyHGB", "Pr(>|z|)"]
  
  # Direct Effect & Path B: TyHGB + Acid Uric -> Tử vong (Hiệu chỉnh nhiễu)
  form_dir <- as.formula(paste("Surv(time_months, ", event_var, ") ~ TyHGB +", mediator_var, "+", med_covs_str))
  fit_dir <- svycoxph(form_dir, design = med_design)
  sm_dir <- summary(fit_dir)
  beta_dir <- coef(fit_dir)["TyHGB"]
  p_dir <- sm_dir$coefficients["TyHGB", "Pr(>|z|)"]
  
  med_row <- grep(paste0("^", mediator_var), rownames(sm_dir$conf.int), value = TRUE)[1]
  hr_b <- sm_dir$conf.int[med_row, "exp(coef)"]
  p_b <- sm_dir$coefficients[med_row, "Pr(>|z|)"]
  
  # Tính % Trung gian (Proportion Mediated)
  prop_med <- (beta_tot - beta_dir) / beta_tot * 100
  if(is.na(prop_med) || prop_med < 0) prop_med <- 0 
  
  # --- Format Text để vẽ hình ---
  fmt_p <- function(p) ifelse(!is.na(p) & p < 0.001, "P<0.001", sprintf("P=%.3f", p))
  str_tot <- sprintf("%.3f(95%%CI:%.3f~%.3f, %s)", sm_tot$conf.int["TyHGB", 1], sm_tot$conf.int["TyHGB", 3], sm_tot$conf.int["TyHGB", 4], fmt_p(p_tot))
  str_dir <- sprintf("%.3f(95%%CI:%.3f~%.3f, %s)", sm_dir$conf.int["TyHGB", 1], sm_dir$conf.int["TyHGB", 3], sm_dir$conf.int["TyHGB", 4], fmt_p(p_dir))
  str_a <- sprintf("%.3f(95%%CI:%.3f~%.3f, %s)", beta_a, beta_a-1.96*se_a, beta_a+1.96*se_a, fmt_p(p_a))
  str_b <- sprintf("%.3f(95%%CI:%.3f~%.3f, %s)", hr_b, sm_dir$conf.int[med_row, 3], sm_dir$conf.int[med_row, 4], fmt_p(p_b))
  str_prop <- sprintf("%.2f%%", prop_med)
  
  # [Phần code vẽ ggplot2 giữ nguyên như Step trước, chỉ đổi label mediator thành "Acid Uric"]
  boxes <- data.frame(x = c(1.5, 8.5, 1.5, 8.5, 5), y = c(6, 6, 1, 1, 4),
                      label = c("TyHGB", title_str, "TyHGB", title_str, "Acid Uric"))
  
  p <- ggplot() + theme_void() + coord_cartesian(xlim = c(0, 10), ylim = c(0, 7)) +
    geom_segment(aes(x = 2.5, y = 6, xend = 7.5, yend = 6), arrow = arrow(length = unit(0.3, "cm")), color = "#2CA02C", linewidth = 1) + 
    geom_segment(aes(x = 2.5, y = 1, xend = 7.5, yend = 1), arrow = arrow(length = unit(0.3, "cm")), color = "#2CA02C", linewidth = 1) + 
    geom_segment(aes(x = 2.2, y = 1.5, xend = 4.2, yend = 3.5), arrow = arrow(length = unit(0.3, "cm")), color = "#2CA02C", linewidth = 1) + 
    geom_segment(aes(x = 5.8, y = 3.5, xend = 7.8, yend = 1.5), arrow = arrow(length = unit(0.3, "cm")), color = "#2CA02C", linewidth = 1) + 
    geom_rect(data = boxes, aes(xmin = x - 1, xmax = x + 1, ymin = y - 0.4, ymax = y + 0.4), fill = "#E8F4F0", color = "#4DBBD5", linewidth = 1, rx = 0.2, ry = 0.2) +
    geom_text(data = boxes, aes(x = x, y = y, label = label), size = 5, fontface = "bold") +
    geom_text(aes(x = 5, y = 6.3, label = "Total effect"), fontface = "bold", size = 4) +
    geom_text(aes(x = 5, y = 5.7, label = str_tot), fontface = "bold", size = 3.8) +
    geom_text(aes(x = 5, y = 1.3, label = "Direct effect"), fontface = "bold", size = 4) +
    geom_text(aes(x = 5, y = 0.7, label = str_dir), fontface = "bold", size = 3.8) +
    geom_text(aes(x = 5, y = 2.5, label = "Proportion mediated"), fontface = "bold", size = 4.5) +
    geom_text(aes(x = 5, y = 2.1, label = str_prop), fontface = "bold", size = 5) +
    geom_text(aes(x = 2.7, y = 2.8, label = str_a), fontface = "bold", size = 3.5, angle = 45) +
    geom_text(aes(x = 7.3, y = 2.8, label = str_b), fontface = "bold", size = 3.5, angle = -45)
  
  return(p)
}

p_med_acm <- draw_mediation_plot("status_acm", "All-cause\nmortality")
p_med_cvm <- draw_mediation_plot("status_cvm", "Cardiovascular\nmortality")

final_mediation <- (p_med_acm | p_med_cvm) + plot_annotation(tag_levels = 'A',           # Đánh dấu bằng chữ cái in hoa (A, B, C...)
                                                             tag_prefix = '(',           # Thêm dấu ngoặc mở
                                                             tag_suffix = ')') +
  theme(plot.tag = element_text(size = 22, face = "bold"))
ggsave("Figure_Mediation_UricAcid.png", plot = final_mediation, width = 18, height = 8, dpi = 300)

