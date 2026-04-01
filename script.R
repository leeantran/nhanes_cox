if(!require(pacman)) install.packages("pacman")
pacman::p_load(survey, survival, gtsummary, survminer, splines, mice, timeROC, 
               dplyr, ggplot2, ggsci, rms, tidycmprsk)

## 1. Chuẩn bị data
raw_df <- read_excel("data.xlsx", sheet = "Modified data")

df <- raw_df %>%
  # Chỉ lấy những ca đủ điều kiện và có data tử vong
  filter(eligstat == 1 & !is.na(mortstat) & !is.na(permth_int)) %>%
  # Tạo biến Tử vong (1: Chết, 0: Sống) và Tử vong tim mạch (CVM)
  mutate(
    all_cause_mortality = ifelse(mortstat == 2, 1, 0),
    # UCOD_LEADING = 1 (Bệnh tim) hoặc 5 (Mạch máu não) quy thành tử vong do tim mạch
    cvm_mortality = ifelse(ucod_leading %in% c(1, 5), 1, 0),
    # Biến trạng thái cho Competing Risk: 0 = Sống, 1 = Chết do CVM, 2 = Chết do NN khác
    competing_risk_stat = case_when(
      mortstat == 1 ~ 0,
      cvm_mortality == 1 ~ 1,
      mortstat == 2 & cvm_mortality == 0 ~ 2
    ),
    # Sử dụng LBXTC (Cholesterol toàn phần), LBXGLU (Glucose huyết tương), LBDHDD (HDL-c)
    chg_index = log((LBXTC * LBXSGL) / (2 * LBDHDD)),
    # Chia tứ phân vị (Quartiles) cho CHG
    CHG_Quartiles = ntile(chg_index, 4),
    CHG_Quartiles = factor(CHG_Quartiles, labels = c("Q1", "Q2", "Q3", "Q4"))
  )

df <- df %>% mutate(across(c(RIAGENDR, RIDRETH1, DMDEDUC2, DIQ010, 
                             BPQ020), as.factor)) # Không có SMQ040


# 2. Xử lý Missing Data (Multiple Imputation)
# Chỉ gán các biến hiệp biến, KHÔNG gán biến kết cục (tử vong) và thời gian
vars_to_impute <- c("LBXTC", "LBXSGL", "LBDHDD", "BMXBMI", "BMXWAIST", "LBXTR", 
                    "LBXSATSI", "LBXSASSI", "LBXGH", "BPXSY1", "BPXDI1")

imputed_data <- mice(df[, vars_to_impute], m=5, method='pmm', seed=123)
df_complete <- df
df_complete[, vars_to_impute] <- complete(imputed_data, 1)

# KHAI BÁO SURVEY DESIGN
nhanes_design <- svydesign(
  id = ~SDMVPSU, 
  strata = ~SDMVSTRA, 
  weights = ~WTMEC2YR, 
  data = df_complete, 
  nest = TRUE
)


## 3. Demographic table
table1 <- nhanes_design %>%
  tbl_svysummary(
    by = CHG_Quartiles,
    include = c(
      # Đặc điểm nền
      RIDAGEYR, RIAGENDR, RIDRETH1, DMDEDUC2, INDFMPIR, 
      SMQ040, ALQ121, BMXBMI, BMXWAIST, BMXHIP,
      # Bệnh nền
      BPXSY1, BPXDI1, DIQ010, BPQ020, MCQ160, MCQ220, KIQ022,
      # Cận lâm sàng
      LBXTC, LBDHDD, LBXTR, LBXGLU, LBXGH, LBXSATSI, LBXSASSI, LBXIN, LBXHSCRP,
      LBXWBCSI, LBXHGB, LBXPLTSI, LBXSCR, LBXSBU,
      # Outcome
      chg_index
    ),
    statistic = list(
      all_continuous() ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{n_unweighted} ({p}%)"
    ),
    missing = "no" 
  ) %>%
  add_p() %>%
  add_overall() %>%
  modify_header(label = "**Characteristics**") %>%
  bold_labels()

table1