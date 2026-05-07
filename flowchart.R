pacman::p_load(readxl, tidyverse, mice, DiagrammeR)

raw_df <- read_excel("data.xlsx", sheet = "RAW DATA") 

# N_1: Tổng số ca ban đầu
n_raw <- nrow(raw_df)

# N_2: Sau khi lọc tuổi
df_age <- raw_df %>% filter(RIDAGEYR >= 20)
n_age <- nrow(df_age)

# N_3: Sau khi lọc các chỉ số cốt lõi và trọng số
df_step1 <- df_age %>% 
  filter(!is.na(LBXTR) & !is.na(LBDSGLSI) & !is.na(LBDHDD) & 
           !is.na(LBXSATSI) & !is.na(LBXSASSI) &                 
           !is.na(BMXBMI) & !is.na(BMXWAIST) & !is.na(BMXWT) & !is.na(BMXHT) & 
           !is.na(WTSAF2YR) & WTSAF2YR > 0) %>% 
  mutate(cvd_hx = ifelse(MCQ160B == 1 | MCQ160C == 1 | MCQ160E == 1 | MCQ160F == 1, 1, 0))
n_labs <- nrow(df_step1)

# --- XỬ LÝ MICE IMPUTATION ĐỂ TÍNH ĐÚNG BỆNH NỀN ---
vars_to_impute <- c("INDFMPIR", "DMDEDUC2", "DMDMARTL", "LBXWBCSI", "LBXSCR", 
                    "LBXSUA", "cvd_hx", "diabetes", "hyperten")
cat_vars <- c("RIAGENDR", "RIDRETH1", "DMDEDUC2", "DMDMARTL", "diabetes", "hyperten", "cvd_hx")
df_step1[cat_vars] <- lapply(df_step1[cat_vars], as.factor)

# Chạy MICE (chỉ 1 vòng cho nhanh)
imputed_data <- mice(df_step1[, vars_to_impute], m = 1, method = 'pmm', seed = 123, printFlag = FALSE)

df_imputed <- df_step1
df_imputed[, vars_to_impute] <- complete(imputed_data, 1)
df_imputed$diabetes <- as.numeric(as.character(df_imputed$diabetes))
df_imputed$hyperten <- as.numeric(as.character(df_imputed$hyperten))

# N_4: Chẩn đoán MASLD và Lọc dữ liệu sinh tồn
df_masld <- df_imputed %>%
  mutate(
    HSI = 8 * (LBXSATSI / LBXSASSI) + BMXBMI + ifelse(RIAGENDR == 2, 2, 0) + ifelse(diabetes == 1, 2, 0),
    steatosis = ifelse(HSI > 36, 1, 0),
    cm_risk = ifelse(BMXBMI >= 25 | (RIAGENDR == 1 & BMXWAIST >= 94) | (RIAGENDR == 2 & BMXWAIST >= 80) | 
                       LBDSGLSI >= 5.55 | diabetes == 1 | hyperten == 1 | LBXTR >= 150 | 
                       (RIAGENDR == 1 & LBDHDD < 40) | (RIAGENDR == 2 & LBDHDD < 50), 1, 0),
    is_masld = ifelse(steatosis == 1 & cm_risk == 1, 1, 0),
    status_cvm = ifelse(mortstat == 1 & ucod_leading %in% c("001", "005"), 1, 0),
    time_months = permth_int
  ) %>%
  filter(is_masld == 1, (is.na(RIDEXPRG) | RIDEXPRG != 1), !is.na(mortstat), !is.na(time_months), eligstat == 1)
n_masld <- nrow(df_masld)

# N_5: Lọc các ca từ chối/không biết trả lời (N cuối cùng)
df_final <- df_masld %>%
  mutate(
    DMDEDUC2 = case_when(
      DMDEDUC2 %in% c(1, 2) ~ "Under High School",
      DMDEDUC2 == 3 ~ "High School to University",
      DMDEDUC2 %in% c(4, 5) ~ "Above University",
      TRUE ~ NA_character_
    ),
    DMDMARTL = case_when(
      DMDMARTL %in% c(1, 6) ~ "Married/Partner", 
      DMDMARTL %in% c(2, 3, 4) ~ "Widowed/Divorced/Separated",
      DMDMARTL == 5 ~ "Never married", 
      TRUE ~ NA_character_
    ),
    RIDRETH1 = case_when(
      RIDRETH1 == 3 ~ "Non-Hispanic White",
      RIDRETH1 == 4 ~ "Non-Hispanic Black",
      TRUE ~ "Other Races"
    )
  ) %>%
  drop_na(DMDEDUC2, DMDMARTL, RIDRETH1)
n_final <- nrow(df_final)

# ==============================================================================
# 2. CHUẨN BỊ TEXT CHO CÁC KHỐI SƠ ĐỒ
# ==============================================================================
node1 <- paste0("NHANES participants\\nin raw dataset\\n(N = ", n_raw, ")")
excl1 <- paste0("Excluded: Age < 20 years\\n(n = ", n_raw - n_age, ")")

node2 <- paste0("Adult participants\\n(N = ", n_age, ")")
excl2 <- paste0("Excluded: Missing core\\nanthropometric/biochemical data\\nor survey weights\\n(n = ", n_age - n_labs, ")")

node3 <- paste0("Participants with complete\\ncore data\\n(N = ", n_labs, ")")
excl3 <- paste0("Excluded: Non-MASLD,\\nPregnant, or missing\\nsurvival data\\n(n = ", n_labs - n_masld, ")")

node4 <- paste0("MASLD participants with\\nsurvival data\\n(N = ", n_masld, ")")
excl4 <- paste0("Excluded: Missing/Refused\\nEducation, Marital Status,\\nor Race\\n(n = ", n_masld - n_final, ")")

node5 <- paste0("FINAL STUDY POPULATION\\nFOR ANALYSIS\\n(N = ", n_final, ")")

# ==============================================================================
# 3. VẼ ĐỒ HỌA GRAPHVIZ (FIX LỖI NGHIÊNG - CHUẨN PRISMA)
# ==============================================================================
graph_code <- paste0("
digraph flowchart {
  # splines = ortho tạo đường mũi tên gấp khúc 90 độ cực kỳ chuyên nghiệp
  graph [layout = dot, rankdir = TB, nodesep = 0.5, ranksep = 0.5, splines = ortho]
  
  # Cài đặt khối chính (Thêm group = 'main' để ép thành 1 cột dọc thẳng đứng)
  node [fontname = Helvetica, shape = box, style = filled, fillcolor = '#E6F2FF', color = '#0055A4', penwidth = 1.5, margin = '0.25,0.15']
  A [label = '", node1, "', group = 'main']
  B [label = '", node2, "', group = 'main']
  C [label = '", node3, "', group = 'main']
  D [label = '", node4, "', group = 'main']
  
  # Khối kết quả cuối cùng
  node [fillcolor = '#D5E8D4', color = '#82B366', fontcolor = black]
  E [label = '", node5, "', group = 'main']

  # Các khối bị loại trừ (Màu hồng đỏ)
  node [fillcolor = '#F8CECC', color = '#B85450', fontcolor = black]
  X1 [label = '", excl1, "']
  X2 [label = '", excl2, "']
  X3 [label = '", excl3, "']
  X4 [label = '", excl4, "']

  # Cài đặt mũi tên
  edge [color = '#333333', penwidth = 1.5, arrowhead = vee]

  # Liên kết trục dọc (Ép trọng số weight = 1000 để luôn luôn thẳng)
  A -> B [weight = 1000]
  B -> C [weight = 1000]
  C -> D [weight = 1000]
  D -> E [weight = 1000]

  # Liên kết nhánh loại trừ (Tự động bẻ góc 90 độ sang ngang)
  A -> X1 
  B -> X2 
  C -> X3 
  D -> X4 

  # Ép các khối màu đỏ nằm song song ngang hàng với khối màu xanh
  { rank = same; A; X1 }
  { rank = same; B; X2 }
  { rank = same; C; X3 }
  { rank = same; D; X4 }
}
")

# 4. VẼ FLOWCHART
grViz(graph_code)