################################################################################
#                    CAPRIA RESURVEY 2002-2024: ANALISI COMPLETE
#                    APPENDICE A.1 - Script R analisi principale
#                    Pipeline completa analisi vegetazionali e dendrometriche
#
#  Autore: Alberto Gozzi
#  Affiliazione: Università di Bologna - Distal
#  Email: alberto.gozzi@tudio.unibo.it
#  Data: Febbraio 2025
#  Versione: 3.0
#
################################################################################

# =============================================================================
# PANORAMICA DELLO SCRIPT
# =============================================================================
#
# OBIETTIVO:
#   Analizzare il cambiamento della vegetazione tra 2002 e 2024 nell'area
#   di studio di Capria (Appennino tosco-emiliano), integrando:
#   - Dati floristici (presenza/assenza e coperture)
#   - Tratti funzionali (forme biologiche, Ellenberg, AFS, corotipi)
#   - Dati dendrometrici (struttura soprassuolo arboreo)
#
# APPROCCIO METODOLOGICO:
#   - Inferenza primaria: presenza/assenza (distanza Jaccard)
#   - Inferenza di supporto: coperture (distanza Bray-Curtis)
#   - Controllo ricollocazione: test appaiati, permutazioni bloccate
#   - Standardizzazione area: 100 m² (con cautela interpretativa)
#
# SEZIONI PRINCIPALI:
#   0) Setup e caricamento pacchetti
#   1) Analisi univariate per plot (ricchezza, habitat, Ellenberg CWM)
#   2) Analisi multivariate (NMDS, PERMANOVA, betadisper, envfit)
#   2.5) Analisi Bray-Curtis supplementare (coperture)
#   3) Indicator Species Analysis (IndVal)
#   4) Analisi dendrometrica (densità, area basimetrica, distribuz. diametrica)
#   5) Generazione grafici (17 grafici PNG a 300 dpi)
#   6) Export tabelle (6 file CSV)
#
# PARAMETRI CHIAVE (riproducibilità):
#   - Seed globale: 123
#   - Permutazioni PERMANOVA: 999 (bloccate per PlotBase)
#   - Permutazioni betadisper: 999
#   - Permutazioni envfit: 999
#   - Permutazioni IndVal: 999 (bloccate per PlotBase)
#   - NMDS: k=2 dimensioni, trymax=200
#   - Classi diametriche: 5 cm
#
# FILE INPUT RICHIESTI (nella cartella DATA_DIR):
#   - plot_specie_coperture_strato.csv   (coperture floristiche)
#   - plot_2002_2024_stazioni.csv        (metadata plot)
#   - traits_completo.csv                (tratti funzionali specie)
#   - Dendrometria.csv                   (rilievi dendrometrici)
#
# OUTPUT GENERATI:
#   - Cartella figures/: 17 grafici PNG (300 dpi)
#   - Cartella tables/: 6 tabelle CSV
#   - Console: statistiche e risultati test
#
# TEMPO ESECUZIONE STIMATO: 2-3 minuti (PC standard)
#
# =============================================================================

################################################################################
# INIZIO CODICE
################################################################################

################################################################################
#                    CAPRIA RESURVEY 2002-2024: ANALISI COMPLETE
#                    Script R per tesi - Pipeline completa v3
#                    Include: Forme biologiche, Corotipi, AFS
################################################################################

# =============================================================================
# 0) SETUP E CARICAMENTO PACCHETTI
# =============================================================================

required_packages <- c("vegan", "dplyr", "tidyr", "ggplot2", "indicspecies", 
                       "permute", "patchwork", "RColorBrewer")

for(pkg in required_packages) {
  if(!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# Directory dati
DATA_DIR <- "C:/Users/alber/Desktop/ANALISIFINALISSISPERA"

# Directory output (necessari prima dei primi ggsave/write.csv)
fig_dir <- file.path(DATA_DIR, "figures")
if(!dir.exists(fig_dir)) dir.create(fig_dir)
tab_dir <- file.path(DATA_DIR, "tables")
if(!dir.exists(tab_dir)) dir.create(tab_dir)


# Funzione per leggere CSV con separatore automatico
read_csv_auto <- function(path) {
  ln <- readLines(path, n = 1, warn = FALSE)
  n_sem <- lengths(regmatches(ln, gregexpr(";", ln)))
  n_com <- lengths(regmatches(ln, gregexpr(",", ln)))
  if (n_sem > n_com) {
    read.csv2(path, stringsAsFactors = FALSE, na.strings = c("", "NA"))
  } else {
    read.csv(path, stringsAsFactors = FALSE, na.strings = c("", "NA"))
  }
}

# Funzione per standardizzare nomi specie
clean_species <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("_", " ", x)
  x <- gsub("\\s+", " ", x)
  sapply(strsplit(x, " "), function(parts) {
    if(length(parts) >= 1) {
      parts[1] <- paste0(toupper(substr(parts[1], 1, 1)), 
                         tolower(substr(parts[1], 2, nchar(parts[1]))))
    }
    if(length(parts) >= 2) {
      parts[2] <- tolower(parts[2])
    }
    paste(parts, collapse = " ")
  })
}

# Funzione per estrarre Anno e PlotBase da PlotID (formato "2002_1")
parse_plotid <- function(plotid) {
  parts <- strsplit(as.character(plotid), "_")
  data.frame(
    Anno = as.integer(sapply(parts, `[`, 1)),
    PlotBase = as.integer(sapply(parts, `[`, 2))
  )
}

cat("=== Setup completato ===\n\n")

# =============================================================================
# 0.1) CARICAMENTO DATI GREZZI
# =============================================================================

cat("--- Caricamento dati ---\n")

# File coperture (PlotID;Species;Stratus;Cover)
coperture <- read_csv_auto(file.path(DATA_DIR, "plot_specie_coperture_strato.csv"))
names(coperture) <- c("PlotID", "Species", "Stratus", "Cover")

# Estrai Anno e PlotBase
coperture <- cbind(coperture, parse_plotid(coperture$PlotID))

# File stazioni
stazioni <- read_csv_auto(file.path(DATA_DIR, "plot_2002_2024_stazioni.csv"))
stazioni <- cbind(stazioni, parse_plotid(stazioni$PlotID))

# File traits
traits <- read_csv_auto(file.path(DATA_DIR, "traits_completo.csv"))

# File dendrometria
dendro <- read_csv_auto(file.path(DATA_DIR, "Dendrometria.csv"))

cat("File caricati:\n")
cat("  - Coperture:", nrow(coperture), "righe\n")
cat("  - Stazioni:", nrow(stazioni), "righe\n")
cat("  - Traits:", nrow(traits), "specie\n")
cat("  - Dendrometria:", nrow(dendro), "alberi\n\n")

# =============================================================================
# 0.2) STANDARDIZZAZIONE TASSONOMICA
# =============================================================================

cat("--- Standardizzazione tassonomica ---\n")

coperture$Species <- clean_species(coperture$Species)
traits$Species <- clean_species(traits$Species)

sp_coperture <- unique(coperture$Species[!is.na(coperture$Species) & coperture$Species != ""])
sp_traits <- unique(traits$Species)
sp_mancanti <- setdiff(sp_coperture, sp_traits)

if(length(sp_mancanti) > 0) {
  cat("ATTENZIONE: Specie in coperture non presenti in traits:\n")
  print(sp_mancanti)
} else {
  cat("Tutte le specie hanno corrispondenza in traits\n")
}
cat("\n")

# =============================================================================
# 0.3) COMBINE COVER - Aggregazione strati per specie/plot
# =============================================================================

cat("--- Combine cover (aggregazione strati) ---\n")

coperture <- coperture[!is.na(coperture$Species) & coperture$Species != "", ]
coperture$Cover <- as.numeric(coperture$Cover)
coperture <- coperture[!is.na(coperture$Cover), ]

coperture_combined <- coperture %>%
  group_by(PlotID, Anno, PlotBase, Species) %>%
  summarise(Cover = sum(Cover, na.rm = TRUE), .groups = "drop") %>%
  mutate(Cover = pmin(Cover, 100))

cat("Righe originali:", nrow(coperture), "\n")
cat("Righe dopo combine:", nrow(coperture_combined), "\n\n")

# =============================================================================
# 0.4) COSTRUZIONE MATRICI COMUNITÀ
# =============================================================================

cat("--- Costruzione matrici comunità ---\n")

comm_raw_wide <- coperture_combined %>%
  select(PlotID, Species, Cover) %>%
  pivot_wider(names_from = Species, values_from = Cover, values_fill = 0)

rownames_comm <- comm_raw_wide$PlotID
comm_raw <- as.matrix(comm_raw_wide[, -1])
rownames(comm_raw) <- rownames_comm

comm_pa <- (comm_raw > 0) * 1

cat("Dimensioni comm_raw:", nrow(comm_raw), "plot x", ncol(comm_raw), "specie\n")
cat("Dimensioni comm_pa:", nrow(comm_pa), "plot x", ncol(comm_pa), "specie\n\n")

# =============================================================================
# 0.5) PREPARAZIONE METADATI PLOT E DESIGN APPAIATO
# =============================================================================

cat("--- Preparazione design appaiato ---\n")

plots_meta <- data.frame(PlotID = rownames(comm_raw), stringsAsFactors = FALSE)
plots_meta <- cbind(plots_meta, parse_plotid(plots_meta$PlotID))

# Verifica coppie appaiate
plot_counts <- plots_meta %>%
  group_by(PlotBase) %>%
  summarise(n_anni = n(), anni = paste(sort(unique(Anno)), collapse = ","), .groups = "drop")

paired_plots <- plot_counts %>% filter(n_anni == 2, anni == "2002,2024")
cat("Plot con coppia completa 2002-2024:", nrow(paired_plots), "\n")
cat("PlotBase appaiati:", paste(paired_plots$PlotBase, collapse = ", "), "\n\n")

# Filtra solo plot appaiati
plots_paired <- plots_meta %>% filter(PlotBase %in% paired_plots$PlotBase)
plots_paired <- plots_paired[order(plots_paired$PlotBase, plots_paired$Anno), ]

comm_raw_p <- comm_raw[plots_paired$PlotID, ]
comm_pa_p <- comm_pa[plots_paired$PlotID, ]

comm_raw_p <- comm_raw_p[, colSums(comm_raw_p) > 0]
comm_pa_p <- comm_pa_p[, colSums(comm_pa_p) > 0]

cat("Dopo filtering paired:\n")
cat("  comm_raw_p:", nrow(comm_raw_p), "plot x", ncol(comm_raw_p), "specie\n")
cat("  comm_pa_p:", nrow(comm_pa_p), "plot x", ncol(comm_pa_p), "specie\n\n")

# Aggiungi metadati stazioni
plots_paired <- plots_paired %>%
  left_join(stazioni %>% select(PlotID, Surface_m2, Arb_cover, Herb_cover), by = "PlotID")

# Converti Surface_m2
plots_paired$Surface_m2 <- as.numeric(gsub(",", ".", plots_paired$Surface_m2))

cat("=== SEZIONE 0 COMPLETATA ===\n\n")

# =============================================================================
# 1) ANALISI UNIVARIATE PER PLOT
# =============================================================================

cat("###############################################\n")
cat("# SEZIONE 1: ANALISI UNIVARIATE\n")
cat("###############################################\n\n")

# -----------------------------------------------------------------------------
# 1.1) Ricchezza specifica (S)
# -----------------------------------------------------------------------------

cat("--- 1.1 Ricchezza specifica ---\n")

plots_paired$Richness <- rowSums(comm_pa_p)

rich_summary <- plots_paired %>%
  group_by(Anno) %>%
  summarise(
    n = n(), mean_S = mean(Richness), sd_S = sd(Richness),
    median_S = median(Richness), min_S = min(Richness), max_S = max(Richness),
    .groups = "drop"
  )

cat("\nRicchezza per anno:\n")
print(rich_summary)

# Test Wilcoxon appaiato
rich_wide <- plots_paired %>%
  select(PlotBase, Anno, Richness) %>%
  pivot_wider(names_from = Anno, values_from = Richness, names_prefix = "Y")

wilcox_rich <- wilcox.test(rich_wide$Y2002, rich_wide$Y2024, paired = TRUE, exact = FALSE)

cat("\nWilcoxon paired test (Ricchezza):\n")
print(wilcox_rich)

rich_wide$delta_S <- rich_wide$Y2024 - rich_wide$Y2002
cat("\nDelta ricchezza (2024 - 2002):\n")
cat("  Media:", round(mean(rich_wide$delta_S), 2), "\n")
cat("  Mediana:", median(rich_wide$delta_S), "\n\n")

# -----------------------------------------------------------------------------
# 1.2) Ricchezza standardizzata a 100 m² (CAUTELA)
# -----------------------------------------------------------------------------

cat("--- 1.2 Ricchezza standardizzata a 100 m² ---\n")
cat("NOTA: Nel 2002 le aree non erano standard (50-120 m²) - interpretare con cautela\n\n")

plots_paired$Richness_100 <- plots_paired$Richness * (100 / plots_paired$Surface_m2)

rich100_summary <- plots_paired %>%
  group_by(Anno) %>%
  summarise(mean_S100 = mean(Richness_100, na.rm = TRUE), 
            sd_S100 = sd(Richness_100, na.rm = TRUE), .groups = "drop")

cat("Ricchezza standardizzata per anno:\n")
print(rich100_summary)

rich100_wide <- plots_paired %>%
  select(PlotBase, Anno, Richness_100) %>%
  pivot_wider(names_from = Anno, values_from = Richness_100, names_prefix = "Y")

if(sum(!is.na(rich100_wide$Y2002)) >= 3 && sum(!is.na(rich100_wide$Y2024)) >= 3) {
  wilcox_rich100 <- wilcox.test(rich100_wide$Y2002, rich100_wide$Y2024, paired = TRUE, exact = FALSE)
  cat("\nWilcoxon paired test (Ricchezza standardizzata):\n")
  print(wilcox_rich100)
}

# -----------------------------------------------------------------------------
# 1.3) Scomposizione: specie APERTO vs CHIUSO
# -----------------------------------------------------------------------------

cat("\n--- 1.3 Specie habitat aperto vs chiuso ---\n")

traits$L_num <- suppressWarnings(as.numeric(traits$L))

sp_aperto <- traits$Species[!is.na(traits$L_num) & traits$L_num >= 7]
sp_chiuso <- traits$Species[!is.na(traits$L_num) & traits$L_num <= 4]

cat("Specie habitat aperto (L >= 7):", length(sp_aperto), "\n")
cat("Specie habitat chiuso (L <= 4):", length(sp_chiuso), "\n")

count_group <- function(comm_mat, sp_list) {
  sp_present <- intersect(colnames(comm_mat), sp_list)
  if(length(sp_present) == 0) return(rep(0, nrow(comm_mat)))
  rowSums(comm_mat[, sp_present, drop = FALSE] > 0)
}

cover_group <- function(comm_mat, sp_list) {
  sp_present <- intersect(colnames(comm_mat), sp_list)
  if(length(sp_present) == 0) return(rep(0, nrow(comm_mat)))
  rowSums(comm_mat[, sp_present, drop = FALSE])
}

plots_paired$N_aperto <- count_group(comm_pa_p, sp_aperto)
plots_paired$N_chiuso <- count_group(comm_pa_p, sp_chiuso)
plots_paired$Cover_aperto <- cover_group(comm_raw_p, sp_aperto)
plots_paired$Cover_chiuso <- cover_group(comm_raw_p, sp_chiuso)

habitat_summary <- plots_paired %>%
  group_by(Anno) %>%
  summarise(
    mean_N_aperto = mean(N_aperto), mean_N_chiuso = mean(N_chiuso),
    mean_Cov_aperto = mean(Cover_aperto), mean_Cov_chiuso = mean(Cover_chiuso),
    .groups = "drop"
  )

cat("\nSpecie per gruppo habitat:\n")
print(habitat_summary)

aperto_wide <- plots_paired %>%
  select(PlotBase, Anno, N_aperto) %>%
  pivot_wider(names_from = Anno, values_from = N_aperto, names_prefix = "Y")

chiuso_wide <- plots_paired %>%
  select(PlotBase, Anno, N_chiuso) %>%
  pivot_wider(names_from = Anno, values_from = N_chiuso, names_prefix = "Y")

cat("\nTest Wilcoxon - Specie aperto:\n")
print(wilcox.test(aperto_wide$Y2002, aperto_wide$Y2024, paired = TRUE, exact = FALSE))

cat("\nTest Wilcoxon - Specie chiuso:\n")
print(wilcox.test(chiuso_wide$Y2002, chiuso_wide$Y2024, paired = TRUE, exact = FALSE))

aperto_wide$delta <- aperto_wide$Y2024 - aperto_wide$Y2002
chiuso_wide$delta <- chiuso_wide$Y2024 - chiuso_wide$Y2002

cat("\nDelta specie aperto: media =", round(mean(aperto_wide$delta), 2), "\n")
cat("Delta specie chiuso: media =", round(mean(chiuso_wide$delta), 2), "\n")

# -----------------------------------------------------------------------------
# 1.4) Ancient Forest Species (AFS)
# -----------------------------------------------------------------------------

cat("\n--- 1.4 Ancient Forest Species (AFS) ---\n")

traits$AFS_num <- suppressWarnings(as.numeric(traits$AFS))
sp_afs <- traits$Species[!is.na(traits$AFS_num) & traits$AFS_num == 1]

cat("Specie AFS nel dataset:", length(sp_afs), "\n")
if(length(sp_afs) > 0) cat("Lista AFS:", paste(sp_afs, collapse = ", "), "\n")

plots_paired$N_AFS <- count_group(comm_pa_p, sp_afs)
plots_paired$Cover_AFS <- cover_group(comm_raw_p, sp_afs)

afs_summary <- plots_paired %>%
  group_by(Anno) %>%
  summarise(mean_N_AFS = mean(N_AFS), mean_Cover_AFS = mean(Cover_AFS), .groups = "drop")

cat("\nAFS per anno:\n")
print(afs_summary)

if(sum(plots_paired$N_AFS) > 0) {
  afs_wide <- plots_paired %>%
    select(PlotBase, Anno, N_AFS) %>%
    pivot_wider(names_from = Anno, values_from = N_AFS, names_prefix = "Y")
  
  cat("\nTest Wilcoxon - N AFS:\n")
  print(wilcox.test(afs_wide$Y2002, afs_wide$Y2024, paired = TRUE, exact = FALSE))
}

# -----------------------------------------------------------------------------
# 1.5) Spettri forme biologiche (Raunkiaer)
# -----------------------------------------------------------------------------

cat("\n--- 1.5 Spettri forme biologiche (Raunkiaer) ---\n")

calc_spectrum_weighted <- function(comm_mat, traits_df, trait_col) {
  sp_comm <- colnames(comm_mat)
  trait_vec <- traits_df[[trait_col]][match(sp_comm, traits_df$Species)]
  
  results <- list()
  for(i in 1:nrow(comm_mat)) {
    abund <- comm_mat[i, ]
    ok <- !is.na(trait_vec) & abund > 0
    if(sum(ok) > 0) {
      df_temp <- data.frame(trait = trait_vec[ok], cover = abund[ok])
      spectrum <- df_temp %>%
        group_by(trait) %>%
        summarise(cover_tot = sum(cover), .groups = "drop") %>%
        mutate(perc = cover_tot / sum(cover_tot) * 100)
      results[[i]] <- spectrum
    } else {
      results[[i]] <- data.frame(trait = NA, cover_tot = 0, perc = 0)
    }
  }
  results
}

spectra_list <- calc_spectrum_weighted(comm_raw_p, traits, "Forma_biologica")

forms <- c("P", "Ch", "H", "G", "T", "Np")

for(f in forms) {
  plots_paired[[paste0("FB_", f)]] <- sapply(spectra_list, function(sp) {
    if(is.null(sp) || all(is.na(sp$trait))) return(0)
    val <- sp$perc[sp$trait == f]
    if(length(val) == 0) return(0)
    return(val)
  })
}

fb_summary <- plots_paired %>%
  group_by(Anno) %>%
  summarise(across(starts_with("FB_"), mean, na.rm = TRUE), .groups = "drop")

cat("\nSpettro forme biologiche (% copertura media):\n")
print(fb_summary)

cat("\nTest appaiati per forma biologica:\n")
for(f in forms) {
  var_name <- paste0("FB_", f)
  fb_wide <- plots_paired %>%
    select(PlotBase, Anno, all_of(var_name)) %>%
    pivot_wider(names_from = Anno, values_from = all_of(var_name), names_prefix = "Y")
  
  if(sum(fb_wide$Y2002, na.rm = TRUE) > 0 || sum(fb_wide$Y2024, na.rm = TRUE) > 0) {
    test_res <- wilcox.test(fb_wide$Y2002, fb_wide$Y2024, paired = TRUE, exact = FALSE)
    delta <- mean(fb_wide$Y2024 - fb_wide$Y2002, na.rm = TRUE)
    cat(sprintf("  %s: Delta = %+.1f%%, p = %.3f\n", f, delta, test_res$p.value))
  }
}

# -----------------------------------------------------------------------------
# 1.6) Spettri corotipi
# -----------------------------------------------------------------------------

cat("\n--- 1.6 Spettri corotipi ---\n")

traits$Corotipo_macro <- case_when(
  grepl("Europ|Eurasiat|Eurosib|Circumbor|Paleotemp", traits$Corotipo, ignore.case = TRUE) ~ "Temperato",
  grepl("Medit|Stenomedit|Eurimedit", traits$Corotipo, ignore.case = TRUE) ~ "Mediterraneo",
  grepl("Endem", traits$Corotipo, ignore.case = TRUE) ~ "Endemico",
  grepl("Orof", traits$Corotipo, ignore.case = TRUE) ~ "Orofitico",
  grepl("Cosmop|Subcosmop", traits$Corotipo, ignore.case = TRUE) ~ "Cosmopolita",
  TRUE ~ "Altro"
)

spectra_coro <- calc_spectrum_weighted(comm_raw_p, traits, "Corotipo_macro")

coro_cats <- c("Temperato", "Mediterraneo", "Endemico", "Orofitico", "Cosmopolita")

for(cc in coro_cats) {
  plots_paired[[paste0("Coro_", cc)]] <- sapply(spectra_coro, function(sp) {
    if(is.null(sp) || all(is.na(sp$trait))) return(0)
    val <- sp$perc[sp$trait == cc]
    if(length(val) == 0) return(0)
    return(val)
  })
}

coro_summary <- plots_paired %>%
  group_by(Anno) %>%
  summarise(across(starts_with("Coro_"), mean, na.rm = TRUE), .groups = "drop")

cat("\nSpettro corotipi (% copertura media):\n")
print(coro_summary)

# -----------------------------------------------------------------------------
# 1.7) Ellenberg CWM
# -----------------------------------------------------------------------------

cat("\n--- 1.7 Ellenberg CWM ---\n")

for(v in c("L", "T", "U", "R", "N", "C", "S")) {
  traits[[paste0(v, "_num")]] <- suppressWarnings(as.numeric(traits[[v]]))
}

calc_cwm <- function(comm_mat, trait_df, trait_col) {
  sp_comm <- colnames(comm_mat)
  trait_vec <- trait_df[[trait_col]][match(sp_comm, trait_df$Species)]
  
  cwm <- numeric(nrow(comm_mat))
  for(i in 1:nrow(comm_mat)) {
    abund <- comm_mat[i, ]
    ok <- !is.na(trait_vec) & abund > 0
    if(sum(ok) > 0) {
      cwm[i] <- sum(abund[ok] * trait_vec[ok]) / sum(abund[ok])
    } else {
      cwm[i] <- NA
    }
  }
  cwm
}

plots_paired$Ell_L <- calc_cwm(comm_raw_p, traits, "L_num")
plots_paired$Ell_U <- calc_cwm(comm_raw_p, traits, "U_num")
plots_paired$Ell_N <- calc_cwm(comm_raw_p, traits, "N_num")

ell_summary <- plots_paired %>%
  group_by(Anno) %>%
  summarise(
    mean_L = mean(Ell_L, na.rm = TRUE), sd_L = sd(Ell_L, na.rm = TRUE),
    mean_U = mean(Ell_U, na.rm = TRUE), sd_U = sd(Ell_U, na.rm = TRUE),
    mean_N = mean(Ell_N, na.rm = TRUE), sd_N = sd(Ell_N, na.rm = TRUE),
    .groups = "drop"
  )

cat("\nEllenberg CWM per anno:\n")
print(ell_summary)

cat("\nTest appaiati Ellenberg:\n")
for(var in c("Ell_L", "Ell_U", "Ell_N")) {
  var_wide <- plots_paired %>%
    select(PlotBase, Anno, all_of(var)) %>%
    pivot_wider(names_from = Anno, values_from = all_of(var), names_prefix = "Y")
  
  test_res <- wilcox.test(var_wide$Y2002, var_wide$Y2024, paired = TRUE, exact = FALSE)
  delta_val <- mean(var_wide$Y2024 - var_wide$Y2002, na.rm = TRUE)
  cat(sprintf("  %s: Delta = %+.3f, p = %.3f\n", var, delta_val, test_res$p.value))
}

cat("\n=== SEZIONE 1 COMPLETATA ===\n\n")

# =============================================================================
# 2) ANALISI MULTIVARIATE
# =============================================================================

cat("###############################################\n")
cat("# SEZIONE 2: ANALISI MULTIVARIATE\n")
cat("###############################################\n\n")

cat("--- 2.1 NMDS su Jaccard binaria ---\n")

dist_jac <- vegdist(comm_pa_p, method = "jaccard", binary = TRUE)

set.seed(123)
nmds <- metaMDS(dist_jac, k = 2, trymax = 200, trace = 0)

cat("Stress NMDS:", round(nmds$stress, 4), "\n")
cat("  < 0.05: eccellente, < 0.10: buono, < 0.20: accettabile\n\n")

nmds_scores <- as.data.frame(scores(nmds, display = "sites"))
nmds_scores$PlotID <- rownames(nmds_scores)
nmds_scores <- cbind(nmds_scores, parse_plotid(nmds_scores$PlotID))

cat("--- 2.2 PERMANOVA ---\n")

perm_design <- how(
  within = Within(type = "free"),
  plots = Plots(strata = factor(plots_paired$PlotBase), type = "none"),
  nperm = 999
)

adonis_res <- adonis2(dist_jac ~ Anno, data = plots_paired, permutations = perm_design)

cat("\nPERMANOVA (permutazioni bloccate per PlotBase):\n")
print(adonis_res)

cat("\n--- 2.3 BETADISPER ---\n")

bd <- betadisper(dist_jac, factor(plots_paired$Anno))
bd_test <- anova(bd)

cat("\nTest dispersione (ANOVA):\n")
print(bd_test)

cat("\nDispersione media per anno:\n")
print(tapply(bd$distances, plots_paired$Anno, mean))

if(bd_test$`Pr(>F)`[1] > 0.05) {
  cat("\n--> Dispersione OMOGENEA: PERMANOVA interpretabile\n")
} else {
  cat("\n--> ATTENZIONE: Dispersione NON omogenea\n")
}

cat("\n--- 2.4 Envfit ---\n")

envfit_vars <- plots_paired %>% select(Ell_L, Ell_U, Ell_N) %>% as.data.frame()

set.seed(123)
ef <- envfit(nmds, envfit_vars, permutations = 999, na.rm = TRUE)

cat("\nEnvfit results:\n")
print(ef)

cat("\n=== SEZIONE 2 COMPLETATA ===\n\n")

# =============================================================================
# 2.5) ANALISI MULTIVARIATE - BRAY-CURTIS (supplementare)
# =============================================================================
# NOTA METODOLOGICA: Le aree di campionamento differiscono tra 2002 (50-120 m²) 
# e 2024 (100 m² standardizzati). Poiché le coperture sono percentuali, 
# Bray-Curtis è meno sensibile a questo bias rispetto a metriche di abbondanza 
# assoluta, ma i risultati vanno interpretati con cautela.
# =============================================================================

cat("\n--- 2.5 NMDS e PERMANOVA su Bray-Curtis (coperture) ---\n")
cat("NOTA: Analisi supplementare - interpretare con cautela per differenze aree campionamento\n\n")

# Distanza Bray-Curtis
dist_bc <- vegdist(comm_raw_p, method = "bray")

# NMDS Bray-Curtis
set.seed(123)
nmds_bc <- metaMDS(dist_bc, k = 2, trymax = 200, trace = 0)

cat("Stress NMDS (Bray-Curtis):", round(nmds_bc$stress, 4), "\n")

# Scores NMDS
nmds_bc_scores <- as.data.frame(scores(nmds_bc, display = "sites"))
nmds_bc_scores$PlotID <- rownames(nmds_bc_scores)
nmds_bc_scores <- cbind(nmds_bc_scores, parse_plotid(nmds_bc_scores$PlotID))

# PERMANOVA Bray-Curtis
perm_design_bc <- how(
  within = Within(type = "free"),
  plots = Plots(strata = factor(plots_paired$PlotBase), type = "none"),
  nperm = 999
)

adonis_bc <- adonis2(dist_bc ~ Anno, data = plots_paired, permutations = perm_design_bc)

cat("\nPERMANOVA Bray-Curtis (permutazioni bloccate per PlotBase):\n")
print(adonis_bc)

# BETADISPER Bray-Curtis
bd_bc <- betadisper(dist_bc, factor(plots_paired$Anno))
bd_bc_test <- anova(bd_bc)

cat("\nTest dispersione Bray-Curtis (ANOVA):\n")
print(bd_bc_test)

cat("\nDispersione media per anno (Bray-Curtis):\n")
print(tapply(bd_bc$distances, plots_paired$Anno, mean))

if(bd_bc_test$`Pr(>F)`[1] > 0.05) {
  cat("\n--> Dispersione OMOGENEA: PERMANOVA interpretabile\n")
} else {
  cat("\n--> ATTENZIONE: Dispersione NON omogenea - interpretare con cautela\n")
}

# Envfit su Bray-Curtis
set.seed(123)
ef_bc <- envfit(nmds_bc, envfit_vars, permutations = 999, na.rm = TRUE)

cat("\nEnvfit results (Bray-Curtis):\n")
print(ef_bc)

# =============================================================================
# CONFRONTO JACCARD vs BRAY-CURTIS
# =============================================================================

cat("\n")
cat("=======================================================\n")
cat("CONFRONTO JACCARD (P/A) vs BRAY-CURTIS (Coperture)\n")
cat("=======================================================\n\n")

comparison_df <- data.frame(
  Metrica = c("Stress NMDS", "PERMANOVA R²", "PERMANOVA p", "Betadisper p"),
  Jaccard = c(
    round(nmds$stress, 3),
    round(adonis_res$R2[1], 3),
    round(adonis_res$`Pr(>F)`[1], 3),
    round(bd_test$`Pr(>F)`[1], 3)
  ),
  Bray_Curtis = c(
    round(nmds_bc$stress, 3),
    round(adonis_bc$R2[1], 3),
    round(adonis_bc$`Pr(>F)`[1], 3),
    round(bd_bc_test$`Pr(>F)`[1], 3)
  )
)

print(comparison_df)

cat("\nInterpretazione:\n")
cat("- Se entrambi significativi: cambiamento robusto in composizione E struttura\n")
cat("- Se solo Jaccard sign.: cambiano le specie rare, non le dominanti\n")
cat("- Se solo Bray-Curtis sign.: cambiano le abbondanze, non la lista specie\n")

cat("\n=== SEZIONE 2.5 COMPLETATA ===\n\n")


# =============================================================================
# GRAFICI BRAY-CURTIS
# =============================================================================

cat("--- Grafici Bray-Curtis ---\n")

# Palette colori
col_2002 <- "#E69F00"
col_2024 <- "#0072B2"

# Prepara dati wide per frecce
nmds_bc_wide <- nmds_bc_scores %>%
  select(PlotBase, Anno, NMDS1, NMDS2) %>%
  pivot_wider(names_from = Anno, values_from = c(NMDS1, NMDS2))

# NMDS Bray-Curtis con traiettorie
p_nmds_bc <- ggplot() +
  geom_segment(data = nmds_bc_wide,
               aes(x = NMDS1_2002, y = NMDS2_2002, 
                   xend = NMDS1_2024, yend = NMDS2_2024),
               alpha = 0.4, color = "grey50") +
  geom_point(data = nmds_bc_scores %>% filter(Anno == 2002),
             aes(x = NMDS1, y = NMDS2), 
             shape = 21, fill = col_2002, size = 3.5, stroke = 0.5) +
  geom_point(data = nmds_bc_scores %>% filter(Anno == 2024),
             aes(x = NMDS1, y = NMDS2), 
             shape = 22, fill = col_2024, size = 3.5, stroke = 0.5) +
  geom_text(data = nmds_bc_scores,
            aes(x = NMDS1, y = NMDS2, label = PlotBase),
            size = 2.2, vjust = -0.7, alpha = 0.7) +
  labs(title = paste0("NMDS (Bray-Curtis) - Stress = ", round(nmds_bc$stress, 3)),
       subtitle = paste0("PERMANOVA R² = ", round(adonis_bc$R2[1], 2), 
                         ", p = ", round(adonis_bc$`Pr(>F)`[1], 3)),
       x = "NMDS1", y = "NMDS2") +
  annotate("point", x = max(nmds_bc_scores$NMDS1) * 0.85, y = max(nmds_bc_scores$NMDS2) * 0.95, 
           shape = 21, fill = col_2002, size = 4) +
  annotate("text", x = max(nmds_bc_scores$NMDS1) * 0.92, y = max(nmds_bc_scores$NMDS2) * 0.95, 
           label = "2002", hjust = 0, size = 3.5) +
  annotate("point", x = max(nmds_bc_scores$NMDS1) * 0.85, y = max(nmds_bc_scores$NMDS2) * 0.80, 
           shape = 22, fill = col_2024, size = 4) +
  annotate("text", x = max(nmds_bc_scores$NMDS1) * 0.92, y = max(nmds_bc_scores$NMDS2) * 0.80, 
           label = "2024", hjust = 0, size = 3.5) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        aspect.ratio = 1)

ggsave(file.path(fig_dir, "15_nmds_braycurtis.png"), p_nmds_bc, width = 8, height = 8, dpi = 300)

# NMDS Bray-Curtis con envfit
ef_bc_vectors <- as.data.frame(scores(ef_bc, display = "vectors"))
ef_bc_vectors$var <- rownames(ef_bc_vectors)
arrow_scale <- 0.8

p_nmds_bc_envfit <- ggplot() +
  geom_point(data = nmds_bc_scores,
             aes(x = NMDS1, y = NMDS2, shape = factor(Anno), fill = factor(Anno)),
             size = 3, stroke = 0.5) +
  scale_shape_manual(values = c("2002" = 21, "2024" = 22)) +
  scale_fill_manual(values = c("2002" = col_2002, "2024" = col_2024)) +
  geom_segment(data = ef_bc_vectors,
               aes(x = 0, y = 0, xend = NMDS1 * arrow_scale, yend = NMDS2 * arrow_scale),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "black", linewidth = 1) +
  geom_text(data = ef_bc_vectors,
            aes(x = NMDS1 * arrow_scale * 1.15, y = NMDS2 * arrow_scale * 1.15, label = var),
            size = 4, fontface = "bold") +
  labs(title = "NMDS Bray-Curtis con gradienti Ellenberg",
       subtitle = paste0("Envfit: L r² = ", round(ef_bc$vectors$r["Ell_L"], 2),
                         ", U r² = ", round(ef_bc$vectors$r["Ell_U"], 2),
                         ", N r² = ", round(ef_bc$vectors$r["Ell_N"], 2)),
       x = "NMDS1", y = "NMDS2", shape = "Anno", fill = "Anno") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        aspect.ratio = 1)

ggsave(file.path(fig_dir, "16_nmds_braycurtis_envfit.png"), p_nmds_bc_envfit, width = 8, height = 7, dpi = 300)

# Betadisper Bray-Curtis
bd_bc_df <- data.frame(
  Anno = factor(plots_paired$Anno),
  Distance = bd_bc$distances,
  PlotBase = plots_paired$PlotBase
)

p_betadisper_bc <- ggplot(bd_bc_df, aes(x = Anno, y = Distance)) +
  geom_boxplot(aes(fill = Anno), alpha = 0.7, outlier.shape = NA, width = 0.5) +
  geom_point(shape = 21, size = 2.5, alpha = 0.8, aes(fill = Anno)) +
  scale_fill_manual(values = c("2002" = col_2002, "2024" = col_2024)) +
  labs(title = "Dispersione multivariata (Betadisper) - Bray-Curtis",
       subtitle = paste0("ANOVA p = ", round(bd_bc_test$`Pr(>F)`[1], 3)),
       x = "Anno", y = "Dispersione (Bray-Curtis)") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold"))

ggsave(file.path(fig_dir, "17_betadisper_braycurtis.png"), p_betadisper_bc, width = 6, height = 5, dpi = 300)

cat("Grafici Bray-Curtis salvati:\n")
cat("  15_nmds_braycurtis.png\n")
cat("  16_nmds_braycurtis_envfit.png\n")
cat("  17_betadisper_braycurtis.png\n")

# =============================================================================
# EXPORT TABELLA CONFRONTO
# =============================================================================

write.csv(comparison_df, file.path(tab_dir, "06_jaccard_vs_braycurtis.csv"), row.names = FALSE)
cat("\nTabella confronto salvata: 06_jaccard_vs_braycurtis.csv\n")

# =============================================================================
# 3) INDICATOR SPECIES ANALYSIS
# =============================================================================

cat("###############################################\n")
cat("# SEZIONE 3: INDICATOR SPECIES ANALYSIS\n")
cat("###############################################\n\n")

groups <- factor(plots_paired$Anno)

isa_ctrl <- how(
  within = Within(type = "free"),
  plots = Plots(strata = factor(plots_paired$PlotBase), type = "none"),
  nperm = 999
)

set.seed(123)
isa_res <- multipatt(comm_raw_p, groups, func = "IndVal.g", control = isa_ctrl, duleg = TRUE)

isa_summary <- isa_res$sign
isa_summary$Species <- rownames(isa_summary)
isa_summary <- isa_summary %>% arrange(p.value) %>% select(Species, everything())

isa_sig <- isa_summary %>% filter(p.value <= 0.05)

cat("Specie indicatrici significative (p <= 0.05):", nrow(isa_sig), "\n\n")

if(nrow(isa_sig) > 0) {
  isa_2002 <- isa_sig %>% filter(s.2002 == 1)
  isa_2024 <- isa_sig %>% filter(s.2024 == 1)
  
  cat("=== Indicatrici 2002 ===\n")
  if(nrow(isa_2002) > 0) {
    isa_2002 <- isa_2002 %>%
      left_join(traits %>% select(Species, L_num, Forma_biologica), by = "Species")
    print(isa_2002 %>% select(Species, stat, p.value, L_num, Forma_biologica))
  } else cat("Nessuna\n")
  
  cat("\n=== Indicatrici 2024 ===\n")
  if(nrow(isa_2024) > 0) {
    isa_2024 <- isa_2024 %>%
      left_join(traits %>% select(Species, L_num, Forma_biologica), by = "Species")
    print(isa_2024 %>% select(Species, stat, p.value, L_num, Forma_biologica))
  } else cat("Nessuna\n")
}

cat("\n=== SEZIONE 3 COMPLETATA ===\n\n")

# =============================================================================
# 4) ANALISI DENDROMETRICA
# =============================================================================

cat("###############################################\n")
cat("# SEZIONE 4: ANALISI DENDROMETRICA\n")
cat("###############################################\n\n")

# Mostra nomi originali per debug
cat("Colonne originali:\n")
print(names(dendro))

# Pulizia nomi colonne - R converte spazi e parentesi in punti
# Quindi "DBH (cm)" diventa "DBH..cm."
# Usiamo gsub per standardizzare
names(dendro) <- gsub("\\.+", "_", names(dendro))  # punti multipli -> singolo underscore
names(dendro) <- gsub("_+", "_", names(dendro))    # underscore multipli -> singolo
names(dendro) <- gsub("^_|_$", "", names(dendro))  # rimuovi underscore iniziali/finali
names(dendro) <- tolower(names(dendro))

cat("\nColonne dopo pulizia:\n")
print(names(dendro))

# Parse PlotID
dendro_parsed <- parse_plotid(dendro$plotid)
dendro$anno <- dendro_parsed$Anno
dendro$plot <- dendro_parsed$PlotBase

cat("\nPlot presenti nei dati dendrometrici:", paste(unique(dendro$plot), collapse = ", "), "\n")
cat("Anni:", paste(unique(dendro$anno), collapse = ", "), "\n")
cat("N. alberi totali:", nrow(dendro), "\n")

# Trova le colonne corrette usando pattern matching
dbh_col <- names(dendro)[grepl("^dbh_cm$|^dbh_cm_$", names(dendro))][1]
if(is.na(dbh_col)) dbh_col <- names(dendro)[grepl("dbh", names(dendro)) & grepl("cm", names(dendro))][1]

ba_col <- names(dendro)[grepl("ba_tree", names(dendro))][1]

size_col <- names(dendro)[grepl("plot_size", names(dendro))][1]

cat("\nColonne identificate:\n")
cat("  DBH:", dbh_col, "\n")
cat("  BA tree:", ba_col, "\n")
cat("  Plot size:", size_col, "\n\n")

# Estrai valori numerici (gestisce sia virgola che punto decimale)
dendro$dbh <- as.numeric(gsub(",", ".", dendro[[dbh_col]]))
dendro$ba_tree <- as.numeric(gsub(",", ".", dendro[[ba_col]]))
dendro$plot_size <- as.numeric(gsub(",", ".", dendro[[size_col]]))

# Fattore conversione per ettaro
dendro$conv_factor <- 10000 / dendro$plot_size

# NOTA: I dati dendrometrici sono disponibili solo per alcuni plot
# Le analisi seguenti si riferiscono solo a questi plot

cat("--- 4.1 Metriche stand ---\n")

stand_metrics <- dendro %>%
  group_by(plot, anno) %>%
  summarise(
    n_trees = n(),
    n_ha = n() * first(conv_factor),
    g_ha = sum(ba_tree, na.rm = TRUE) * first(conv_factor),
    dbh_mean = mean(dbh, na.rm = TRUE),
    dbh_sd = sd(dbh, na.rm = TRUE),
    dq = sqrt(mean(dbh^2, na.rm = TRUE)),
    dbh_max = max(dbh, na.rm = TRUE),
    .groups = "drop"
  )

cat("\nMetriche stand per plot e anno:\n")
print(stand_metrics)

# Calcolo variazioni per plot appaiati (se disponibili entrambi gli anni)
stand_wide <- stand_metrics %>%
  pivot_wider(
    names_from = anno, 
    values_from = c(n_trees, n_ha, g_ha, dbh_mean, dq, dbh_max),
    names_sep = "_"
  )

# Verifica quali plot hanno entrambi gli anni
plots_with_both <- stand_wide %>%
  filter(!is.na(n_ha_2002) & !is.na(n_ha_2024))

if(nrow(plots_with_both) > 0) {
  plots_with_both <- plots_with_both %>%
    mutate(
      delta_n_ha = n_ha_2024 - n_ha_2002,
      delta_g_ha = g_ha_2024 - g_ha_2002,
      delta_dq = dq_2024 - dq_2002
    )
  
  cat("\nVariazioni 2002-2024 per plot con dati completi:\n")
  print(plots_with_both %>% select(plot, delta_n_ha, delta_g_ha, delta_dq))
  
  cat("\nVariazioni medie:\n")
  cat("  Delta N/ha:", round(mean(plots_with_both$delta_n_ha, na.rm = TRUE), 1), "\n")
  cat("  Delta G (m²/ha):", round(mean(plots_with_both$delta_g_ha, na.rm = TRUE), 2), "m²/ha\n")
  cat("  Delta Dq (cm):", round(mean(plots_with_both$delta_dq, na.rm = TRUE), 2), "cm\n")
} else {
  cat("\nNOTA: Nessun plot con dati dendrometrici completi per entrambi gli anni\n")
}

cat("\n--- 4.2 Distribuzioni diametriche ---\n")

# Classi diametriche da 5 cm
dendro$dbh_class <- floor(dendro$dbh / 5) * 5 + 2.5

diam_dist <- dendro %>%
  group_by(anno, dbh_class) %>%
  summarise(
    n = n(),
    n_ha = sum(conv_factor),
    .groups = "drop"
  )

diam_wide <- diam_dist %>%
  select(anno, dbh_class, n_ha) %>%
  pivot_wider(names_from = anno, values_from = n_ha, values_fill = 0, names_prefix = "Y")

cat("\nDistribuzione N/ha per classe diametrica:\n")
print(diam_wide)

cat("\n--- 4.3 Composizione specifica del soprassuolo ---\n")

# Standardizza nomi specie
dendro$species_clean <- trimws(dendro$species)

species_comp <- dendro %>%
  group_by(anno, species_clean) %>%
  summarise(
    n = n(),
    n_ha = sum(conv_factor),
    g_ha = sum(ba_tree * conv_factor, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(anno, desc(g_ha))

cat("\nComposizione per specie (ordinate per G):\n")

cat("\n=== 2002 ===\n")
print(species_comp %>% filter(anno == 2002))

cat("\n=== 2024 ===\n")
print(species_comp %>% filter(anno == 2024))

# Riassunto per anno
species_summary <- species_comp %>%
  group_by(anno) %>%
  summarise(
    n_species = n_distinct(species_clean),
    tot_n_ha = sum(n_ha),
    tot_g_ha = sum(g_ha),
    .groups = "drop"
  )

cat("\nRiassunto totale per anno:\n")
print(species_summary)

cat("\n--- 4.4 Struttura verticale (se altezze disponibili) ---\n")

# Verifica colonna altezza
alt_col <- names(dendro)[grepl("altezza|height|h_m", names(dendro), ignore.case = TRUE)][1]

if(!is.na(alt_col)) {
  dendro$altezza <- as.numeric(gsub(",", ".", dendro[[alt_col]]))
  
  if(sum(!is.na(dendro$altezza) & dendro$altezza > 0) > 0) {
    alt_summary <- dendro %>%
      filter(!is.na(altezza) & altezza > 0) %>%
      group_by(anno) %>%
      summarise(
        n_with_height = n(),
        h_mean = mean(altezza, na.rm = TRUE),
        h_max = max(altezza, na.rm = TRUE),
        h_dom = mean(sort(altezza, decreasing = TRUE)[1:min(5, n())]),
        .groups = "drop"
      )
    
    cat("\nStatistiche altezze:\n")
    print(alt_summary)
  } else {
    cat("Dati di altezza non sufficienti\n")
  }
} else {
  cat("Colonna altezza non trovata\n")
}

cat("\n=== SEZIONE 4 COMPLETATA ===\n\n")


# =============================================================================
# 5) GRAFICI 
# =============================================================================

cat("###############################################\n")
cat("# SEZIONE 5: GRAFICI \n")
cat("###############################################\n\n")

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)


# Palette colori coerente
col_2002 <- "#E69F00"  # arancione
col_2024 <- "#0072B2"  # blu

# =============================================================================
# 5.1) RICCHEZZA SPECIFICA (con linee appaiate ed etichette)
# =============================================================================

cat("--- 5.1 Ricchezza specifica ---\n")

p_rich <- ggplot(plots_paired, aes(x = factor(Anno), y = Richness)) +
  geom_boxplot(aes(fill = factor(Anno)), alpha = 0.7, outlier.shape = NA, width = 0.5) +
  geom_line(aes(group = PlotBase), alpha = 0.4, color = "grey40") +
  geom_point(aes(fill = factor(Anno)), shape = 21, size = 2.5, alpha = 0.8) +
  geom_text(aes(label = PlotBase), size = 2, hjust = -0.3, vjust = 0.5, alpha = 0.7) +
  scale_fill_manual(values = c("2002" = col_2002, "2024" = col_2024)) +
  labs(title = "Ricchezza specifica 2002 vs 2024",
       subtitle = paste("Wilcoxon paired test: p =", round(wilcox_rich$p.value, 3)),
       x = "Anno", y = "Numero specie") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank())

ggsave(file.path(fig_dir, "01_richness_paired.png"), p_rich, width = 7, height = 6, dpi = 300)

# =============================================================================
# 5.2) RICCHEZZA STANDARDIZZATA A 100 m²
# =============================================================================

cat("--- 5.2 Ricchezza standardizzata ---\n")

p_rich100 <- ggplot(plots_paired, aes(x = factor(Anno), y = Richness_100)) +
  geom_boxplot(aes(fill = factor(Anno)), alpha = 0.7, outlier.shape = NA, width = 0.5) +
  geom_point(aes(fill = factor(Anno)), shape = 21, size = 2.5, alpha = 0.8) +
  scale_fill_manual(values = c("2002" = col_2002, "2024" = col_2024)) +
  labs(title = "Ricchezza standardizzata a 100 m²",
       subtitle = "NOTA: aree 2002 non standard (50-120 m²)",
       x = "Anno", y = "Ricchezza standardizzata a 100 m²") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold"))

ggsave(file.path(fig_dir, "02_richness_100m2.png"), p_rich100, width = 6, height = 5, dpi = 300)

# =============================================================================
# 5.3) SPECIE HABITAT APERTO vs CHIUSO
# =============================================================================

cat("--- 5.3 Habitat aperto vs chiuso ---\n")

habitat_long <- plots_paired %>%
  select(PlotBase, Anno, N_aperto, N_chiuso) %>%
  pivot_longer(cols = c(N_aperto, N_chiuso), names_to = "Habitat", values_to = "N") %>%
  mutate(Habitat = ifelse(Habitat == "N_aperto", "Aperto (L>=7)", "Chiuso (L<=4)"))

p_habitat <- ggplot(habitat_long, aes(x = factor(Anno), y = N, fill = Habitat)) +
  geom_boxplot(alpha = 0.7, position = position_dodge(0.8), width = 0.6) +
  scale_fill_manual(values = c("Aperto (L>=7)" = col_2002, "Chiuso (L<=4)" = "#009E73")) +
  labs(title = "Specie per tipo di habitat",
       x = "Anno", y = "Numero specie") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "right")

ggsave(file.path(fig_dir, "03_habitat_species.png"), p_habitat, width = 7, height = 5, dpi = 300)

# =============================================================================
# 5.4) SPETTRO FORME BIOLOGICHE
# =============================================================================

cat("--- 5.4 Forme biologiche ---\n")

fb_long <- plots_paired %>%
  select(PlotBase, Anno, starts_with("FB_")) %>%
  pivot_longer(cols = starts_with("FB_"), names_to = "Forma", values_to = "Perc") %>%
  mutate(Forma = gsub("FB_", "", Forma),
         Forma = factor(Forma, levels = c("P", "Np", "Ch", "H", "G", "T")))

p_fb <- ggplot(fb_long, aes(x = Forma, y = Perc, fill = factor(Anno))) +
  geom_boxplot(alpha = 0.7, position = position_dodge(0.8), width = 0.6) +
  scale_fill_manual(values = c("2002" = col_2002, "2024" = col_2024), name = "Anno") +
  labs(title = "Spettro forme biologiche",
       x = "Forma biologica", y = "% copertura") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(fig_dir, "04_life_forms.png"), p_fb, width = 9, height = 5, dpi = 300)

# =============================================================================
# 5.5) ELLENBERG CWM - Boxplot con linee appaiate
# =============================================================================

cat("--- 5.5 Ellenberg CWM ---\n")

# Prepara dati wide per calcolo delta
ell_L_wide <- plots_paired %>%
  select(PlotBase, Anno, Ell_L) %>%
  pivot_wider(names_from = Anno, values_from = Ell_L, names_prefix = "Y") %>%
  mutate(delta = Y2024 - Y2002)

ell_U_wide <- plots_paired %>%
  select(PlotBase, Anno, Ell_U) %>%
  pivot_wider(names_from = Anno, values_from = Ell_U, names_prefix = "Y") %>%
  mutate(delta = Y2024 - Y2002)

ell_N_wide <- plots_paired %>%
  select(PlotBase, Anno, Ell_N) %>%
  pivot_wider(names_from = Anno, values_from = Ell_N, names_prefix = "Y") %>%
  mutate(delta = Y2024 - Y2002)

# Boxplot L
p_ell_L <- ggplot(plots_paired, aes(x = factor(Anno), y = Ell_L)) +
  geom_boxplot(aes(fill = factor(Anno)), alpha = 0.7, outlier.shape = NA, width = 0.5) +
  geom_line(aes(group = PlotBase), alpha = 0.3, color = "grey40") +
  geom_point(aes(fill = factor(Anno)), shape = 21, size = 2, alpha = 0.8) +
  scale_fill_manual(values = c("2002" = col_2002, "2024" = col_2024)) +
  labs(x = "Anno", y = "Ellenberg L (weighted mean)") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none")

# Boxplot U
p_ell_U <- ggplot(plots_paired, aes(x = factor(Anno), y = Ell_U)) +
  geom_boxplot(aes(fill = factor(Anno)), alpha = 0.7, outlier.shape = NA, width = 0.5) +
  geom_line(aes(group = PlotBase), alpha = 0.3, color = "grey40") +
  geom_point(aes(fill = factor(Anno)), shape = 21, size = 2, alpha = 0.8) +
  scale_fill_manual(values = c("2002" = col_2002, "2024" = col_2024)) +
  labs(x = "Anno", y = "Ellenberg U (weighted mean)") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none")

# Boxplot N
p_ell_N <- ggplot(plots_paired, aes(x = factor(Anno), y = Ell_N)) +
  geom_boxplot(aes(fill = factor(Anno)), alpha = 0.7, outlier.shape = NA, width = 0.5) +
  geom_line(aes(group = PlotBase), alpha = 0.3, color = "grey40") +
  geom_point(aes(fill = factor(Anno)), shape = 21, size = 2, alpha = 0.8) +
  scale_fill_manual(values = c("2002" = col_2002, "2024" = col_2024)) +
  labs(x = "Anno", y = "Ellenberg N (weighted mean)") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none")

# Combina boxplot
p_ell_box <- p_ell_L + p_ell_N + p_ell_U +
  plot_annotation(title = "Ellenberg CWM",
                  theme = theme(plot.title = element_text(face = "bold", size = 14)))

ggsave(file.path(fig_dir, "05_ellenberg_cwm_boxplot.png"), p_ell_box, width = 12, height = 4, dpi = 300)

# Istogrammi Delta
p_delta_L <- ggplot(ell_L_wide, aes(x = delta)) +
  geom_histogram(bins = 10, fill = "grey60", color = "black", alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Distribuzione ΔL", x = "ΔL (2024 - 2002)", y = "Frequency") +
  theme_minimal(base_size = 11)

p_delta_U <- ggplot(ell_U_wide, aes(x = delta)) +
  geom_histogram(bins = 10, fill = "grey60", color = "black", alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Distribuzione ΔU", x = "ΔU (2024 - 2002)", y = "Frequency") +
  theme_minimal(base_size = 11)

p_delta_N <- ggplot(ell_N_wide, aes(x = delta)) +
  geom_histogram(bins = 10, fill = "grey60", color = "black", alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Distribuzione ΔN", x = "ΔN (2024 - 2002)", y = "Frequency") +
  theme_minimal(base_size = 11)

p_ell_delta <- p_delta_L + p_delta_U + p_delta_N

ggsave(file.path(fig_dir, "06_ellenberg_delta_histograms.png"), p_ell_delta, width = 12, height = 4, dpi = 300)

# =============================================================================
# 5.6) NMDS CON TRAIETTORIE
# =============================================================================

cat("--- 5.6 NMDS con traiettorie ---\n")

# Prepara dati wide per frecce
nmds_wide <- nmds_scores %>%
  select(PlotBase, Anno, NMDS1, NMDS2) %>%
  pivot_wider(names_from = Anno, values_from = c(NMDS1, NMDS2))

p_nmds <- ggplot() +
  # Frecce traiettorie
  geom_segment(data = nmds_wide,
               aes(x = NMDS1_2002, y = NMDS2_2002, 
                   xend = NMDS1_2024, yend = NMDS2_2024),
               arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
               alpha = 0.5, color = "grey50") +
  # Punti 2002
  geom_point(data = nmds_scores %>% filter(Anno == 2002),
             aes(x = NMDS1, y = NMDS2), 
             shape = 21, fill = col_2002, size = 3.5, stroke = 0.5) +
  # Punti 2024
  geom_point(data = nmds_scores %>% filter(Anno == 2024),
             aes(x = NMDS1, y = NMDS2), 
             shape = 22, fill = col_2024, size = 3.5, stroke = 0.5) +
  # Etichette
  geom_text(data = nmds_scores,
            aes(x = NMDS1, y = NMDS2, label = PlotBase),
            size = 2.5, vjust = -0.8, alpha = 0.8) +
  labs(title = paste0("NMDS con traiettorie 2002 → 2024"),
       subtitle = paste0("Stress = ", round(nmds$stress, 3)),
       x = "NMDS1", y = "NMDS2") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        aspect.ratio = 1) +
  # Legenda manuale
  annotate("point", x = Inf, y = Inf, shape = 21, fill = col_2002, size = 3, vjust = 2, hjust = 1.5) +
  annotate("text", x = Inf, y = Inf, label = "2002", vjust = 2.2, hjust = 0.5, size = 3) +
  annotate("point", x = Inf, y = Inf, shape = 22, fill = col_2024, size = 3, vjust = 4, hjust = 1.5) +
  annotate("text", x = Inf, y = Inf, label = "2024", vjust = 4.2, hjust = 0.5, size = 3)

ggsave(file.path(fig_dir, "07_nmds_trajectories.png"), p_nmds, width = 8, height = 8, dpi = 300)

# =============================================================================
# 5.7) BETADISPER
# =============================================================================

cat("--- 5.7 Betadisper ---\n")

bd_df <- data.frame(
  Anno = factor(plots_paired$Anno),
  Distance = bd$distances,
  PlotBase = plots_paired$PlotBase
)

p_betadisper <- ggplot(bd_df, aes(x = Anno, y = Distance)) +
  geom_boxplot(aes(fill = Anno), alpha = 0.7, outlier.shape = NA, width = 0.5) +
  geom_point(shape = 21, size = 2.5, alpha = 0.8, aes(fill = Anno)) +
  scale_fill_manual(values = c("2002" = col_2002, "2024" = col_2024)) +
  labs(title = "Dispersione multivarata (Betadisper)",
       subtitle = paste0("ANOVA p = ", round(bd_test$`Pr(>F)`[1], 3)),
       x = "Anno", y = "Dispersione (Jaccard)") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold"))

ggsave(file.path(fig_dir, "08_betadisper.png"), p_betadisper, width = 6, height = 5, dpi = 300)

# =============================================================================
# 5.8) NMDS CON ENVFIT
# =============================================================================

cat("--- 5.8 NMDS con Envfit ---\n")

# Estrai vettori envfit
ef_vectors <- as.data.frame(scores(ef, display = "vectors"))
ef_vectors$var <- rownames(ef_vectors)

# Scala vettori per visualizzazione
arrow_scale <- 0.8

p_nmds_envfit <- ggplot() +
  # Punti
  geom_point(data = nmds_scores,
             aes(x = NMDS1, y = NMDS2, shape = factor(Anno), fill = factor(Anno)),
             size = 3, stroke = 0.5) +
  scale_shape_manual(values = c("2002" = 21, "2024" = 22)) +
  scale_fill_manual(values = c("2002" = col_2002, "2024" = col_2024)) +
  # Vettori envfit
  geom_segment(data = ef_vectors,
               aes(x = 0, y = 0, xend = NMDS1 * arrow_scale, yend = NMDS2 * arrow_scale),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "black", linewidth = 1) +
  geom_text(data = ef_vectors,
            aes(x = NMDS1 * arrow_scale * 1.1, y = NMDS2 * arrow_scale * 1.1, label = var),
            size = 4, fontface = "bold") +
  labs(title = "NMDS con gradienti Ellenberg",
       subtitle = paste0("Envfit: L r² = ", round(ef$vectors$r["Ell_L"], 2),
                         ", U r² = ", round(ef$vectors$r["Ell_U"], 2),
                         ", N r² = ", round(ef$vectors$r["Ell_N"], 2)),
       x = "NMDS1", y = "NMDS2", shape = "Anno", fill = "Anno") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        aspect.ratio = 1)

ggsave(file.path(fig_dir, "09_nmds_envfit.png"), p_nmds_envfit, width = 8, height = 7, dpi = 300)

# =============================================================================
# 5.9) DENDROMETRIA - Traiettorie per plot
# =============================================================================

cat("--- 5.9 Dendrometria traiettorie ---\n")

# Verifica che stand_metrics esista e abbia dati per entrambi gli anni
if(exists("stand_metrics") && nrow(stand_metrics) > 0) {
  
  # Filtra solo plot con entrambi gli anni
  plots_dendro_paired <- stand_metrics %>%
    group_by(plot) %>%
    filter(n() == 2) %>%
    ungroup()
  
  if(nrow(plots_dendro_paired) > 0) {
    
    # Traiettoria Dq
    p_dq <- ggplot(plots_dendro_paired, aes(x = anno, y = dq, color = factor(plot), group = plot)) +
      geom_line(linewidth = 1.2) +
      geom_point(size = 3) +
      scale_color_brewer(palette = "Set2", name = "Plot") +
      labs(title = "Dg (cm)", x = "Anno", y = "Dg (cm)") +
      theme_minimal(base_size = 11) +
      theme(legend.position = "top")
    
    # Traiettoria G
    p_g <- ggplot(plots_dendro_paired, aes(x = anno, y = g_ha, color = factor(plot), group = plot)) +
      geom_line(linewidth = 1.2) +
      geom_point(size = 3) +
      scale_color_brewer(palette = "Set2", name = "Plot") +
      labs(title = "G (m²/ha)", x = "Anno", y = "G (m²/ha)") +
      theme_minimal(base_size = 11) +
      theme(legend.position = "top")
    
    # Traiettoria N
    p_n <- ggplot(plots_dendro_paired, aes(x = anno, y = n_ha, color = factor(plot), group = plot)) +
      geom_line(linewidth = 1.2) +
      geom_point(size = 3) +
      scale_color_brewer(palette = "Set2", name = "Plot") +
      labs(title = "N (stems/ha)", x = "Anno", y = "N (stems/ha)") +
      theme_minimal(base_size = 11) +
      theme(legend.position = "top")
    
    p_dendro_traj <- p_dq + p_g + p_n +
      plot_layout(guides = "collect") +
      plot_annotation(title = "Traiettorie dendrometriche 2002→2024 per plot",
                      theme = theme(plot.title = element_text(face = "bold.italic", size = 14)))
    
    ggsave(file.path(fig_dir, "10_dendro_trajectories.png"), p_dendro_traj, width = 14, height = 5, dpi = 300)
    
  } else {
    cat("  Nessun plot con dati dendrometrici appaiati\n")
  }
} else {
  cat("  Dati dendrometrici non disponibili\n")
}

# =============================================================================
# 5.10) DENDROMETRIA - Area basimetrica per classe e anno
# =============================================================================

cat("--- 5.10 BA per classe diametrica ---\n")

if(exists("dendro") && nrow(dendro) > 0) {
  
  # Calcola BA per classe
  ba_by_class <- dendro %>%
    group_by(anno, dbh_class) %>%
    summarise(ba_ha = sum(ba_tree * conv_factor, na.rm = TRUE), .groups = "drop")
  
  p_ba_class <- ggplot(ba_by_class, aes(x = dbh_class, y = ba_ha, color = factor(anno))) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3) +
    scale_color_manual(values = c("2002" = "firebrick", "2024" = "steelblue"), name = "Anno") +
    labs(title = "Area basimetrica per classe diametrica (aggregato tra plot)",
         x = "Classe DBH (cm)", y = "BA (m²/ha)") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "right")
  
  ggsave(file.path(fig_dir, "11_ba_by_dbh_class.png"), p_ba_class, width = 9, height = 5, dpi = 300)
}

# =============================================================================
# 5.11) DENDROMETRIA - BA per classe e specie (stacked)
# =============================================================================

cat("--- 5.11 BA per classe e specie ---\n")

if(exists("dendro") && nrow(dendro) > 0) {
  
  ba_by_class_sp <- dendro %>%
    group_by(anno, dbh_class, species_clean) %>%
    summarise(ba_ha = sum(ba_tree * conv_factor, na.rm = TRUE), .groups = "drop")
  
  p_ba_sp <- ggplot(ba_by_class_sp, aes(x = dbh_class, y = ba_ha, fill = species_clean)) +
    geom_bar(stat = "identity", position = "stack", width = 4) +
    facet_wrap(~ anno, ncol = 2) +
    scale_fill_brewer(palette = "Set3", name = "Specie") +
    labs(title = "Area basimetrica per classe diametrica e specie (aggregato)",
         x = "Classe DBH (cm)", y = "BA (m²/ha)") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "right")
  
  ggsave(file.path(fig_dir, "12_ba_by_species_class.png"), p_ba_sp, width = 14, height = 6, dpi = 300)
}

# =============================================================================
# 5.12) DENDROMETRIA - Densità per classe e specie (stacked)
# =============================================================================

cat("--- 5.12 Densità per classe e specie ---\n")

if(exists("dendro") && nrow(dendro) > 0) {
  
  n_by_class_sp <- dendro %>%
    group_by(anno, dbh_class, species_clean) %>%
    summarise(n_ha = sum(conv_factor), .groups = "drop")
  
  p_n_sp <- ggplot(n_by_class_sp, aes(x = dbh_class, y = n_ha, fill = species_clean)) +
    geom_bar(stat = "identity", position = "stack", width = 4) +
    facet_wrap(~ anno, ncol = 2) +
    scale_fill_brewer(palette = "Set3", name = "Specie") +
    labs(title = "Densità per classe diametrica e specie (aggregato)",
         x = "Classe DBH (cm)", y = "Densità (stems/ha)") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "right")
  
  ggsave(file.path(fig_dir, "13_density_by_species_class.png"), p_n_sp, width = 14, height = 6, dpi = 300)
}

# =============================================================================
# 5.13) NMDS COMBINATO CON PERMANOVA INFO
# =============================================================================

cat("--- 5.13 NMDS con info PERMANOVA ---\n")

p_nmds_perm <- ggplot() +
  geom_segment(data = nmds_wide,
               aes(x = NMDS1_2002, y = NMDS2_2002, 
                   xend = NMDS1_2024, yend = NMDS2_2024),
               alpha = 0.4, color = "grey50") +
  geom_point(data = nmds_scores %>% filter(Anno == 2002),
             aes(x = NMDS1, y = NMDS2), 
             shape = 21, fill = col_2002, size = 3.5, stroke = 0.5) +
  geom_point(data = nmds_scores %>% filter(Anno == 2024),
             aes(x = NMDS1, y = NMDS2), 
             shape = 22, fill = col_2024, size = 3.5, stroke = 0.5) +
  geom_text(data = nmds_scores,
            aes(x = NMDS1, y = NMDS2, label = PlotBase),
            size = 2.2, vjust = -0.7, alpha = 0.7) +
  labs(title = paste0("NMDS (Jaccard) - Stress = ", round(nmds$stress, 3)),
       subtitle = paste0("PERMANOVA R² = ", round(adonis_res$R2[1], 2), 
                         ", p = ", round(adonis_res$`Pr(>F)`[1], 3)),
       x = "NMDS1", y = "NMDS2") +
  annotate("point", x = max(nmds_scores$NMDS1) * 0.9, y = max(nmds_scores$NMDS2) * 0.95, 
           shape = 21, fill = col_2002, size = 4) +
  annotate("text", x = max(nmds_scores$NMDS1) * 0.95, y = max(nmds_scores$NMDS2) * 0.95, 
           label = "2002", hjust = 0, size = 3.5) +
  annotate("point", x = max(nmds_scores$NMDS1) * 0.9, y = max(nmds_scores$NMDS2) * 0.85, 
           shape = 22, fill = col_2024, size = 4) +
  annotate("text", x = max(nmds_scores$NMDS1) * 0.95, y = max(nmds_scores$NMDS2) * 0.85, 
           label = "2024", hjust = 0, size = 3.5) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        aspect.ratio = 1)

ggsave(file.path(fig_dir, "14_nmds_permanova.png"), p_nmds_perm, width = 8, height = 8, dpi = 300)

cat("\n=== GRAFICI SALVATI IN:", fig_dir, "===\n")
cat("File generati:\n")
cat("  01_richness_paired.png\n")
cat("  02_richness_100m2.png\n")
cat("  03_habitat_species.png\n")
cat("  04_life_forms.png\n")
cat("  05_ellenberg_cwm_boxplot.png\n")
cat("  06_ellenberg_delta_histograms.png\n")
cat("  07_nmds_trajectories.png\n")
cat("  08_betadisper.png\n")
cat("  09_nmds_envfit.png\n")
cat("  10_dendro_trajectories.png\n")
cat("  11_ba_by_dbh_class.png\n")
cat("  12_ba_by_species_class.png\n")
cat("  13_density_by_species_class.png\n")
cat("  14_nmds_permanova.png\n")

cat("\n=== SEZIONE 5 COMPLETATA ===\n\n")

# =============================================================================
# 6) EXPORT TABELLE
# =============================================================================

cat("\n###############################################\n")
cat("# SEZIONE 6: EXPORT TABELLE\n")
cat("###############################################\n\n")


summary_table <- plots_paired %>%
  group_by(Anno) %>%
  summarise(
    n_plots = n(), S_mean = mean(Richness), S_sd = sd(Richness),
    N_aperto_mean = mean(N_aperto), N_chiuso_mean = mean(N_chiuso),
    N_AFS_mean = mean(N_AFS),
    Ell_L_mean = mean(Ell_L, na.rm = TRUE),
    Ell_U_mean = mean(Ell_U, na.rm = TRUE),
    Ell_N_mean = mean(Ell_N, na.rm = TRUE),
    .groups = "drop"
  )
write.csv(summary_table, file.path(tab_dir, "01_summary_by_year.csv"), row.names = FALSE)

delta_table <- data.frame(
  PlotBase = rich_wide$PlotBase,
  delta_S = rich_wide$delta_S,
  delta_aperto = aperto_wide$delta,
  delta_chiuso = chiuso_wide$delta
)
write.csv(delta_table, file.path(tab_dir, "02_delta_by_plot.csv"), row.names = FALSE)

write.csv(isa_summary, file.path(tab_dir, "03_indicator_species.csv"), row.names = FALSE)

permanova_df <- data.frame(
  Term = rownames(adonis_res),
  Df = adonis_res$Df, SumOfSqs = adonis_res$SumOfSqs,
  R2 = adonis_res$R2, F_val = adonis_res$F, p_value = adonis_res$`Pr(>F)`
)
write.csv(permanova_df, file.path(tab_dir, "04_permanova.csv"), row.names = FALSE)

write.csv(fb_summary, file.path(tab_dir, "05_life_forms_spectrum.csv"), row.names = FALSE)

cat("Tabelle salvate in:", tab_dir, "\n")

cat("\n###############################################\n")
cat("# ANALISI COMPLETATE\n")
cat("###############################################\n")