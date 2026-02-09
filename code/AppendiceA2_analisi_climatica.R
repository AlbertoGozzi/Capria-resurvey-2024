################################################################################
#                    ANALISI CLIMATICA - CLIMOGRAMMI DI WALTER
#                    APPENDICE A.2 - Script R analisi climatica
#                    Stazioni di Corniolo (590m) e Campigna (1068m)
#                    Confronto periodo Salvatori vs periodo attuale
#
#  Autore: Alberto Gozzi
#  Affiliazione: Università di Bologna - Distal
#  Email: alberto.gozzi@tudio.unibo.it
#  Data: Febbraio 2025
#  Versione: 1.0
#
################################################################################

# =============================================================================
# PANORAMICA DELLO SCRIPT
# =============================================================================
#
# OBIETTIVO:
#   Caratterizzare il clima dell'area di studio attraverso climogrammi di
#   Walter per due stazioni altitudinali (Corniolo 590m, Campigna 1068m)
#   confrontando il periodo dello studio Salvatori (1990-2010) con il
#   periodo attuale (2010-2020).
#
# METODO:
#   - Climogramma di Walter (diagramma ombrotermico)
#   - Scala: 1°C = 2mm precipitazione
#   - Compressione per P > 100mm (scala 1:10)
#   - Identificazione periodi secchi (P < 2T)
#
# INPUT RICHIESTI:
#   Temperature:
#     - media_campigna.xlsx (medie mensili 2000-2010, 2010-2020)
#     - media_corniolo.xlsx (medie mensili 2000-2010, 2010-2020)
#   
#   Precipitazioni:
#     - precip_media_campigna_90_02.csv (giornaliere 1990-2002)
#     - precip_media_corniolo_90_02.csv (giornaliere 1990-2002)
#     - Precip_campigna_01_25.csv (giornaliere 2001-2025)
#     - Precip_Corniolo_01_25.csv (giornaliere 2001-2025)
#
# OUTPUT GENERATI:
#   - Climogramma_Salvatori.png (periodo 1990-2010)
#   - Climogramma_Attuale.png (periodo 2010-2020)
#   - Statistiche climatiche in console
#   - Tabella confronto periodi
#
# TEMPO ESECUZIONE: ~30 secondi
#
# =============================================================================

################################################################################
# INIZIO CODICE
################################################################################

################################################################################
#                    ANALISI CLIMATICA - CLIMOGRAMMI DI WALTER                  
#                    Stazioni di Corniolo (590m) e Campigna (1068m)            
#                    Confronto periodo Salvatori vs periodo attuale            
################################################################################

# ==============================================================================
# LIBRERIE NECESSARIE
# ==============================================================================
# install.packages(c("readxl", "dplyr", "ggplot2"))
library(readxl)
library(dplyr)
library(ggplot2)

# ==============================================================================
# 1. CARICAMENTO E PREPARAZIONE DATI
# ==============================================================================

# --- Temperature (medie decennali già calcolate) ---
# File: media_campigna.xlsx e media_corniolo.xlsx
# Struttura: Mese | Media 2000-2010 | Media 2010-2020

temp_campigna <- read_excel("media_campigna.xlsx", skip = 1, col_names = FALSE)
colnames(temp_campigna) <- c("mese", "T_2000_2010", "T_2010_2020")

temp_corniolo <- read_excel("media_corniolo.xlsx", skip = 1, col_names = FALSE)
colnames(temp_corniolo) <- c("mese", "T_2000_2010", "T_2010_2020")

# Estrazione vettori temperature
temp_camp_2000_2010 <- as.numeric(temp_campigna$T_2000_2010)
temp_camp_2010_2020 <- as.numeric(temp_campigna$T_2010_2020)
temp_corn_2000_2010 <- as.numeric(temp_corniolo$T_2000_2010)
temp_corn_2010_2020 <- as.numeric(temp_corniolo$T_2010_2020)

# --- Precipitazioni periodo Salvatori (1990-2002) ---
# File: precip_media_campigna_90_02.csv e precip_media_corniolo_90_02.csv

load_precip_90_02 <- function(filepath) {
  df <- read.csv2(filepath, skip = 3, header = FALSE, 
                  fileEncoding = "latin1", stringsAsFactors = FALSE)
  colnames(df) <- c("inizio", "fine", "precip")
  df$data <- as.Date(substr(df$inizio, 1, 10))
  df$anno <- as.numeric(format(df$data, "%Y"))
  df$mese <- as.numeric(format(df$data, "%m"))
  df$precip <- as.numeric(gsub(",", ".", df$precip))
  
  # Filtro 1990-2002
  df <- df[df$anno >= 1990 & df$anno <= 2002, ]
  
  # Somme mensili per anno
  somme <- df %>%
    group_by(anno, mese) %>%
    summarise(tot = sum(precip, na.rm = TRUE), .groups = "drop")
  
  # Media per mese
  medie <- somme %>%
    group_by(mese) %>%
    summarise(media = mean(tot, na.rm = TRUE), .groups = "drop")
  
  return(medie$media)
}

precip_camp_90_02 <- load_precip_90_02("precip_media_campigna_90_02.csv")
precip_corn_90_02 <- load_precip_90_02("precip_media_corniolo_90_02.csv")

# --- Precipitazioni periodo attuale (2001-2019/2017) ---
# File: Precip_campigna_01_25.csv e Precip_Corniolo_01_25.csv

load_precip_attuale <- function(filepath, anno_start, anno_end) {
  df <- read.csv2(filepath, skip = 2, header = FALSE,
                  fileEncoding = "latin1", stringsAsFactors = FALSE)
  colnames(df) <- c("inizio", "fine", "precip")
  df$data <- as.Date(substr(df$inizio, 1, 10))
  df$anno <- as.numeric(format(df$data, "%Y"))
  df$mese <- as.numeric(format(df$data, "%m"))
  df$precip <- as.numeric(gsub(",", ".", df$precip))
  
  df <- df[df$anno >= anno_start & df$anno <= anno_end, ]
  
  somme <- df %>%
    group_by(anno, mese) %>%
    summarise(tot = sum(precip, na.rm = TRUE), .groups = "drop")
  
  medie <- somme %>%
    group_by(mese) %>%
    summarise(media = mean(tot, na.rm = TRUE), .groups = "drop")
  
  return(medie$media)
}

precip_camp_attuale <- load_precip_attuale("Precip_campigna_01_25.csv", 2001, 2019)
precip_corn_attuale <- load_precip_attuale("Precip_Corniolo_01_25.csv", 2001, 2017)

# ==============================================================================
# 2. STATISTICHE CLIMATICHE
# ==============================================================================

calc_stats <- function(temp, precip, nome) {
  cat("\n", nome, "\n")
  cat(paste(rep("-", 50), collapse = ""), "\n")
  cat("Temperatura media annua:", round(mean(temp), 1), "°C\n")
  cat("T min (mese più freddo):", round(min(temp), 1), "°C\n")
  cat("T max (mese più caldo):", round(max(temp), 1), "°C\n")
  cat("Escursione termica:", round(max(temp) - min(temp), 1), "°C\n")
  cat("Precipitazioni annue:", round(sum(precip)), "mm\n")
  
  # Classificazione continentalità
  esc <- max(temp) - min(temp)
  if (esc > 20) {
    tipo <- "CONTINENTALE"
  } else if (esc > 15) {
    tipo <- "SUBCONTINENTALE"
  } else if (esc > 10) {
    tipo <- "SUBOCEANICO"
  } else {
    tipo <- "OCEANICO"
  }
  cat("Tipo climatico:", tipo, "\n")
  
  # Periodi secchi (P < 2T)
  mesi_secchi <- which(precip < 2 * temp)
  if (length(mesi_secchi) > 0) {
    mesi_nomi <- c("Gen", "Feb", "Mar", "Apr", "Mag", "Giu", 
                   "Lug", "Ago", "Set", "Ott", "Nov", "Dic")
    cat("Mesi secchi (P < 2T):", paste(mesi_nomi[mesi_secchi], collapse = ", "), "\n")
  } else {
    cat("Mesi secchi: nessuno\n")
  }
}

cat("\n========== PERIODO SALVATORI (T: 2000-2010, P: 1990-2002) ==========\n")
calc_stats(temp_camp_2000_2010, precip_camp_90_02, "CAMPIGNA (1068 m)")
calc_stats(temp_corn_2000_2010, precip_corn_90_02, "CORNIOLO (590 m)")

cat("\n========== PERIODO ATTUALE (T: 2010-2020, P: 2001-2019) ==========\n")
calc_stats(temp_camp_2010_2020, precip_camp_attuale, "CAMPIGNA (1068 m)")
calc_stats(temp_corn_2010_2020, precip_corn_attuale, "CORNIOLO (590 m)")

# ==============================================================================
# 3. FUNZIONE PER CLIMOGRAMMA DI WALTER
# ==============================================================================

# Scala Walter: per P > 100mm la scala si comprime a 1/10
scala_walter <- function(p) {
  ifelse(p <= 100, p, 100 + (p - 100) / 10)
}

# Funzione per creare climogramma base R
create_climogram <- function(temp, precip, titolo, quota) {
  
  mesi <- c("G", "F", "M", "A", "M", "G", "L", "A", "S", "O", "N", "D")
  x <- 1:12
  
  # Scala Walter
  temp_scaled <- temp * 2  # 1°C = 2mm
  precip_scaled <- scala_walter(precip)
  
  # Limiti Y
  ymax <- max(c(temp_scaled, precip_scaled)) * 1.1
  
  # Plot
  plot(x, precip_scaled, type = "n", 
       xlim = c(0.5, 12.5), ylim = c(0, ymax),
       xlab = "", ylab = "mm / °C×2",
       main = paste0(titolo, " m ", quota),
       xaxt = "n", las = 1)
  axis(1, at = x, labels = mesi)
  
  # Interpolazione per riempimento
  x_fine <- seq(1, 12, length.out = 200)
  temp_fine <- approx(x, temp_scaled, x_fine)$y
  precip_fine <- approx(x, precip_scaled, x_fine)$y
  
  # Area umida (P > 2T) - tratteggio
  for (i in 1:(length(x_fine)-1)) {
    if (precip_fine[i] >= temp_fine[i]) {
      polygon(c(x_fine[i], x_fine[i+1], x_fine[i+1], x_fine[i]),
              c(temp_fine[i], temp_fine[i+1], precip_fine[i+1], precip_fine[i]),
              col = "gray90", border = NA, density = 30, angle = 90)
    }
  }
  
  # Area secca (P < 2T) - puntini
  for (i in 1:(length(x_fine)-1)) {
    if (temp_fine[i] > precip_fine[i]) {
      polygon(c(x_fine[i], x_fine[i+1], x_fine[i+1], x_fine[i]),
              c(precip_fine[i], precip_fine[i+1], temp_fine[i+1], temp_fine[i]),
              col = "gray80", border = NA, density = 10)
    }
  }
  
  # Linee
  lines(x, temp_scaled, lwd = 2, col = "black")
  lines(x, precip_scaled, lwd = 2, col = "black")
  
  # Asse Y con valori reali precipitazione
  if (max(precip) > 100) {
    at_vals <- c(0, 50, 100)
    lab_vals <- c("0", "50", "100")
    for (v in seq(150, 500, 50)) {
      pos <- scala_walter(v)
      if (pos <= ymax) {
        at_vals <- c(at_vals, pos)
        lab_vals <- c(lab_vals, as.character(v))
      }
    }
    axis(2, at = at_vals, labels = lab_vals, las = 1)
  }
}

# ==============================================================================
# 4. GENERAZIONE GRAFICI
# ==============================================================================

# --- Periodo Salvatori ---
png("Climogramma_Salvatori.png", width = 1200, height = 600, res = 120)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 2))

create_climogram(temp_corn_2000_2010, precip_corn_90_02, "CORNIOLO", 590)
create_climogram(temp_camp_2000_2010, precip_camp_90_02, "CAMPIGNA", 1068)

mtext("Diagramma ombrotermico - Periodo 1990-2010", 
      side = 1, outer = TRUE, line = -1, cex = 1.2)
dev.off()

# --- Periodo Attuale ---
png("Climogramma_Attuale.png", width = 1200, height = 600, res = 120)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 2))

create_climogram(temp_corn_2010_2020, precip_corn_attuale, "CORNIOLO", 590)
create_climogram(temp_camp_2010_2020, precip_camp_attuale, "CAMPIGNA", 1068)

mtext("Diagramma ombrotermico - Periodo 2010-2020", 
      side = 1, outer = TRUE, line = -1, cex = 1.2)
dev.off()

# ==============================================================================
# 5. TABELLA CONFRONTO
# ==============================================================================

cat("\n\n========== TABELLA CONFRONTO PERIODI ==========\n\n")

df_confronto <- data.frame(
  Stazione = c("Corniolo", "Corniolo", "Campigna", "Campigna"),
  Periodo = c("1990-2010", "2010-2020", "1990-2010", "2010-2020"),
  T_media = c(mean(temp_corn_2000_2010), mean(temp_corn_2010_2020),
              mean(temp_camp_2000_2010), mean(temp_camp_2010_2020)),
  T_gen = c(temp_corn_2000_2010[1], temp_corn_2010_2020[1],
            temp_camp_2000_2010[1], temp_camp_2010_2020[1]),
  T_lug = c(temp_corn_2000_2010[7], temp_corn_2010_2020[7],
            temp_camp_2000_2010[7], temp_camp_2010_2020[7]),
  Escursione = c(max(temp_corn_2000_2010) - min(temp_corn_2000_2010),
                 max(temp_corn_2010_2020) - min(temp_corn_2010_2020),
                 max(temp_camp_2000_2010) - min(temp_camp_2000_2010),
                 max(temp_camp_2010_2020) - min(temp_camp_2010_2020)),
  P_annua = c(sum(precip_corn_90_02), sum(precip_corn_attuale),
              sum(precip_camp_90_02), sum(precip_camp_attuale))
)

print(df_confronto, digits = 1)

cat("\n\nAnalisi completata. File generati:\n")
cat("- Climogramma_Salvatori.png\n")
cat("- Climogramma_Attuale.png\n")
