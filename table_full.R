library(flextable)
library(officer)

# --- Deine Dataframes ---
frequency <- c("alpha","beta low","beta high","delta","gamma low","gamma mid", 
               "gamma high", "theta")

# Stability
p_value1 <- c(1.4101e-4, 0.0016, 0.2085, 0.0022, 0.3903, 0.0012, 0.1215, 9.0308e-5)
eta_square1 <- c(0.1705, 0.1207, 0.0202, 0.1142, 0.0095, 0.1260, 0.0305, 0.1794)
t_statistics1 <- c(-4.7399, -5.2914, -3.259, -5.2543, -2.2601, 5.961, 2.1404, -5.0532)
mean_s_1 <- c(0.041869995768794, 1.523100244513761, 0.205113060430492, 
              0.00628549207572607, 5.94804520338910, 10.0788106233599, 
              7.05073894963155, 0.0298218846479384)
mean_s_2 <- c(0.0493780940816792, 1.58276034639412, 0.227940595332729, 
              0.00671766997649478, 6.016272786091022, 9.644579219495318,
              6.887590727956090, 0.032061436388776)
sd_s_1 <- c(0.067194390437282, 0.208373771439243, 1.426219231752598,
            0.006432279920756, 4.939460290674851, 7.670914178552150,
            4.924776936808299, 0.035955313191309)
sd_s_2 <- c(0.068388978118039, 0.218808581653669, 1.361984108948688,
            0.006723483459609, 4.891301238248282, 7.506546685826247,
            4.936042726366757, 0.037577613890711)
d1 <- data.frame(frequency, p = p_value1, eta2 = eta_square1, t = t_statistics1)
dfmean_1 <- data.frame(frequency, M1=mean_s_1, M2=mean_s_2, sd1=sd_s_1, sd2=sd_s_2)


# Transition Energy
p_value2 <- c(2.2377e-6, 6.7584e-6, 2.8921e-4, 0.0092, 0.4405, 1.9454e-6, 1.6250e-4, 0.0146)
eta_square2 <- c(0.2508, 0.2300, 0.1559, 0.0839, 0.0076, 0.2534, 0.1676, 0.0740)
t_statistics2 <- c(-3.9776, -4.8086, -3.5841, -5.1109, -2.1418, 0.4682, -1.1134, -3.9412)
mean_t_1 <- c(89.827978880622070, 11.989352646175742, 1.242609995093793,
              2.354968985514337e+02, 0.263902151050484, 0.179916509067566,
              0.293768968286261, 83.247444783200650)
mean_t_2 <- c(1.704883853391501e+02, 17.387306800605360, 1.552139133347837,
              3.053955695993724e+02, 0.291392674743452, 0.172369182053099,
              0.353060544757536, 1.320825528576605e+02)
sd_t_1 <- c(1.327307335647365e+02, 16.089314548945808, 1.456430612981913,
            1.487145408376230e+02, 0.209838562686656, 0.148843973382010,
            0.376087616764598, 1.141257512714554e+02)
sd_t_2 <- c(2.720351649341163e+02, 21.580211991862626, 1.873211589500833, 
            2.336893622798761e+02, 0.264677430248447, 0.230771246907799,
            0.765464679076269, 2.063026440064077e+02)
dfmean_2 <- data.frame(frequency, M1=mean_t_1, M2=mean_t_2, sd1=sd_t_1, sd2=sd_t_2)
d2 <- data.frame(frequency, p = p_value2, eta2 = eta_square2, t = t_statistics2)

# --- Funktion zum Formatieren der p-Werte ---
format_pvalue <- function(x) {
  stars <- ifelse(is.na(x), "",
                  ifelse(x < 0.001, "***",
                         ifelse(x < 0.01,  "**",
                                ifelse(x < 0.05,  "*", ""))))
  
  if (is.na(x)) return("")
  
  if (x < 0.001) {
    paste0("<.001", stars)
  } else {
    paste0(formatC(x, format = "f", digits = 3), stars)
  }
}

# --- Kombiniertes Dataframe mit leeren Spalten ---
df_combined <- data.frame(
  Frequency = d1$frequency,
  p_Stability = sapply(d1$p, format_pvalue),
  p_Transition = sapply(d2$p, format_pvalue),
  empty1 = "",                       # leere Spalte zwischen p und eta²
  eta2_Stability = round(d1$eta2, 2),
  eta2_Transition = round(d2$eta2, 2),
  empty2 = "",                       # leere Spalte zwischen eta² und t
  t_Stability = round(d1$t, 2),
  t_Transition = round(d2$t, 2),
  check.names = FALSE
)

# --- Flextable erstellen ---
ft <- flextable(df_combined)
ft <- theme_vanilla(ft)
ft <- set_caption(ft, "Table 1\nStatistical Results for Stability and Transition Energy (APA Style)")

# Gruppen-Header
ft <- add_header_row(
  ft,
  values = c("", "p", "", "η²", "", "t"),
  colwidths = c(1, 2, 1, 2, 1, 2)   # Summe = 9 Spalten
)

# Header-Labels, leere Spalten mit "" überschreiben
ft <- set_header_labels(ft,
                        p_Stability = "Stability",
                        p_Transition = "Transition Energy",
                        empty1 = "",                # leere Spalte
                        eta2_Stability = "Stability",
                        eta2_Transition = "Transition Energy",
                        empty2 = "",                # leere Spalte
                        t_Stability = "Stability",
                        t_Transition = "Transition Energy")


# Zentrieren und Spaltenbreite
ft <- align(ft, align = "center", part = "all")
ft <- width(ft, j = c("Frequency"), width = 1.5)           # erste Spalte
ft <- width(ft, j = c("p_Stability", "p_Transition"), width = 1.2)  # p-Werte
ft <- width(ft, j = c("empty1"), width = 0.3)             # leere Spalte
ft <- width(ft, j = c("eta2_Stability", "eta2_Transition"), width = 1.2)  # eta²
ft <- width(ft, j = c("empty2"), width = 0.3)             # leere Spalte
ft <- width(ft, j = c("t_Stability", "t_Transition"), width = 1.2)        # t-Werte

# Dünne Linie definieren
thin_border <- fp_border(color = "black", width = 1)  # width in points

# Linie unter Header
ft <- hline(ft, i = 1, part = "header", border = thin_border)

# APA-Footer
ft <- add_footer_lines(ft,
                       "Note. p = p-value; *, **, *** indicate p < .05, .01, .001; η² = eta squared.")

# --- Export nach Word ---
doc <- read_docx()
doc <- body_add_flextable(doc, value = ft)
print(doc, target = "APA_Table.docx")

## Tabelle für die Mittelwerte

# --- Neues Dataframe mit Paarung M und SD ---
df_means_grouped <- data.frame(
  Frequency = dfmean_1$frequency,
  M1_s = round(dfmean_1$M1, 2),
  SD1_s = round(dfmean_1$sd1, 2),
  M2_s = round(dfmean_1$M2, 2),
  SD2_s = round(dfmean_1$sd2, 2),
  M1_t = round(dfmean_2$M1, 2),
  SD1_t = round(dfmean_2$sd1, 2),
  M2_t = round(dfmean_2$M2, 2),
  SD2_t = round(dfmean_2$sd2, 2),
  check.names = FALSE
)

# --- Flextable erstellen ---
ft_means2 <- flextable(df_means_grouped)
ft_means2 <- theme_vanilla(ft_means2)
ft_means2 <- set_caption(ft_means2,
                         "Table 2\nMeans (M) and Standard Deviations (SD) for Stability and Transition Energy"
)

# Header: Zwei Zeilen -> obere Zeile für Gruppierung, zweite Zeile für M/SD
ft_means2 <- add_header_row(
  ft_means2,
  values = c("", "Stability M1", "Stability M2", "Transition M1", "Transition M2"),
  colwidths = c(1, 2, 2, 2, 2)   # Summe = 9 Spalten
)

# Untere Headerlabels (M/SD)
ft_means2 <- set_header_labels(ft_means2,
                               M1_s = "M", SD1_s = "SD",
                               M2_s = "M", SD2_s = "SD",
                               M1_t = "M", SD1_t = "SD",
                               M2_t = "M", SD2_t = "SD"
)

# Zentrieren und Breite anpassen
ft_means2 <- align(ft_means2, align = "center", part = "all")
ft_means2 <- width(ft_means2, j = "Frequency", width = 1.5)
ft_means2 <- width(ft_means2, j = 2:9, width = 1.0)

# Linie unter oberem Header
thin_border <- fp_border(color = "black", width = 1)
ft_means2 <- hline(ft_means2, i = 1, part = "header", border = thin_border)

# --- Export nach Word ---
doc <- read_docx()
doc <- body_add_flextable(doc, value = ft_means2)
print(doc, target = "APA_Table_Means_Grouped.docx")
