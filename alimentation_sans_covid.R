# -------------------------------
# 1) PARTIE ENTRAINEMENT : jusqu'à fin 2017
train_data <- data %>%
  filter(Date <= as.Date("2017-12-31"))

# 2) PARTIE TEST : uniquement l'année 2018
test_data <- data %>%
  filter(Date >= as.Date("2018-01-01") & Date <= as.Date("2018-12-31"))

# ----------------------------------------------------
# Création des séries temporelles au format ts
# ----------------------------------------------------
# Extraire l'année et le mois de la toute première observation d'entraînement
start_year  <- as.numeric(format(min(train_data$Date), "%Y"))
start_month <- as.numeric(format(min(train_data$Date), "%m"))

to_ts <- function(vec){
  ts(vec, start = c(start_year, start_month), frequency = 12)
}

# Exemple : on cible la colonne "Alimentation"
ts_Alimentation_train <- to_ts(train_data$Alimentation)

# ----------------------------------------------------
# Ajuster un modèle SARIMA (via auto.arima)
# ----------------------------------------------------
fit_alim <- auto.arima(ts_Alimentation_train)

# h = 12 pour prédire 12 mois (janv 2018 à déc 2018)
fc_alim <- forecast(fit_alim, h = 12)

# On récupère le vecteur des prévisions "ponctuelles"
pred_alim <- fc_alim$mean

# ----------------------------------------------------
# Comparer les prédictions à la réalité sur 2018
# ----------------------------------------------------
# test_data contient (normalement) 12 lignes (Janv 2018 à Déc 2018)
# On récupère la série observée (Alimentation) pour 2018
actual_alim_2018 <- test_data$Alimentation

# (Optionnel) Si votre forecast a exactement 12 points, 
# et que test_data a 12 lignes, pas besoin d'ajouter NA
# actual_alim_2018 <- c(actual_alim_2018, NA) 
# => seulement si vous aviez un décalage. 

# Idem pour les Dates
dates_2018 <- test_data$Date

# Construire la table de comparaison
comparison <- data.frame(
  Date      = dates_2018,
  Observed  = actual_alim_2018,
  Predicted = pred_alim
)

# Calcul des erreurs
comparison$Error <- comparison$Observed - comparison$Predicted

print(comparison)

# ----------------------------------------------------
# Indicateurs d'évaluation : MSE, RMSE, MAE, MAPE
# ----------------------------------------------------
# Extraire les erreurs
errors <- comparison %>%
  filter(!is.na(Error)) %>%
  pull(Error)

# MSE
mse <- mean(errors^2)
# RMSE
rmse <- sqrt(mse)
# MAE
mae <- mean(abs(errors))

# MAPE (en pourcent)
obs <- comparison %>%
  filter(!is.na(Observed)) %>%
  pull(Observed)

mape <- mean(abs(errors / obs) * 100)

# Affichage
mse
rmse
mae
mape

# ----------------------------------------------------
# (Optionnel) Visualiser
# ----------------------------------------------------
library(ggplot2)

ggplot(comparison, aes(x = Date)) +
  geom_line(aes(y = Observed,  color = "Observé"), size = 1) +
  geom_line(aes(y = Predicted, color = "Prédit"),  size = 1) +
  scale_color_manual(
    name = "Légende",
    values = c("Observé" = "blue", "Prédit" = "red")
  ) +
  labs(
    title = "Comparaison prévision SARIMA vs Observé (Alimentation, 2018)",
    x = "Date",
    y = "Valeur"
  ) +
  theme_minimal()



