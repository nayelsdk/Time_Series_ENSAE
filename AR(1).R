# Dictionnaire pour stocker les résultats du modèle AR(1)
results_ar1 <- list()

# Dataframe pour stocker les erreurs du modèle AR(1)
error_metrics_ar1 <- data.frame(Variable = character(), MSE = numeric(), RMSE = numeric(), MAE = numeric(), MAPE = numeric(), stringsAsFactors = FALSE)

# ---------------------- Boucle pour AR(1) ----------------------
for (var in variables) {
  # Création de la série temporelle
  ts_train <- to_ts(train_data[[var]])

  # ---------------------- AR(1) ----------------------
  fit_ar1 <- Arima(ts_train, order = c(1, 0, 0)) # Modèle AR(1)
  forecast_ar1 <- forecast(fit_ar1, h = 12)
  predictions_ar1 <- forecast_ar1$mean

  # ---------------------- Comparaison ----------------------
  actual_values <- test_data[[var]]

  # Table de comparaison AR(1)
  comparison_ar1 <- data.frame(
    Date          = test_data$Date,
    Observed      = actual_values,
    Predicted_AR1 = predictions_ar1
  )
  comparison_ar1$Error_AR1 <- comparison_ar1$Observed - comparison_ar1$Predicted_AR1

  # ---------------------- Calcul des erreurs ----------------------
  errors_ar1 <- comparison_ar1 %>%
    filter(!is.na(Error_AR1)) %>%
    pull(Error_AR1)
  mse_ar1 <- mean(errors_ar1^2, na.rm = TRUE)
  rmse_ar1 <- sqrt(mse_ar1)
  mae_ar1 <- mean(abs(errors_ar1), na.rm = TRUE)
  obs_ar1 <- comparison_ar1 %>%
    filter(!is.na(Observed)) %>%
    pull(Observed)
  mape_ar1 <- mean(abs(errors_ar1 / obs_ar1) * 100, na.rm = TRUE)

  # ---------------------- Stocker les résultats ----------------------
  results_ar1[[var]] <- list(
    "Model" = fit_ar1,
    "Forecast" = forecast_ar1,
    "Comparison" = comparison_ar1,
    "MSE" = mse_ar1,
    "RMSE" = rmse_ar1,
    "MAE" = mae_ar1,
    "MAPE" = mape_ar1
  )

  # Ajouter les résultats au dataframe des erreurs
  error_metrics_ar1 <- rbind(error_metrics_ar1, data.frame(Variable = var, MSE = mse_ar1, RMSE = rmse_ar1, MAE = mae_ar1, MAPE = mape_ar1))

  # ---------------------- Affichage des résultats ----------------------
  print(paste("Résultats AR(1) pour :", var))
  print(comparison_ar1)
  print(paste("MSE:", mse_ar1, "RMSE:", rmse_ar1, "MAE:", mae_ar1, "MAPE:", mape_ar1, "%"))
}

# ---------------------- Comparaison globale ----------------------
print("Comparaison des performances entre SARIMA (existant) et AR(1) :")
print("AR(1):")
print(error_metrics_ar1)
