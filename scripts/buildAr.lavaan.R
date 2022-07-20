library(lavaan)

buildAr <- function(waves, varName) {
  ## Create AR Variables
  arVars <- NULL
  for (i in waves) {
    arVars <- paste0(
      arVars,
      paste0(
        "wx",
        i,
        " =~ 1*",
        varName,
        i,
        "\n"
      )
    )
  }
  stabilities <- NULL
  for (i in waves[-1]) {
    stabilities <- paste0(
      stabilities,
      paste0(
        "wx",
        i,
        " ~ a*wx",
        (i - 1),
        "\n"
      )
    )
  }
  variances <- paste0(
    "wx",
    waves[1],
    " ~~ arx * wx",
    waves[1],
    "\n"
  )
  for (i in waves[-1]) {
    variances <- paste0(
      variances,
      paste0(
        "wx",
        i,
        " ~~ arx2 * wx",
        i,
        "\n"
      )
    )
  }
  finalModel <- paste0(
    arVars, stabilities, variances,
    paste0(
      "## Constraints\n",
      "arx2 == 1 - (arx * a^2)"
    )
  )
  return(finalModel)
}

buildAr2 <- function(waves, varName) {
  stabilities <- NULL
  for (i in 2:length(waves)) {
    stabilities <- paste0(
      stabilities,
      paste0(
        varName,
        waves[i],
        " ~ ",
        varName,
        waves[i - 1],
        "\n"
      )
    )
  }
  variances <- NULL
  for (i in 1:length(waves)) {
    variances <- paste0(
      variances,
      paste0(
        varName,
        waves[i],
        " ~~ ",
        varName,
        waves[i],
        "\n"
      )
    )
  }
  finalModel <- paste0(
    stabilities, variances
  )
  return(finalModel)
}
