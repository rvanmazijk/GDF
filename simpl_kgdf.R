kGDF <- function(k = 1,
                 y,
                 fm,
                 data,
                 distr = "normal",
                 mod = "GLM",
                 noise = 0.25,
                 perf,
                 pass,
                 method = "jackknife") {

    # Sample size n
    n <- length(y)

    # Set noise
    if (distr == "normal") {
        noise <- noise * sd(y)
    }

    # dataframes for ye and yehat
    ye <- yehat <- data.frame(rep(NA, n))

    if (method == "jackknife") {
        # Jack-knife makes sure that each data point is perturbed per pass

        # loop to re-run model until each data point has been perturbed
        i <- 1  # subscription
        for (t in 1:pass) {

            # Check for NAs
            check.NA <- rep(NA, n)

            while (sum(is.na(check.NA)) > 0) {

                # Overwrite y2 each loop with original y
                y2 <- y

                # Check which data point still needs to be perturbed
                to.fill <- which(is.na(check.NA))

                # Randomise subset of size k to perturb during this loop
                if ((length(to.fill) == 1) & (k == 1)) {
                    # If there is only one NA to be disturbed,
                    # and sample size = k = 1,
                    # perturb last remaining point
                    y2[to.fill] <- y[to.fill] + rnorm(k, 0, noise)
                    check.NA[to.fill] <- 1  # mark perturbed subset
                } else{
                    if ((length(to.fill) >= k) & (length(to.fill) > 1)) {
                        # Randomly select k data points
                        to.use <- sample(
                            to.fill,
                            size = k,
                            replace = FALSE
                        )
                    }
                    if ((length(to.fill) < k)) {
                        # If there are fewer remaining NAs
                        # than data points to be disturbed,
                        # fill up sample of size k
                        fill.up <- sample(
                            1:n,
                            size = k - length(to.fill),
                            replace = FALSE
                        )
                        to.use <- c(to.fill, fill.up)
                    }
                    if (distr == "normal") {
                        y2[to.use] <- y2[to.use] + rnorm(k, 0, noise)
                    }
                    if (distr == "binary") {
                        y2[to.use] <- !y2[to.use]
                    }
                    check.NA[to.use] <- 1  # mark perturbed subset
                }

                # Provide "y2" in the global environment
                assign("y2", y2, envir = .GlobalEnv)

                # Accommodate model with perturbed data
                if (mod == "GLM" |
                    mod == "GAM" |
                    mod == "ANN") {
                    fm.updated <- update(fm, y2 ~ ., data = data)
                }
                if (mod == "rF") {
                    if (fm$type == "regression") {
                        fm.updated <- update(fm, y2 ~ ., data = data)
                    }
                    if (fm$type == "classification") {
                        fm.updated <- update(fm, as.factor(y2) ~ ., data = data)
                    }
                }
                if (mod == "BRT") {
                    fm.updated <- update(fm, as.formula(y2 ~ .), data = data)
                }

                # Store intermediates
                ye[, i] <- y2
                if (mod == "rF") {
                    if (fm$type == "regression") {
                        yehat[, i] <- predict(fm.updated, type = "response")
                    }
                    if (fm$type == "classification") {
                        yehat[, i] <- predict(fm.updated, type = "prob")[, 2]
                    }
                }
                if (mod == "ANN") {
                    yehat[, i] <- predict(fm.updated, type = "raw")
                }
                if (mod == "GLM" |
                    mod == "GAM") {
                    yehat[, i] <- predict(fm.updated, type = "response")
                }
                if (mod == "BRT") {
                    yehat[, i <- predict(
                        fm.updated,
                        type = "response",
                        n.trees = perf
                    )]
                }
                # Move one column further in the data frames
                i <- i + 1
                # Remove y2 from global environment
                rm(y2, envir = .GlobalEnv)
            }

        }

        cat("i =", i - 1)

    }

    if (method == "bootstrap") {
        # Meaning that the k perturbed data points are selected in every pass
        # from the complete dataset, possibly leaving out a few
        for (i in 1:pass) {

            # Perturbation
            y2 <- y  # Perturb copy of y
            to.perturb <- sample(n, size = k, replace = FALSE)
            if (distr == "normal") {
                y2[to.perturb] <- y2[to.perturb] + rnorm(k, 0, noise)
            }
            if (distr == "binary") {
                y2[to.perturb] <- !y2[to.perturb]
            }

            # For model update needed, apparently
            assign("y2", y2, envir = .GlobalEnv)

            # Accommodate model with perturbed data
            if (mod == "GLM" |
                mod == "GAM" |
                mod == "ANN") {
                fm.updated <- update(fm, y2 ~ ., data = data)
            }
            if (mod == "rF" &
                distr == "normal") {
                update(fm, y2 ~ ., data = data)
            }
            if (mod == "rF" &
                distr == "binary") {
                update(fm, as.factor(y2) ~ ., data = data)
            }
            if (mod == "BRT") {
                fm.updated <- update(fm, as.formula(y2 ~ .), data = data)
            }

            # Store intermediates
            ye[, i] <- y2
            if (mod == "rF") {
                if (fm$type == "regression") {
                    yehat[, i] <- predict(fm.updated, type = "response")
                }
                if (fm$type == "classification") {
                    yehat[, i] <- predict(fm.updated, type = "prob")[, 2]
                }
            }
            if (mod == "ANN") {
                yehat[, i] <- predict(fm.updated, type = "raw")
            }
            if (mod == "GLM" |
                mod == "GAM") {
                yehat[, i] <- predict(fm.updated, type = "response")
            }
            if (mod == "BRT") {
                yehat[, i <- predict(
                    fm.updated,
                    type = "response",
                    n.trees = perf
                )]
            }

            rm(y2, envir = .GlobalEnv)

        }
    }

    # LR: Horizontal Method:
    ye <- as.matrix(ye)
    yehat <- as.matrix(yehat)
    #Fit a LR to yehat[i,] vs ye[,i]
    fm.LR <-
        lapply(
            seq_along(ye[, 1]),
            function(x) lm(yehat[x, ] ~ ye[x, ])
        )
    #Calculate sum of slopes m[i] as GDF
    GDF.hor <-
        sum(
            sapply(
                fm.LR,
                function(i) coefficients(i)[2]
            ),
            na.rm = TRUE
        )

    return(GDF.hor)

}
