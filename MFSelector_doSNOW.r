mfselector <- function (data, nsc, stageord = F, stagename = F, type = 1, nline = T, dline = T, pdf = 1:100, cmp = 0, permut = 0, svdenoise = 0.03, svdetimes = 0, cores = detectCores())
{
    classes <- classes_all <- length(nsc)
    samples <- samples_all <- sum(nsc)
    nsc_all <- nsc

    if (stageord[1] == F) {
        stageord <- 1:samples
    }
    stageord_all <- stageord

    if (cmp != 0) {
        classes <- length(nsc) - cmp
        stageord <- stageord[1:sum(nsc[1:classes])]
        nsc <- nsc[1:classes]
        samples <- sum(nsc)
    }

    if (stagename[1] == F) {
        stagename <- paste("Class", seq(1, classes))
    }

    all_IDs <- as.character(data[, 1])
    data_0 <- data[, -1]
    mydata <- data_0[, stageord_all]

    nor <- function(x) {
        min <- min(x)
        max <- max(x)
        fun <- (x - min)/(max - min)
        return(fun)
    }

    nor_train <- t(apply(mydata, 1, nor))
    Class_Tag <- vector(mode = "list", length = classes)
    Class_Tag[[1]] <- c(1:nsc[1])
    for (j in 2:classes) {
        Class_Tag[[j]] <- c(c(Class_Tag[[j - 1]][nsc[j - 1]] +
            1):c(Class_Tag[[j - 1]][nsc[j - 1]] + nsc[j]))
    }

    t_data <- data.frame(t(mydata[, 1:samples]))

    mainfun <- function(dataX, samples, classes, Class_Tag, type) {
        vec_min_bct <- numeric()
        point <- numeric()
        flags <- rep(0, samples)
        for (vs in 1:(classes - 1)) {
            one <- 1:sum(nsc[1:vs])
            all <- (sum(nsc[1:vs]) + 1):samples
            vec_bct <- numeric()
            for (s in 1:length(one)) {
                if (type == 1) {
                  temp_bct_one <- sum(flags[one][which(dataX[s] > dataX[one])] == 0)
                  temp_bct_all <- sum(flags[all][which(dataX[s] < dataX[all])] == 0)
                }
                else if (type == 2) {
                  temp_bct_one <- sum(flags[one][which(dataX[s] < dataX[one])] == 0)
                  temp_bct_all <- sum(flags[all][which(dataX[s] > dataX[all])] == 0)
                }
                temp_bct <- temp_bct_one + temp_bct_all
                vec_bct <- c(vec_bct, temp_bct)
            }
            min_idx <- org_min_idx <- rep(NA, vs)
            org_min_idx[1] <- which.min(vec_bct[Class_Tag[[1]]])
            for (i in 1:vs) {
                min_idx[i] <- which.min(vec_bct[Class_Tag[[i]]])
                if (i != 1) {
                  org_min_idx[i] <- sum(nsc[1:(i - 1)]) + min_idx[i]
                }
            }
            min_in_which_group <- order(vec_bct[org_min_idx], decreasing = T)[vs]
            if (min_in_which_group == 1) {
                min_idx_in_higher_class <- which.min(vec_bct)
            }
            else {
                min_idx_in_higher_class <- sum(nsc[1:(min_in_which_group - 1)]) + which.min(vec_bct[Class_Tag[[min_in_which_group]]])
            }
            vec_min_bct <- c(vec_min_bct, vec_bct[min_idx_in_higher_class])
            point <- c(point, min_idx_in_higher_class)
            if (type == 1) {
                flags[one] <- flags[one] + (dataX[min_idx_in_higher_class] > dataX[one])
                flags[all] <- flags[all] + (dataX[min_idx_in_higher_class] < dataX[all])
            }
            else if (type == 2) {
                flags[one] <- flags[one] + (dataX[min_idx_in_higher_class] < dataX[one])
                flags[all] <- flags[all] + (dataX[min_idx_in_higher_class] > dataX[all])
            }
        }
        vec_BCT <- sum(vec_min_bct)
        return(c(vec_BCT, point))
    }

    all_info <- mclapply(t_data, mainfun, samples, classes, Class_Tag, type, mc.cores = cores)

    nlines <- numeric(length = length(all_info))
    all_DE <- vector(length = length(all_info))
    for (i in 1:length(all_info)) {
        nlines[i] <- (classes - 1) - sum(duplicated(as.numeric(mydata[i, all_info[[i]][-1]])))
        all_DE[i] <- all_info[[i]][1]
    }

    if (nline == T) {
        check_idx <- which(nlines == (classes - 1))
        withline <- rep("TRUE", length(check_idx))
    }
    else if (nline == F) {
        check_idx <- which(nlines != (classes - 1))
        withline <- rep("FALSE", length(check_idx))
    }
    else if (nline == "both") {
        check_idx <- 1:length(all_info)
        withline <- character(length = length(all_info))
        check_idx1 <- which(nlines == (classes - 1))
        check_idx2 <- which(nlines != (classes - 1))
        withline[check_idx1] <- "TRUE"
        withline[check_idx2] <- "FALSE"
    }

    if (permut != F) {
        arranged_data <- mydata[, 1:samples]
        for (times in 1:permut) {
            class_order_permu <- sample(sum(nsc))
            arranged_data_permu <- arranged_data[, class_order_permu]
            permu_data <- data.frame(t(arranged_data_permu))
            permu_all_info <- t(as.data.frame(mclapply(permu_data, mainfun, samples, classes, Class_Tag, type, mc.cores = cores)))
            all_DE <- cbind(all_DE, permu_all_info[, 1])
        }
        DE_pvalue <- rep(NA, nrow(permu_all_info))
        DE_qvalue <- rep(NA, nrow(permu_all_info))
        for (i in 1:nrow(permu_all_info)) {
            DE_pvalue[i] <- length(which(all_DE[, -1] <= all_DE[i, 1]))/(nrow(permu_all_info) * permut)
            DE_qvalue[i] <- min(length(which(all_DE[, -1] <= all_DE[i, 1]))/(length(which(all_DE[, 1] <= all_DE[i, 1])) * permut), 1)
        }
        DE <- all_DE[, 1]
    }
    else {
        DE <- all_DE
    }

    if (svdetimes != F) {
        diff_sqr <- numeric()
        for (j in 1:svdetimes) {
            data_noise <- nor_train[, 1:samples] + runif(dim(nor_train)[1] * samples, min = 0, max = 1) * svdenoise * sample(c(1, -1), dim(nor_train)[1] * samples, replace = T)
            t_nor_train <- as.data.frame(apply(data_noise, 1, nor))
            all_info_noise <- t(as.data.frame(mclapply(t_nor_train, mainfun, samples, classes, Class_Tag, type, mc.cores = cores)))
            diff_sqr <- cbind(diff_sqr, (DE - all_info_noise[, 1])^2)
        }
        conf <- function(x) {
            sum(x)/svdetimes
        }
        SVDE <- apply(diff_sqr, 1, conf)
    }

    DE_max <- max(DE)
    DE_min <- min(DE)
    if (permut == F) {
        new_raw <- NULL
        if (svdetimes == F) {
            new_raw_temp <- cbind(all_IDs[check_idx], DE[check_idx], withline)
            colnames(new_raw_temp) <- c("ID", "DE", "with N-1 distinct lines")
            sort_idx <- order(as.numeric(new_raw_temp[, 2]), decreasing = F)
            new_raw <- new_raw_temp[sort_idx, ]
        }
        else {
            new_raw_temp <- cbind(all_IDs[check_idx], DE[check_idx], SVDE[check_idx], withline)
            colnames(new_raw_temp) <- c("ID", "DE", "SVDE", "with N-1 distinct lines")
            sort_idx <- order(as.numeric(new_raw_temp[, 2]), decreasing = F)
            new_raw_temp <- new_raw_temp[sort_idx, ]
            for (i in DE_min:DE_max) {
                sort_idx <- order(as.numeric(new_raw_temp[which(new_raw_temp[, 2] == i), 3]), decreasing = F)
                new_raw <- rbind(new_raw, new_raw_temp[which(new_raw_temp[, 2] == i)[sort_idx], ])
            }
        }
    }
    else {
        new_raw <- NULL
        if (svdetimes == F) {
            new_raw_temp <- cbind(all_IDs[check_idx], DE[check_idx], DE_pvalue[check_idx], DE_qvalue[check_idx], withline)
            colnames(new_raw_temp) <- c("ID", "DE", "p-value", "q-value", "with N-1 distinct lines")
            sort_idx <- order(as.numeric(new_raw_temp[, 2]), decreasing = F)
            new_raw <- new_raw_temp[sort_idx, ]
        }
        else {
            new_raw_temp <- cbind(all_IDs[check_idx], DE[check_idx], DE_pvalue[check_idx], DE_qvalue[check_idx], SVDE[check_idx], withline)
            colnames(new_raw_temp) <- c("ID", "DE", "p-value", "q-value", "SVDE", "with N-1 distinct lines")
            sort_idx <- order(as.numeric(new_raw_temp[, 2]), decreasing = F)
            new_raw_temp <- new_raw_temp[sort_idx, ]
            for (i in DE_min:DE_max) {
                sort_idx <- order(as.numeric(new_raw_temp[which(new_raw_temp[, 2] == i), 5]), decreasing = F)
                new_raw <- rbind(new_raw, new_raw_temp[which(new_raw_temp[, 2] == i)[sort_idx], ])
            }
        }
    }
    sort_idx <- match(as.character(new_raw[, 1]), all_IDs)
    if (exists("txtname") == F) {
        clock_txt <- paste(substring(Sys.time(), 12, 13), substring(Sys.time(), 15, 16), substring(Sys.time(), 18, 19), "txt", sep = ".")
        txtname <- paste("mfselector", Sys.Date(), clock_txt, sep = "_")
    }
    write.table(new_raw, txtname, sep = "\t", row.names = F, quote = F)
    color_vector <- rep(1, nsc_all[1])
    for (i in 2:classes_all) {
        color_vector <- c(color_vector, rep(i, nsc_all[i]))
    }
    ylab_name <- "Relative Gene Expression Values"
    xlab_name <- "Sample Index"
    if (samples_all <= 12) {
        x_axis <- 1.5
    }
    else {
        x_axis <- 0.5
    }
    if (type == 1) {
        y_axis <- 0.3
    }
    else if (type == 2) {
        y_axis <- 0.97
    }
    if (exists("pdfname") == F) {
        clock_pdf <- paste(substring(Sys.time(), 12, 13), substring(Sys.time(), 15, 16), substring(Sys.time(), 18, 19), "pdf", sep = ".")
        pdfname <- paste("mfselector", Sys.Date(), clock_pdf, sep = "_")
    }
    pdf(pdfname)
    if (pdf[1] == "all") {
        for (i in 1:dim(new_raw)[1]) {
            if ((svdetimes != F) && (permut != F)) {
                fname <- paste("Top ", i, ": ", new_raw[i, 1], ", DE: ", new_raw[i, 2], ", p-value: ", round(new_raw[i, 3], 5), ", SVDE: ", new_raw[i, 5], sep = "")
            }
            else if ((svdetimes != F) && (permut == F)) {
                fname <- paste("Top ", i, ": ", new_raw[i, 1], ", DE: ", new_raw[i, 2], ", SVDE: ", new_raw[i, 5], sep = "")
            }
            else if ((svdetimes == F) && (permut != F)) {
                fname <- paste("Top ", i, ": ", new_raw[i, 1], ", DE: ", new_raw[i, 2], ", p-value: ", round(new_raw[i, 3], 5), sep = "")
            }
            else {
                fname <- paste("Top ", i, ": ", new_raw[i, 1], ", DE: ", new_raw[i, 2], sep = "")
            }
            plot(nor_train[sort_idx[i], ], col = color_vector, pch = color_vector, xlab = xlab_name, ylab = ylab_name, cex = 1, lwd = 1)
            legend(x_axis, y_axis, stagename, pch = 1:classes_all, col = 1:classes_all, cex = 1, pt.lwd = 1.5)
            title(main = list(fname, cex = 1))
            if (dline == T) {
                for (j in 2:classes_all) {
                  abline(h = nor_train[sort_idx[i], all_info[sort_idx][[i]][j]], col = j - 1, lty = 2)
                }
            }
        }
    }
    else {
        if (length(pdf) == 1) {
            if ((svdetimes != F) && (permut != F)) {
                fname <- paste("Top ", pdf, ": ", new_raw[pdf, 1], ", DE: ", new_raw[pdf, 2], ", p-value: ", round(as.numeric(new_raw[pdf, 3]), 6), ", SVDE: ", new_raw[pdf, 5], sep = "")
            }
            else if ((svdetimes != F) && (permut == F)) {
                fname <- paste("Top ", pdf, ": ", new_raw[pdf, 1], ", DE: ", new_raw[pdf, 2], ", SVDE: ", new_raw[pdf, 3], sep = "")
            }
            else if ((svdetimes == F) && (permut != F)) {
                fname <- paste("Top ", pdf, ": ", new_raw[pdf, 1], ", DE: ", new_raw[pdf, 2], ", p-value: ", round(as.numeric(new_raw[pdf, 3]), 6), sep = "")
            }
            else {
                fname <- paste("Top ", pdf, ": ", new_raw[pdf, 1], ", DE: ", new_raw[pdf, 2], sep = "")
            }
            plot(nor_train[sort_idx[pdf], ], col = color_vector, pch = color_vector, xlab = xlab_name, ylab = ylab_name, cex = 1, lwd = 1, cex.lab = 1.2, cex.axis = 1.2, font.lab = 1.2, font.axis = 1.2)
            legend(x_axis, y_axis, stagename, pch = 1:classes_all, col = 1:classes_all, cex = 1, pt.lwd = 1.5)
            title(main = list(fname, cex = 1))
            if (dline == T) {
                for (j in 2:classes_all) {
                  abline(h = nor_train[sort_idx[pdf], all_info[sort_idx][[pdf]][j]], col = j - 1, lty = 2)
                }
            }
        }
        else {
            if (length(pdf) <= dim(new_raw)[1]) {
                for (i in pdf[1]:length(pdf)) {
                  if ((svdetimes != F) && (permut != F)) {
                    fname <- paste("Top ", i, ": ", new_raw[i, 1], ", DE: ", new_raw[i, 2], ", p-value: ", round(as.numeric(new_raw[i, 3]), 6), ", SVDE: ", new_raw[i, 5], sep = "")
                  }
                  else if ((svdetimes != F) && (permut == F)) {
                    fname <- paste("Top ", i, ": ", new_raw[i, 1], ", DE: ", new_raw[i, 2], ", SVDE: ", new_raw[i, 3], sep = "")
                  }
                  else if ((svdetimes == F) && (permut != F)) {
                    fname <- paste("Top ", i, ": ", new_raw[i, 1], ", DE: ", new_raw[i, 2], ", p-value: ", round(as.numeric(new_raw[i, 3]), 6), sep = "")
                  }
                  else {
                    fname <- paste("Top ", i, ": ", new_raw[i, 1], ", DE: ", new_raw[i, 2], sep = "")
                  }
                  plot(nor_train[sort_idx[i], ], col = color_vector, pch = color_vector, xlab = xlab_name, ylab = ylab_name, cex = 1, lwd = 1, cex.lab = 1.2, cex.axis = 1.2, font.lab = 1.2, font.axis = 1.2)
                  legend(x_axis, y_axis, stagename, pch = 1:classes_all, col = 1:classes_all, cex = 1, pt.lwd = 1.5)
                  title(main = list(fname, cex = 1))
                  if (dline == T) {
                    for (j in 2:classes_all) {
                      abline(h = nor_train[sort_idx[i], all_info[sort_idx][[i]][j]], col = j - 1, lty = 2)
                    }
                  }
                }
            }
            else {
                for (i in 1:dim(new_raw)[1]) {
                  if ((svdetimes != F) && (permut != F)) {
                    fname <- paste("Top ", i, ": ", new_raw[i, 1], ", DE: ", new_raw[i, 2], ", p-value: ", round(new_raw[i, 3], 5), ", SVDE: ", new_raw[i, 5], sep = "")
                  }
                  else if ((svdetimes != F) && (permut == F)) {
                    fname <- paste("Top ", i, ": ", new_raw[i, 1], ", DE: ", new_raw[i, 2], ", SVDE: ", new_raw[i, 5], sep = "")
                  }
                  else if ((svdetimes == F) && (permut != F)) {
                    fname <- paste("Top ", i, ": ", new_raw[i, 1], ", DE: ", new_raw[i, 2], ", p-value: ", round(new_raw[i, 3], 5), sep = "")
                  }
                  else {
                    fname <- paste("Top ", i, ": ", new_raw[i, 1], ", DE: ", new_raw[i, 2], sep = "")
                  }
                  plot(nor_train[sort_idx[i], ], col = color_vector, pch = color_vector, xlab = xlab_name, ylab = ylab_name, cex = 1, lwd = 1)
                  legend(x_axis, y_axis, stagename, pch = 1:classes_all, col = 1:classes_all, cex = 1, pt.lwd = 1.5)
                  title(main = list(fname, cex = 1))
                  if (dline == T) {
                    for (j in 2:classes_all) {
                      abline(h = nor_train[sort_idx[i], all_info[sort_idx][[i]][j]], col = j - 1, lty = 2)
                    }
                  }
                }
            }
        }
    }
    dev.off()
}
