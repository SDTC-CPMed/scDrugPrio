#
#
# CORRELATIONS OF CELL TYPE HETEROGENEITY
#
# BY SAMUEL SCHAEFER
######################################################################

fp <- getwd()
setwd(paste(fp, "/..",sep=""))
dir.create(path = "Output/Comparison_heterogenity")

# Load RA
ra <- as.matrix(read.table(file = "Output/Cell_typing/OLD_RA_MOUSE_%_of_cells_of_individual_samples_per_cluster.txt",sep="\t",header = T))
rownames(ra) <- ra[,1]
ra <- ra[,-1]

# Load CD
cd <- as.matrix(read.table(file = "Input/CD GSE134809/Individual_patients/CD_%_of_cells_of_individual_samples_per_cluster.txt",sep="\t",header = T))
rownames(cd) <- cd[,1]
# alternative to above: cd <- read.table(file = "Individual_patients/CD_%_of_cells_of_individual_samples_per_cluster.txt", sep="\t", header = T)
colnames(cd)[1] <- "GEO_ID"
# rename values in column GEO_ID based on "inflammed" or "uninflammed" sample
temp <- cbind(c("GSM3972009","GSM3972010","GSM3972011","GSM3972012",
                "GSM3972013","GSM3972014","GSM3972016","GSM3972015",
                "GSM3972017","GSM3972018","GSM3972020","GSM3972019",
                "GSM3972022","GSM3972021","GSM3972024","GSM3972023",
                "GSM3972026","GSM3972025","GSM3972028","GSM3972027",
                "GSM3972030","GSM3972029"),c("Infl","Uninfl",
                                             "Infl","Uninfl","Infl","Uninfl","Infl","Uninfl",
                                             "Infl","Uninfl","Infl","Uninfl","Infl","Uninfl",
                                             "Infl","Uninfl","Infl","Uninfl","Infl","Uninfl",
                                             "Infl","Uninfl"))
cd[,1] <- temp[match(temp[,1], cd[,1]),2]
cd <- cbind(cd, NA)
cd[,ncol(cd)] <- paste("Patient", c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11),sep=" ")
colnames(cd)[ncol(cd)] <- "Patient"

# works not great as only smoothed lines are fitted
#install.packages("PerformanceAnalytics")
#library("PerformanceAnalytics")
#my_data <- mtcars[, c(1,3,4,5,6,7)]
#chart.Correlation(my_data, histogram=TRUE, pch=19)

# fits linear regression lines
chart.Correlation.linear <-   function (R, histogram = TRUE, method=c("pearson", "kendall", "spearman"), ...)
  { # @author R Development Core Team
    # @author modified by Peter Carl & Marek Lahoda
    # Visualization of a Correlation Matrix. On top the (absolute) value of the correlation plus the result 
    # of the cor.test as stars. On botttom, the bivariate scatterplots, with a linear regression fit. 
    # On diagonal, the histograms with probability, density and normal density (gaussian) distribution.
    
    x = checkData(R, method="matrix")
    
    if(missing(method)) method=method[1] #only use one
    cormeth <- method
    
    # Published at http://addictedtor.free.fr/graphiques/sources/source_137.R
    panel.cor <- function(x, y, digits=2, prefix="", use="pairwise.complete.obs", method=cormeth, cex.cor, ...)
    {
      usr <- par("usr"); on.exit(par(usr))
      par(usr = c(0, 1, 0, 1))
      r <- cor(x, y, use=use, method=method) # MG: remove abs here
      txt <- format(c(r, 0.123456789), digits=digits)[1]
      txt <- paste(prefix, txt, sep="")
      if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
      
      test <- cor.test(as.numeric(x),as.numeric(y), method=method)
      # borrowed from printCoefmat
      Signif <- symnum(test$p.value, corr = FALSE, na = FALSE,
                       cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                       symbols = c("***", "**", "*", ".", " "))
      # MG: add abs here and also include a 30% buffer for small numbers
      text(0.5, 0.5, txt, cex = cex * (abs(r) + .3) / 1.3)
      text(.8, .8, Signif, cex=cex, col=2)
    }
    
    #remove method from dotargs
    dotargs <- list(...)
    dotargs$method <- NULL
    rm(method)
    
    hist.panel = function (x, ...=NULL ) {
      par(new = TRUE)
      hist(x,
           col = "light gray",
           probability = TRUE,
           axes = FALSE,
           main = "",
           breaks = "FD")
      lines(density(x, na.rm=TRUE),
            col = "red",
            lwd = 1)
      # adding line representing density of normal distribution with parameters correponding to estimates of mean and standard deviation from the data 
      ax.x = seq(min(x), max(x), 0.1)                                                  # ax.x containts points corresponding to data range on x axis
      density.est = dnorm(ax.x, mean = mean(x), sd = sd(x))   # density corresponding to points stored in vector ax.x 
      lines(ax.x, density.est, col = "blue", lwd = 1, lty = 1)                                # adding line representing density into histogram
      rug(x)
    }
    
    # Linear regression line fit over points
    reg <- function(x, y, ...) {
      points(x,y, ...)
      abline(lm(y~x), col = "red") 
    }
    
    # Draw the chart
    if(histogram)
      pairs(x, gap=0, lower.panel=reg, upper.panel=panel.cor, diag.panel=hist.panel)
    else
      pairs(x, gap=0, lower.panel=reg, upper.panel=panel.cor) 
  }

# For RA sick
my_data <- t(ra[grepl(rownames(ra), pattern = "Sick"),])
mode(my_data) <- "numeric"

library("gplots")
balloonplot(x = as.table(as.matrix(t(my_data))), main = "Patients", xlab = "", ylab = "", label = F, show.margins = F)
#mosaicplot(as.table(as.matrix(my_data)), shade = T, las = 2, main = "Patients")

chisq.test(my_data) # P = 0.9919

# Correlation plot
pdf(file = "Output/Comparison_heterogenity/Sick_RA_mice_pearsson_correlation_cellular_pattern_performance_analytics.pdf",paper = "a4r")
chart.Correlation.linear(my_data, histogram=F, pch=19, method = "pearson")
dev.off()


# For CD sick
my_data <- cd[cd[,1]%in% "Infl",]
rownames(my_data) <- my_data[,ncol(my_data)]
my_data <- my_data[,-c(1,ncol(my_data))]
my_data <- t(my_data)
mode(my_data) <- "numeric"

#pdf(file = "Output/Comparison_heterogenity/Sick_CD_samples_cellular_pattern_baloon_plot.pdf", width = 8, height = 20)
#  balloonplot(x = as.table(as.matrix(t(my_data))), xlab = "", ylab = "", label = F, show.margins = F)
#dev.off()
#pdf(file = "Output/Comparison_heterogenity/Sick_CD_samples_cellular_pattern_mosaic_plots.pdf", width = 20, height = 10)
    #mosaicplot(as.table(as.matrix(my_data)), shade = T, las = 2, main = "Patients",  cex.axis = 1, off = 8, dir = c("h","v"))
#dev.off()

chisq <- chisq.test(my_data) # P = 2.2e-16
library(corrplot)
corrplot(t(chisq$residuals), is.cor = F, method = "color")


# Correlation plot
pdf(file = "Output/Comparison_heterogenity/Sick_samples_CD_patients_pearsson_correlation_cellular_pattern_performance_analytics.pdf",paper = "a4r")
chart.Correlation.linear(my_data, histogram=F, pch=19, method = "pearson")
dev.off()


