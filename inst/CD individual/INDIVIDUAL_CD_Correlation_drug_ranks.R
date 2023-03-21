
#######################################################################################
# Final drug ranking evaluation INDIVIDUAL CD rank vs rank plots patients
#
# BY SAMUEL SCHAEFER
#######################################################################################

fp <- getwd()
setwd(paste(fp, "/..", sep=""))
library(doParallel)
library(readxl)
#library(data.table)
library(ggplot2)

# INPUT DATA 
#######################################################################################

pat <- paste("Patient_",1:11,sep="")
drugs <- read.table(file = "Input/Drugs/drugs_with_identical_drug_target_combinations_literature_PPI.txt", sep="\t", header = T, stringsAsFactors = F)
drugs <- unique(c(drugs[,1], drugs[,2]))
in_path <- "Output/CD/Individual_patients/DCA_MAST_DEGs_predictions/"
cd_drugs <- unique(read.delim(file = "Input/Drugs/CD_drugs_from_DrugBank_201015.txt", sep="\t", header = T)[,1])

inter_cents <- foreach(i = c(1:length(pat)), .combine = "cbind") %do% {
  temp <- as.matrix(read_xlsx(path = paste(in_path, pat[i], "/Drug_ranking/FINAL_drug_ranking_evaluated.xlsx",sep=""),sheet = 1, col_names = T))
  temp <- temp[,c("DrugBank", "combined_centrality_score")]
  temp <- temp[match(drugs, temp[,1]),2]
  return(temp)
}
rownames(inter_cents) <- drugs
colnames(inter_cents) <- pat
inter_cents <- inter_cents[rowSums(!is.na(inter_cents))>0,]

drug_ranks <- foreach(i =c(1:ncol(inter_cents)), .combine = "cbind") %do% {
  return(rank( - as.numeric(inter_cents[,i]), na.last = "keep"))
}
dimnames(drug_ranks) <- dimnames(inter_cents)
drug_ranks <- as.data.frame(drug_ranks)

drug_ranks <- cbind(drug_ranks, CD_drug = factor(x = rownames(drug_ranks)%in% cd_drugs, levels = c(F, T)))

my_cols <- c("grey50", "red")
my_shapes <- c(19,17)
my_alpha <- c(0.5,1)
my_cex <- c(0.5,1.5)

# Panel background color
library(RColorBrewer)
cols = brewer.pal(11, "RdBu")   # goes from red to white to blue
pal = colorRampPalette(cols)
cor_colors = data.frame(correlation = seq(-1,1,0.01), 
                        correlation_color = pal(201)[1:201])  # assigns a color for each r correlation value
cor_colors$correlation_color = as.character(cor_colors$correlation_color)



# Plot data
#######################################################################################
pairs(drug_ranks[,1:11], 
      pch = my_shapes[drug_ranks$CD_drug],
      cex = my_cex[drug_ranks$CD_drug],
      lower.panel = NULL,
      col = my_cols[drug_ranks$CD_drug])





# Correlation panel
panel.cor <- function(x, y){
  
  #usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  if(length(intersect(which(!is.na(x)), which(!is.na(y))))>5){
    r <- round(cor(x, y, use = "complete.obs"), digits=2)
    p <- cor.test(x, y, use = "complete.obs", alternative = "two.sided")$p.value < 0.05
  } else {
    r = NA
    p = F
  }
  #txt <- paste0("r = ", r)
  if(p){
    txt <- paste(r,expression("*"),sep="")
  } else {
    txt <- r
  }
  
  # PANEL BACKGROUND COLOR
  u <- par('usr') 
  names(u) <- c("xleft", "xright", "ybottom", "ytop")
  bgcolor = cor_colors[2+(-r+1)*100,2]    # converts correlation into a specific color
  do.call(rect, c(col = bgcolor, as.list(u))) # colors the correlation box
  
  # TEXT
  #cex.cor <- 10 / strwidth(txt)
  cex.cor <- 2
  text(0.5, 0.5, txt, cex = cex.cor)
  
}
# Customize upper panel
upper.panel<-function(x, y){
  points(x,y,
         pch = my_shapes[drug_ranks$CD_drug],
         cex = my_cex[drug_ranks$CD_drug],
         col = alpha(my_cols[drug_ranks$CD_drug], my_alpha[drug_ranks$CD_drug])
         )
  
  if(length(intersect(which(!is.na(x)), which(!is.na(y))))>5){
    abline(lm(y ~ x, na.action = "na.omit"), col = "blue", lwd = 2)
    
    newx <- seq(min(x, na.rm = T), max(x, na.rm = T), by = 1)
    preds <- predict(object = lm(y ~ x), interval = "confidence", newdata = data.frame(x = newx))
    polygon(c(rev(newx), newx), c(rev(preds[,3]), preds[,2]), col = "blue", density = 15, border = NA, lwd = 0.3)
    lines(newx, y = preds[ ,3], lty = 'dashed', col = "blue", lwd = 1)
    lines(newx, preds[ ,2], lty = 'dashed', col = "blue", lwd = 1)
  }
}


# Create the plot
pdf(file = paste(in_path, "Correlation_between_drug_ranks_all_individual_patients.pdf", sep=""), width = 10, height = 10)
pairs(x = drug_ranks[,1:11],
      lower.panel = panel.cor,
      upper.panel = upper.panel)
dev.off()

png(filename = paste(in_path, "Correlation_between_drug_ranks_all_individual_patients.png", sep=""), width = 25, height = 25, res = 600, units = "cm")
pairs(x = drug_ranks[,1:11],
      lower.panel = panel.cor,
      upper.panel = upper.panel)
dev.off()

write.table(drug_ranks, file = paste(in_path, "Correlation_between_drug_ranks_all_individual_patients.txt",sep=""), sep="\t", col.names = NA, row.names = T)


