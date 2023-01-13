#
# Parse XML files
# DrugBank XML to table
############################
"
XML levels:
doc = document
 root = rootnode
  child = child
   subn = subnode
  child
 root
doc

XML functions for a given node:

Function Description
xmlName() name of the node
xmlSize() number of subnodes
xmlAttrs() named character vector of all attributes
xmlGetAttr() value of a single attribute
xmlValue() contents of a leaf node
xmlParent() name of parent node
xmlAncestors() name of ancestor nodes
getSibling() siblings to the right or to the left
xmlNamespace() the namespace (if there?s one)
"

extract_drug_drug_interaction_information <- function(cores = 3/4, target_info = F, ddi_info = F){
  
  setwd("~/Desktop/Barabasi project/Input/DrugBank")
  
  library(XML)
  library(methods)
  library(doParallel)
  
  # for parallel running
  x <- detectCores()*cores
  registerDoParallel(cores = x)
  rm(x)
  
  # load DrugBank
  drugBank <<- data <- xmlParse(file = "full database.xml")
  
  rootnode <- xmlRoot(data)
  rootsize <- xmlSize(rootnode)
  
  # Builds table with 6 columns (drugID, drug name, approved/not approved, targetID, target name, target organism)
  if(target_info == T){
    out <- foreach(d = c(1:rootsize), .combine = rbind) %dopar% { # for every drug
      
      # extract drugBank ID
      t1 <- xmlValue(getElement(rootnode[[d]], name = "drugbank-id"))
      
      # extract drug name
      t2 <- xmlValue(getElement(rootnode[[d]], name = "name"))
      
      # extract group (approved or not)
      t3 <- getElement(rootnode[[d]], name = "groups")
      x <- xmlSize(t3)
      if(x > 1){
        temp <- vector()
        for(i in c(1:x)){
          temp[i] <- xmlValue(getElement(t3[i], name = "group"))
        }
        temp <- unique(as.character(temp))
        t3 <- temp[1]
        for(i in c(2:length(temp))){
          t3 <- paste(t3, temp[i], sep = ", ")  
        }
        rm(temp, i)
      } else if (x==1) {
        t3 <- xmlValue(t3)
      } else { # x == 0 -> information not existent for drug
        t3 <- NA
      }
      rm(x)
      
      # extract targets (id, name, gene-name, action and organism)
      temp <- getElement(rootnode[[d]], name = "targets")
      x <- xmlSize(temp)
      if(x>1){
        t4 <- vector()
        t5 <- vector()
        t6 <- vector()
        t7 <- vector()
        t8 <- vector()
        for(i in c(1:x)){
          temp2 <- getElement(temp[i], name = "target")
          t8[i] <- xmlValue(getElement(temp2, name = "organism"))
          temp3 <- getElement(temp2, name = "actions")
          if(xmlSize(temp3)>1){
            t7[i] <- xmlValue(getElement(temp3[1], name = "action"))
            for(j in c(1:(xmlSize(temp3)-1))){
              t7[i] <- paste(t7[i], xmlValue(getElement(temp3[j+1], name = "action")), sep =", ")
            }
          } else if (xmlSize(temp3)==1){
            t7[i] <- xmlValue(temp3)
          } else if (xmlSize(temp3)==0){
            t7[i] <- NA
          }

          temp2 <- getElement(temp2, name = "polypeptide")
          if(xmlSize(temp2)>0){
            temp3 <- xmlGetAttr(temp2, name = "id")
            if(is.null(temp3)){
              t4[i] <- NA
            } else {
              t4[i] <- temp3
            }
            t5[i] <- xmlValue(getElement(temp2, name = "name"))
            t6[i] <- xmlValue(getElement(temp2, name = "gene-name"))
          } else {
            t4[i] <- NA
            t5[i] <- NA
            t6[i] <- NA
          }
        }
        rm(temp2, temp3, i)
      } else if(x == 1){
        temp <- getElement(temp, name = "target")
        t8 <- xmlValue(getElement(temp, name = "organism"))
        temp3 <- getElement(temp, name = "actions")
        if(xmlSize(temp3)>1){
          t7 <- xmlValue(getElement(temp3[1], name = "action"))
          for(i in c(1:(xmlSize(temp3)-1))){
            t7 <- paste(t7, xmlValue(getElement(temp3[i+1], name = "action")), sep =", ")
          }
        } else if (xmlSize(temp3)==1){
          t7 <- xmlValue(temp3)
        } else if (xmlSize(temp3)==0){
          t7 <- NA
        }
      
        temp <- getElement(temp, name = "polypeptide")
        if(xmlSize(temp)>0){
          t4 <- xmlGetAttr(temp, name = "id")
          if(is.null(t4)){
            t4 <- NA
          }
          t5 <- xmlValue(getElement(temp, name = "name"))
          t6 <- xmlValue(getElement(temp, name = "gene-name"))
          if(is.null(t6)){
            t6 <- NA
          }
        } else {
          t4 <- NA
          t5 <- NA
          t6 <- NA
        }
      } else if(x == 0) { # x == 0 -> target information not existent for drug
        t4 <- NA
        t5 <- NA
        t6 <- NA
        t7 <- NA
        t8 <- NA
      }
      rm(temp, x)
      
      x <- length(t4)
      temp <- cbind(rep(t1, times = x), rep(t2, times = x), rep(t3, times = x), t4, t5, t6, t7, t8)
      
      temp
    } # end foreach
    colnames(out) <- c("drugID", "drug_name", "status", "Swiss_target_ID", "target_name", "gene_symbol", "target_organism", "drug_action")
    drugBank <<- out
  }
  
  # Builds table with X columns (drugID, drug name, interacting drugID, interacting drug name, interacting drug description)
  if(ddi_info == T){
    out <- foreach(d = c(1:rootsize), .combine = rbind) %dopar% { # for every drug
      # extract drugBank ID
      t1 <- xmlValue(getElement(rootnode[[d]], name = "drugbank-id"))
      
      # extract drug name
      t2 <- xmlValue(getElement(rootnode[[d]], name = "name"))
      
      # extract DDI (drug id, name and interaction description)
      temp <- getElement(rootnode[[d]], name = "drug-interactions")
      x <- xmlSize(temp)
      if(x>1){
        t3 <- vector()
        t4 <- vector()
        t5 <- vector()
        for(i in c(1:x)){
          temp2 <- getElement(temp[i], name = "drug-interaction")
          t3[i] <- xmlValue(getElement(temp2, name = "drugbank-id"))
          t4[i] <- xmlValue(getElement(temp2, name = "name"))
          t5[i] <- xmlValue(getElement(temp2, name = "description"))
        }
        rm(temp2,i)
      } else if(x == 1) {
        temp <- getElement(temp, name = "drug-interaction")
        t3 <- xmlValue(getElement(temp, name = "drugbank-id"))
        t4 <- xmlValue(getElement(temp, name = "name"))
        t5 <- xmlValue(getElement(temp, name = "description"))
      } else { # x = 0 -> DDI not existent for this drug 
        t3 <- NA
        t4 <- NA
        t5 <- NA
      }
      rm(x, temp)
      
      x <- length(t3)
      temp <- cbind(rep(t1, times = x), rep(t2, times = x), t3, t4, t5)
      
      temp # output
    } # end foreach
    colnames(out) <- c("drugID_1", "drug_name_1", "drugID_2", "drug_name_2", "interaction")
    write.table(out, file = "DDI_info_drugBank.txt", sep = "\t", col.names = T, row.names = F)
    DDI_info <<- out
  }

}



out <- extract_drug_drug_interaction_information(target_info = T)
write.table(out, file = "Input/Drugs/all_drug_targets_drug_bank.txt", sep="\t", col.names = T, row.names = F)


