# Used for primary drug candidate selection based on dc < 1 and zc < -1.64 (that is P < 0.05 for network proximity)
#' Selection of drug candidates that fulfill network distance criteria dc < 1 and zc < -1.64 (that is P < 0.05 for network distance compared to random)
#'
#' @param input Matrix. Column 4 = dc and column 8 = P value.
#'
#' @return Returns filtered matrix (excluding candidates that do not fulfill criteria). If all candidates are excluded this function returns NULL.
#' @keywords internal
#'
drug_candidate_selection <- function(input){
  if(sum(as.numeric(input[,8])<0.05) > 0){
    if(sum(as.numeric(input[,8])<0.05) > 1){
      input <- input[as.numeric(input[,8])<0.05,] # 1 drug significantly closer than random drug to random disease
    } else {
      input <- input[as.numeric(input[,8])<0.05,] # several drugs significantly closer than random drug to random disease
      input <- matrix(input, nrow = 1)
    }
  } else {
    #print("NO DRUGS matching additional filtering criteria: P < 0.05")
    return(NULL)
  }

  if(sum(as.numeric(input[,4])<1) > 0){
    if(sum(as.numeric(input[,4])<1) > 1){
      input <- input[as.numeric(input[,4])<1,] # 1 drug significantly closer than random drug to random disease
    } else {
      input <- input[as.numeric(input[,4])<1,] # sveral drugs significantly closer than random drug to random disease
      input <- matrix(input, nrow = 1)
    }
  } else {
    #print("NO DRUGS matching additional filtering criteria: d(c) < 1")
    return(NULL)
  }
}
