GetError <- function(data, reference){
  deviation <- abs(data - reference)
  percentage <- deviation / reference
  avg_deviation_percentage <- mean(percentage)
}
