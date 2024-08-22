newScales <- function(moderator) {
  # determine number of levels
  levels_count <- length(levels(moderator))
  # create Diagonal Matrix
  newscales <- diag(levels_count)
  # set last cells to 4 (RoB =4= low risk studies)
  newscales[, levels_count] <- 4
  # convert the matrix to a list of vectors
  newscales_list <- split(newscales, row(newscales))
  return(newscales_list)
}
