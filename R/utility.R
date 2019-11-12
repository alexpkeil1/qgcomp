construction <- function(i="", j=""){
  msg ="This function is not yet fully implemented"
  if(j != "") msg = paste0(msg, ": ", j)
  if(i %in% c("msg", "message", "warning", "wrn", "warn")) warning(msg)
  else stop(msg)
}

# fake cox family function
cox <- function(){
  obj = binomial(link="log")
  obj$family="cox"
  obj
}
