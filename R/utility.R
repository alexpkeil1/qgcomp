construction <- function(i=""){
  msg = "This function is not yet functional"
  if(i %in% c("msg", "message", "warning", "wrn", "warn")) warning(msg)
  else stop(msg)
}
