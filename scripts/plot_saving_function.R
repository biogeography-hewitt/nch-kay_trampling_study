ggsave_many <- function(filename,plot,...){
  filename <- file.path("figure",c("jpg","png","svg") ,paste0(filename,c(".jpg",".png",".svg")) )
  for (i in filename)ggsave(i,plot = plot,...)
  return(" -- Done -- ")
}