myggplot <- function(res_sizes) {
  p=res_sizes %>%
    ggplot(aes(x = ratio, y = effect)) +
    geom_density_ridges2(fill="orange")+scale_x_continuous(limits = c(-.25,1.25))
  p=p + geom_dotplot(binaxis='y', stackdir='center', dotsize=1)
  # # Box plot with jittered points
  # # 0.2 : degree of jitter in x direction
  p=p + geom_jitter(shape=16, position=position_jitter(0,0.01))
  p=p+ggtitle(paste("rad=",dimnames(ratio)[[1]][i_rad],"thr=",dimnames(ratio)[[2]][i_thr],sep=" "))+scale_y_discrete(limits = rev(levels(res_sizes$effect)))
  p
}