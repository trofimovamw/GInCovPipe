#library('latex2exp')
library(ggplot2)
library(grid)
source("dateFormatRoutines.r")

plotSplineDerivative<- function(spline.table, outputpath) {
  p_spline_d1 <- ggplot(data=spline.table, aes(x=t, y=value)) +
    #geom_point(aes(tt.df[,1], d1_gam), alpha=0.5, size=2) +
    geom_line(color="red", size = 2, alpha=0.7)+
    xlab("") +
    ylab(expression(paste("derivative of est.", theta))) +
    theme(axis.text = element_text(size=14),
          axis.title = element_text(size=14,face="bold"),
          legend.text = element_text(size=12),
          legend.title = element_text(size=12, face="bold"))

  outputFile = paste0(outputpath, "/esti_theta_derivSplineFit.pdf")
  ggsave(p_spline_d1,
         height = 8,
         width = 16,
         dpi = 220,
         device = "pdf",
         file = outputFile)
}

plotRatio <- function(ratios, outputFile) {
  data <- data.frame(t=seq(from = 0, to = length(ratios)-1,by=1),x=ratios)
  plot_r <- ggplot(data=data, aes(x=t,y=ratios)) +
    geom_line(color="blue", size=2, alpha=0.7)+
    xlab("") +
    ylab("ratio")
  ggsave(plot_r,
         height = 8,
         width = 16,
         dpi = 220,
         device = "pdf",
         file = outputFile)
}

# plotRatioThetaVsCases <- function(spline.table){
#   minX <- min(spline.table$value_trueN)
#   maxX <- max(spline.table$value_trueN)
#   minY <- min(spline.table$value)
#   maxY <- max(spline.table$value)
#   minDate <-
#
#   casesRescaled <- rescale(spline.table$value_trueN, minX, maxX, minY, maxY)
#
#   ratio <- spline.table$value/casesRescaled
#
#   theta_01 <- (spline.table$value- min(spline.table$value))/max(spline.table$value -  min(spline.table$value)) +0.5
#   cases_01 <- (spline.table$value_trueN- min(spline.table$value_trueN))/max(spline.table$value_trueN -  min(spline.table$value_trueN)) +0.5
#
#   ratio = (theta_01 -cases_01)/theta_01
#   ratio <- theta_01/cases_01
#   ggplot() +
#     geom_line(aes(x=day.as.Date(spline.table$t, )[ratio != -Inf], y=ratio[ratio != -Inf]))
#
#
# }

rescale <- function(x, minX, maxX, minY, maxY) {
  return((maxY-minY)*(x - minX)/(maxX-minX) + minY)
}

plotSplineWithNewCases <-function(data.table, meta.table, spline.table, date_m, outputFile) {

  minX <- max(min(spline.table$value_trueN),0)
  maxX <- max(spline.table$value_trueN)
  minY <- max(min(spline.table$value),0)
  maxY <- max(spline.table$value)
  minDate <- min(as.Date(data.table$meanBinDate))
  maxDate <- max(as.Date(data.table$meanBinDate))
  mycolors <- c("estid"="darkred", "trued"="darkblue","esti"="red", "true"="blue")
  rel_vs_true_ratio = max(spline.table$value)/max(spline.table$value_trueN)
  p_spline_esti_realN <- ggplot() +
    #geom_point(aes(doy.as.Date(data.table$doy), data.table$trueN*rel_vs_true_ratio), size=2, color=mycolors["trued"], alpha=0.5)+
    geom_point(aes(days.as.Date(data.table$t, minDate), rescale(data.table$trueN, minX, maxX, minY, maxY)),size=2, color=mycolors["trued"], alpha=0.5)+
    geom_point(aes(days.as.Date(data.table$t, minDate), data.table$value), size=2, color=mycolors["estid"], alpha=0.5)+
    #geom_line(aes(x=doy.as.Date(spline.table$doy), y=spline.table$value_trueN*rel_vs_true_ratio), color=mycolors["true"], size = 2, alpha=0.5)+
    geom_line(aes(x=days.as.Date(spline.table$t, minDate), y=rescale(spline.table$value_trueN, minX, maxX, minY, maxY)), color=mycolors["true"], size = 2, alpha=0.5)+
    geom_line(aes(x=days.as.Date(spline.table$t, minDate), y=spline.table$value), color=mycolors["esti"], size = 2, alpha=0.5)+
    geom_vline(xintercept=as.Date(date_m), linetype = "dashed", color="darkgray") +
    scale_x_date(date_breaks = "months" , date_labels = "%B %Y") +
    scale_y_continuous(
      expression(paste("est. ", theta)),
      #sec.axis = sec_axis(~ . * 1/rel_vs_true_ratio, name = "new cases")
      sec.axis = sec_axis(~ rescale(., minY, maxY, minX, maxX ), name = "active cases")
    ) +
    xlab("")+
    theme(
      axis.text.x=element_text(angle = 0, vjust = 1, hjust=0),
      axis.title.y = element_text(color = mycolors["estid"]),
      axis.text.y = element_text(color = mycolors["estid"]),
      axis.title.y.right = element_text(color = mycolors["trued"]),
      axis.text.y.right = element_text(color = mycolors["trued"]),
      axis.text = element_text(size=14),
      axis.title = element_text(size=15),
      legend.text = element_text(size=14),
      legend.title = element_text(size=15, face="bold")
    )
    p_sample_size <- ggplot() +
    stat_ecdf(aes(days.as.Date(data.table$t, minDate)), geom = "step", size=0.7, color="darkgray")+
    stat_ecdf(aes(days.as.Date(data.table$t, minDate)), geom = "step", size=0.7, color="darkgray", alpha="0.0") +
    geom_vline(xintercept=as.Date(date_m), linetype = "dashed", color="darkgray") +
    scale_y_continuous(
      expression("cumulative frequency on date"),
      sec.axis = sec_axis(~ . , name = derive())) +
    xlab("") +
    xlim(minDate, maxDate) +
    scale_x_date(date_breaks = "months" , date_labels = "%B %Y") +
    theme(
      axis.text.x=element_text(angle = 0, vjust = 1, hjust=0),
      axis.title.y = element_text(color = "black"),
      axis.text.y = element_text(color = "black"),
      axis.title.y.right = element_text(color = "black"),
      axis.text.y.right = element_text(color = "black"),
      axis.text = element_text(size=14),
      axis.title = element_text(size=15),
      legend.text = element_text(size=14),
      legend.title = element_text(size=15, face="bold")
    )

    pdf(outputFile, height = 15, width = 13, paper = "special")
    grid.draw(rbind(ggplotGrob(p_spline_esti_realN), ggplotGrob(p_sample_size), size = "last"))
    dev.off()

}

plotCummulativeSampleSize <- function(data.table,outputFile) {


  minDate <- min(as.Date(data.table$meanBinDate))
  mycolors <- c("estid"="darkred", "trued"="darkblue","esti"="red", "true"="blue")
  p_sample_size <- ggplot() +
    #geom_point(aes(doy.as.Date(data.table$doy), data.table$trueN*rel_vs_true_ratio), size=2, color=mycolors["trued"], alpha=0.5)+
    geom_line(aes(days.as.Date(data.table$t, minDate), cumsum(data.table$sampleSize)),size=2, color=mycolors["trued"], alpha=0.5)+
    #geom_line(aes(x=doy.as.Date(spline.table$doy), y=spline.table$value_trueN*rel_vs_true_ratio), color=mycolors["true"], size = 2, alpha=0.5)+
    scale_y_continuous("sample size") + xlab("")+
    theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12, face="bold"))
    #outputFile = paste0(outputpath, "/esti_vs_real_PopulationSize_splineFit.pdf")
    ggsave(p_sample_size,
           height = 8,
           width = 16,
           dpi = 220,
           device = "pdf",
           file = outputFile)

}

plotSpline <- function(data.table, spline.table, outputFile) {
  minDate <- min(as.Date(data.table$meanBinDate))
  p_spline <- ggplot() +
    geom_point(aes(days.as.Date(data.table$t, minDate), y=data.table$value), alpha=0.5, size=2) +
    geom_line(aes(x=days.as.Date(spline.table$t,minDate), spline.table$value), color="red", size = 2, alpha=0.7) +
    xlab("") +
    scale_x_date(date_breaks = "months" , date_labels = "%Y-%m-%d") +
    ylab(expression(paste("estim. ", theta)))+
    theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12, face="bold"))

  ggsave(p_spline,
       height = 8,
       width = 16,
       dpi = 220,
       device = "pdf",
       file = outputFile)
}

plotNucleotideDiversity_perSite <- function(nuclDiversity_perSite_perBin) {
  #[nuclDiversity_perSite_perBin$nuclDiv !=0, ]
  p_pS <- ggplot(nuclDiversity_perSite_perBin, aes(x=date, y=nuclDiv, col=binning, group = interaction(date, binning))) +
    #geom_boxplot(aes(fill=binning), alpha=0.3, outlier.shape=NA, size=0.8)+
    #geom_point(position=position_jitterdodge(),alpha=0.3)+
    #TODO: warum ist das dodging der points pro Datum nicht gleich dem Boxplot??
    geom_pointrange(aes(ymin = nuclDiv - sqrt(variance), ymax = nuclDiv + sqrt(variance), col=binning), size=0.4, alpha = 0.7,
                      #binaxis='y',
                      position=position_jitterdodge(),
                      linetype='solid') +
    #TODO: anzahl der samples mit hilfe von punktgrößen?
    # scale_size(range = c(0.1, 1),
    #          breaks = c(0, 5, 10, 50, 100),
    #          labels = c("<5", "5", "10", "50", "100+"),
    #          guide = "legend")+
    #TODO wie macht man ein expectation value in ggplot (mathbb)
    #scale_y_continuous(trans = 'log10')+
    ylab(parse(text = TeX('$E(D_i)$'))) +
    xlab("")
    #scale_color_manual(values = c("#00AFBB", "#E7B800"))

  return(p_pS)
}


plotNucleotideDiversity_perBin <- function(nuclDiversity_df) {
  p_pB <- ggplot(nuclDiversity_df, aes(x=date, y=nuclDiv, col=binning)) +
    # geom_pointrange(aes(ymin = nuclDiv - variance_sample, ymax = nuclDiv + variance_sample),#, size=N),
    #                 alpha=0.7,
    #                 #position=position_dodge(width = 5),
    #                 linetype='solid') +
    #TODO: ymin = max(0, nc-var) geht irgendwie nicht
    #geom_ribbon(aes(ymin = max(0,nuclDiv - sqrt(variance_sample)), ymax = nuclDiv + sqrt(variance_sample), fill=binning), alpha=0.2)+
    geom_line(size=1, alpha=0.7) +
    geom_point(aes(size=N), alpha=0.7) +
  scale_size(range = c(1, 10),
                      #breaks = c(0, 5, 10, 50, 100),
                      #labels = c("<5", "5", "10", "50", "100+"),
             guide = FALSE) +
    theme(legend.position="bottom") +
    ylab(parse(text = TeX('$\\pi$'))) +
    xlab("")

  return(p_pB)
}

plotChangeInNucleotideDiversity_perBin <- function(nuclDiversity_df) {
  p_ncC <- ggplot(nuclDiversity_df[!is.nan(nuclDiversity_df$delta_nd),], aes(x=date, y=delta_nd, col=binning)) +
    # geom_pointrange(aes(ymin = nuclDiv - variance_sample, ymax = nuclDiv + variance_sample),#, size=N),
    #                 alpha=0.7,
    #                 #position=position_dodge(width = 5),
    #                 linetype='solid') +
    #TODO: ymin = max(0, nc-var) geht irgendwie nicht
    #geom_ribbon(aes(ymin = max(0,nuclDiv - sqrt(variance_sample)), ymax = nuclDiv + sqrt(variance_sample), fill=binning), alpha=0.2)+
    geom_line(size=1, alpha=0.7) +
    geom_point(aes(size=delta_var), alpha=0.7) +
    scale_size(range = c(1, 10),
               #breaks = c(0, 5, 10, 50, 100),
               #labels = c("<5", "5", "10", "50", "100+"),
               guide = FALSE) +
    theme(legend.position="bottom") +
    ylab(parse(text = TeX('$\\frac{\\Delta\\pi}{\\Delta t}$'))) +
    xlab("")

  return(p_ncC)
}

plotChangeInNDvsPopSize_perBin <- function(nuclDiversity_df) {
  p_nd_vs_pop <- ggplot(nuclDiversity_df[!is.nan(nuclDiversity_df$delta_nd),], aes(x=delta_N, y=delta_nd, col=binning)) +
    # geom_pointrange(aes(ymin = nuclDiv - variance_sample, ymax = nuclDiv + variance_sample),#, size=N),
    #                 alpha=0.7,
    #                 #position=position_dodge(width = 5),
    #                 linetype='solid') +
    #TODO: ymin = max(0, nc-var) geht irgendwie nicht
    #geom_ribbon(aes(ymin = max(0,nuclDiv - sqrt(variance_sample)), ymax = nuclDiv + sqrt(variance_sample), fill=binning), alpha=0.2)+
    #geom_line(size=1, alpha=0.7) +
    geom_point(size=2, alpha=0.7) +
    theme(legend.position="bottom") +
    ylab(parse(text = TeX('$\\Delta\\pi$'))) +
    xlab(parse(text = TeX('$\\Delta N_p$')))

  return(p_nd_vs_pop)
}
