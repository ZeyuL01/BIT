#' display_tables
#' @description To show the ranking table by inspecting the results of Gibbs sampler.
#' @param file_path path to the saved BIT Gibbs sampling results.
#' @param output_path path to save the rank table.
#' @param burnin number of samples used for burn-in. If not specify, BIT will use the half of the iterations as burn in.
#' @return a data.frame object contains TR names, theta_i, BIT scores and 95 CIs for each TR.
#' @export
display_tables<-function(file_path, output_path, burnin=NULL){
  dat<-readRDS(file_path)

  TR_names <- dat[["TR_names"]]
  theta_i_mat<-dat$theta_i

  if(is.null(burnin)){
    burnin=dim(theta_i_mat)[2]%/%2
  }

  theta_i_mat<-theta_i_mat[!duplicated(TR_names),(dim(theta_i_mat)[2]-burnin):dim(theta_i_mat)[2]]
  TR_names<-TR_names[!duplicated(TR_names)]
  tr_results_i<-rowMeans(theta_i_mat)
  BIT_score<-logistic(tr_results_i)
  CI_intervals<-t(apply(theta_i_mat, 1, function(x) quantile(x, probs = c(0.025, 0.975))))
  results_theta_i=data.frame(TR=TR_names,Theta_i=tr_results_i,lower=CI_intervals[,1],upper=CI_intervals[,2],
                             BIT_score=BIT_score,BIT_score_lower=logistic(CI_intervals[,1]),BIT_score_upper=logistic(CI_intervals[,2]))
  results_theta_i=results_theta_i[order(-results_theta_i$Theta_i),]
  results_theta_i$Rank=rank(-results_theta_i$Theta_i)
  row.names(results_theta_i) <- NULL

  write.csv(results_theta_i,paste0(tools::file_path_as_absolute(output_path),"/",
                                   tools::file_path_sans_ext(basename(file_path)),"_rank_table.csv"))

  return(results_theta_i)
}

#' rank_plot
#' @description To draw a barplot for the top n TRs.
#' @param file_path path to the saved BIT Gibbs sampling results.
#' @param output_path path to save the barplot.
#' @param burnin number of samples used for burn-in. If not specify, BIT will use the half of the iterations as burn in.
#' @param n top n TRs will show in the barplot, default: 10.
#' @param colors colors for each bar, default "NPG" for n<=10, has to be manually specified if n>10.
#' @param main main title for the barplot, default: NULL.
#' @param xlab x axis label, default: BIT score.
#' @param ylab y axis label, default: TR symbols.
#' @return a data.frame object contains TR names, theta_i, BIT scores and 95 CIs for each TR.
#' @export
rank_plot<-function(file_path=NULL, output_path, burnin=NULL, n=10, colors="NPG", main=NULL, xlab="BIT score", ylab="TR symbols"){
  dat<-readRDS(file_path)

  TR_names <- dat[["TR_names"]]
  theta_i_mat<-dat$theta_i

  if(is.null(burnin)){
    burnin=dim(theta_i_mat)[2]%/%2
  }

  theta_i_mat<-theta_i_mat[!duplicated(TR_names),(dim(theta_i_mat)[2]-burnin):dim(theta_i_mat)[2]]
  TR_names<-TR_names[!duplicated(TR_names)]
  tr_results_i<-rowMeans(theta_i_mat)
  BIT_score<-logistic(tr_results)
  CI_intervals<-t(apply(theta_i_mat, 1, function(x) quantile(x, probs = c(0.025, 0.975))))

  results_theta_i=data.frame(TR=TR_names,Theta_i=tr_results_i,BIT_score=BIT_score,lower=CI_intervals[,1],upper=CI_intervals[,2])
  results_theta_i=results_theta_i[order(-results_theta_i$Theta_i),]
  results_theta_i$Rank=rank(-results_theta_i$Theta_i)
  row.names(results_theta_i) <- NULL

  output_file_path_complete<-paste0(tools::file_path_as_absolute(output_path),"/",
                                      tools::file_path_sans_ext(basename(file_path)),".pdf")

  TR_names_used<-results_theta_i[n:1,"TR"]
  BIT_score_used<-results_theta_i[n:1,"BIT_score"]
  BIT_score_lower_used<-results_theta_i[n:1,"BIT_score_lower"]
  BIT_score_upper_used<-results_theta_i[n:1,"BIT_score_upper"]

  if(colors=="NPG"){
    colors=ggsci::pal_npg()(n)
  }

  pdf(output_file_path_complete)
  par(mfrow=c(1,1),oma = c(1,1,1,1) + 0.1,mar = c(4,6.5,2,1) + 0.1,mgp=c(3, 0.5, 0),font.axis=2)
  bp<-barplot(BIT_score_used,horiz=TRUE,yaxt="n",cex.axis=1.5,xlim=c(0,0.13),col=colors,tcl=-0.2,main=main,cex.main=1.3,font.axis=1)
  text(BIT_score_used-max(BIT_score_used)/10,bp,round(BIT_score_used,3),font=1,cex=1.2)
  points(BIT_score_used,bp,pch=16)
  segments(BIT_score_used,bp,BIT_score_upper_used,bp,col="black",lwd=2)
  segwidth<-(bp[2]-bp[1])/4
  segments(BIT_score_upper_used,bp+(bp[2]-bp[1])/4,BIT_score_upper_used,bp-(bp[2]-bp[1])/4,lwd=2)

  axis(side=2, las=1, at=bp, labels=TR_names_used,tcl=-0.2,cex.axis=1.2,font.axis=1)
  title(xlab=xlab,line = 3, cex.lab=1.4,font.lab=2)
  title(ylab=ylab,line = 5.5, cex.lab=1.4,font.lab=2)
  title(main=main,cex.main=1.4,font.main=2)
  box()
  dev.off()

  return()
}

#' rank_plot
#' @description To draw a barplot for the top n TRs.
#' @param file1_path path to the saved BIT Gibbs sampling results of input 1.
#' @param file2_path path to the saved BIT Gibbs sampling results of input 2.
#' @param output_path path to save the barplot.
#' @param burnin number of samples used for burn-in. If not specify, BIT will use the half of the iterations as burn in.
#' @return NULL
#' @export
compare_scatter_plot<-function(file1_path, file2_path, output_path, burnin=NULL){
  dat1<-readRDS(file1_path)
  dat2<-readRDS(file2_path)

  TR_names_dat1 <- dat1[["TR_names"]]
  dat1_theta_i_mat<-dat1$theta_i

  TR_names_dat2 <- dat2[["TR_names"]]
  dat2_theta_i_mat<-dat2$theta_i

  if(is.null(burnin)){
    burnin=dim(dat1_theta_i_mat)[2]%/%2
  }

  output_file_path_complete<-paste0(tools::file_path_as_absolute(output_path),"/",
                                    tools::file_path_sans_ext(basename(file1_path)),"_",
                                    tools::file_path_sans_ext(basename(file2_path)),"_compare.pdf")

  theta_i_mat_dat1<-dat1_theta_i_mat[!duplicated(TR_names_dat1),(dim(dat1_theta_i_mat)[2]-burnin):dim(dat1_theta_i_mat)[2]]
  TR_names_dat1<-TR_names_dat1[!duplicated(TR_names_dat1)]
  tr_results_i_dat1<-rowMeans(theta_i_mat_dat1)
  BIT_score<-logistic(tr_results_i_dat1)
  CI_intervals<-t(apply(theta_i_mat_dat1, 1, function(x) quantile(x, probs = c(0.025, 0.975))))

  results_theta_i_dat1=data.frame(TR=TR_names_dat1,Theta_i=tr_results_i_dat1,BIT_score=BIT_score,lower=CI_intervals[,1],upper=CI_intervals[,2])
  results_theta_i_dat1=results_theta_i_dat1[order(-results_theta_i_dat1$Theta_i),]
  results_theta_i_dat1$Rank=rank(-results_theta_i_dat1$Theta_i)
  row.names(results_theta_i_dat1) <- NULL

  theta_i_mat_dat2<-dat2_theta_i_mat[!duplicated(TR_names_dat2),(dim(dat2_theta_i_mat)[2]-burnin):dim(dat2_theta_i_mat)[2]]
  TR_names_dat2<-TR_names_dat2[!duplicated(TR_names_dat2)]
  tr_results_i_dat2<-rowMeans(theta_i_mat_dat2)
  BIT_score<-logistic(tr_results_i_dat2)
  CI_intervals<-t(apply(theta_i_mat_dat2, 1, function(x) quantile(x, probs = c(0.025, 0.975))))

  results_theta_i_dat2=data.frame(TR=TR_names_dat2,Theta_i=tr_results_i_dat2,BIT_score=BIT_score,lower=CI_intervals[,1],upper=CI_intervals[,2])
  results_theta_i_dat2=results_theta_i_dat2[order(-results_theta_i_dat2$Theta_i),]
  results_theta_i_dat2$Rank=rank(-results_theta_i_dat2$Theta_i)
  row.names(results_theta_i_dat2) <- NULL

  file1_values<-results_theta_i_dat1[,"BIT_score"]
  file2_values<-results_theta_i_dat2[match(results_theta_i_dat1[,"TR"],results_theta_i_dat2[,"TR"]),"BIT_score"]

  aligned_tables<-data.frame(TR=results_theta_i_dat1[,"TR"],file1_values=file1_values,file2_values=file2_values)

  TR_names1<-results_theta_i_dat1[1:10,"TR"]
  TR_names2<-results_theta_i_dat2[1:10,"TR"]
  TR_union<-union(TR_names1,TR_names2)

  file1_coords<-aligned_tables[which(aligned_tables[,"TR"]%in%TR_union),"file1_values"]
  file2_coords<-aligned_tables[which(aligned_tables[,"TR"]%in%TR_union),"file2_values"]
  TR_plot_names<-aligned_tables[which(aligned_tables[,"TR"]%in%TR_union),"TR"]

  scale<-c()
  scale1_vec<-(file1_values/max(file1_values))
  scale2_vec<-(file2_values/max(file2_values))
  for(i in 1:length(file1_values)){
    scale<-c(scale,max(c(scale1_vec[i],scale2_vec[i])))
  }

  size<-1-median(scale)+scale

  colors<-rgb(file1_values/max(file1_values),0.4,file2_values/max(file2_values),scale)

  pdf(output_file_path_complete)
  par(mfrow=c(1,1),oma = c(1,1,1,1) + 0.1,mar = c(4,4,2,1) + 0.1,mgp=c(3, 0.5, 0),font.axis=2)
  plot(file1_values,file2_values,type="p",pch=19,col=colors,cex=size, xlab="Region set 1 BIT score",ylab="Region set 2 BIT score",font.lab=2,cex.lab=1.4)
  basicPlotteR::addTextLabels(file1_coords, file2_coords, TR_plot_names, cex.label=1.4,lwd=2)
  box()
  dev.off()

  return()
}

logistic <- function(x) {
  1 / (1 + exp(-x))
}
