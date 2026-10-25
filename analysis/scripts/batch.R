suppressMessages({library(sva);library(limma);library(jsonlite)})
g<-commandArgs(TRUE)[1]
x<-as.matrix(read.delim(gzfile(sprintf("%s_rma.tsv.gz",g)),row.names=1,check.names=FALSE))
B<-read.delim("batches.tsv"); B<-B[B$cohort==g,]; rownames(B)<-B$gsm; B<-B[colnames(x),]
batch<-factor(B$batch); grp<-factor(B$grp)
probes<-c('205476_at','207356_at','212531_at','205916_at','202917_s_at','203535_at','203691_at','41469_at','220322_at','204470_at','202859_x_at','211506_s_at','220054_at','207901_at','1552912_a_at','1561853_a_at','221165_s_at','222974_at','221111_at','208402_at','216876_s_at','234408_at','214218_s_at','221728_x_at','224588_at','224589_at','224590_at','227671_at','243712_at','201909_at','205000_at','206700_s_at')
pcr2<-function(m){v<-apply(m,1,var); top<-order(-v)[1:5000]; p<-prcomp(t(m[top,]),center=TRUE)
  ve<-(p$sdev^2/sum(p$sdev^2))[1:5]
  rb<-sapply(1:5,function(i) summary(lm(p$x[,i]~batch))$r.squared)
  rg<-sapply(1:5,function(i) summary(lm(p$x[,i]~grp))$r.squared)
  list(var_expl=round(ve,3),r2_batch=round(rb,3),r2_group=round(rg,3),weighted_r2_batch=round(sum(ve*rb)/sum(ve),3),weighted_r2_group=round(sum(ve*rg)/sum(ve),3))}
partial<-function(m){X0<-model.matrix(~grp); X1<-model.matrix(~grp+batch)
  rss<-function(X){f<-lm.fit(X,t(m)); colSums(f$residuals^2)}
  tss<-rowSums((m-rowMeans(m))^2); r0<-rss(X0); r1<-rss(X1)
  pr<-(r0-r1)/r0; round(c(median=median(pr),q90=unname(quantile(pr,.9)),frac_gt_0.2=mean(pr>0.2)),3)}
de<-function(m){ref<-if("NL"%in%levels(grp)) "NL" else "HC"; keep<-grp%in%c("LS",ref); gg<-droplevels(grp[keep])
  gg<-relevel(gg,ref=ref); d<-model.matrix(~gg); f<-eBayes(lmFit(m[,keep],d)); tt<-topTable(f,coef=2,number=Inf)
  c(n_fdr05_fc2=sum(tt$adj.P.Val<0.05 & abs(tt$logFC)>=1), n_fdr05=sum(tt$adj.P.Val<0.05))}
res<-list(cohort=g,n=ncol(x),n_batches=nlevels(batch))
mats<-list(none=x)
mats$combat_with_group<-ComBat(x,batch=batch,mod=model.matrix(~grp))
mats$combat_without_group<-ComBat(x,batch=batch,mod=NULL)
for(k in names(mats)){m<-mats[[k]]
  res[[k]]<-list(pca=pcr2(m),partial_r2_batch_given_group=partial(m),de_LS_vs_ref=de(m))
  s<-m[intersect(probes,rownames(m)),]; write.table(round(s,5),sprintf("batch_%s_%s_modprobes.tsv",g,k),sep="\t",quote=FALSE,col.names=NA)
  cat(g,k,"done\n")}
write(toJSON(res,auto_unbox=TRUE,pretty=TRUE,digits=NA),sprintf("batch_%s.json",g))
