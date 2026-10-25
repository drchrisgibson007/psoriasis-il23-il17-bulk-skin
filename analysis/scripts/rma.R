suppressMessages({library(affy);library(affyPLM);library(affyio);library(jsonlite)})
args<-commandArgs(TRUE); g<-args[1]; doqc<-length(args)<2 || args[2]!="noqc"
raw<-file.path(Sys.getenv("RAW_DIR","../raw"),g)
gsms<-readLines(sprintf("../rma/%s_gsms.txt",g))
f<-list.files(raw,pattern="CEL.gz$",full.names=TRUE,ignore.case=TRUE)
gsm<-sub("^(GSM[0-9]+).*","\\1",basename(f))
keep<-gsm %in% gsms; f<-f[keep]; gsm<-gsm[keep]
cat(g,"files used",length(f),"of listed",length(gsms),"\n")
# headers: scan date, cdf
hd<-t(sapply(f,function(x){h<-read.celfile.header(x,info="full");c(cdf=h$cdfName,date=ifelse(is.null(h$ScanDate),NA,h$ScanDate))}))
write.table(data.frame(gsm=gsm,cdf=hd[,1],scandate=hd[,2]),sprintf("../rma/%s_headers.tsv",g),sep="\t",quote=FALSE,row.names=FALSE)
e<-justRMA(filenames=basename(f),celfile.path=raw)
x<-exprs(e); colnames(x)<-gsm
write.table(round(x,5),gzfile(sprintf("../rma/%s_rma.tsv.gz",g)),sep="\t",quote=FALSE,col.names=NA)
rm(e,x);gc()
if(doqc){
  ab<-ReadAffy(filenames=f); sampleNames(ab)<-gsm
  plm<-fitPLM(ab,output.param=list(weights=FALSE,residuals=FALSE,varcov="none",resid.SE=TRUE),verbos=0)
  rm(ab);gc()
  n<-NUSE(plm,type="stats"); r<-RLE(plm,type="stats")
  q<-data.frame(gsm=colnames(n),NUSE_median=n["median",],NUSE_IQR=n["IQR",],RLE_median=r["median",],RLE_IQR=r["IQR",])
  write.table(q,sprintf("../rma/%s_qc.tsv",g),sep="\t",quote=FALSE,row.names=FALSE)
}
cat(g,"done\n")
