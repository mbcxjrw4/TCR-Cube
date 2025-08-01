scatter_one <- function(tumour, geneA, geneB){
    df <- sheet[type=="normal"|tissue.cancer==tumour]
    col.keep <- c("Sample.ID", "type", geneA, geneB)
    df <- df[, ..col.keep]
    data.table::setnames(df, old=c(geneA, geneB), new=c("geneA", "geneB"))
    
    fig <- ggplot(df, aes(x=geneA, y=geneB, colour=type))+
           geom_point(size=0.7, position=position_jitter(h=0.05, w=0.05))+
           scale_colour_manual(name="", values=c("normal"="#5A66AF", "cancer"="#CD1310"))+
           labs(x=paste0(geneA, " log(TPM)"), y=paste0(geneB, " log(TPM)"))+
           ggtitle(tumour)+
           theme_bw()+
           theme(plot.title = element_text(hjust = 0.5),
                 legend.position = "bottom")
    
      return(fig)
}
