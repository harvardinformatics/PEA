MSnProt.ms <- readMSnSet2("ProteinNormalization_matrix.csv", ecol=1)
MSnProt.ms$featureNames <- rep(f)
print(MSnProt.ms)


tiff(paste("Figures/NormalizationBoxPlot.tiff"), width = 4, height = 4, units = 'in', res=600)
.plot(resmir5a6vsn.mss)
dev.off()

pd <- phenoData(resmir5a6vsn.mss)$TreatmentGroup
names(pd) <- sampleNames(resmir5a6vsn.mss)


plotPCA_sc_v2 <- function(m, pdat, component, title='') { # select components
  ## component: either 1 (comp1vscomp2) or 2 (comp2vscomp3)
  pca <- prcomp(m, scale=FALSE)
  df <- as.data.frame(pca$rotation[, 1:4])
  df <- namerows(df, col.name='Samples')
  
  spl <- df$Samples
  cl <- pdat[match(spl, names(pdat))]
  spl <- ifelse(cl==1, pcaCtl, pcaTreat)
  df$Samples <- spl
  
  if (component=='1') {
    p <- ggplot(df, aes(PC1, PC2)) + geom_point(shape=21, size = 2, stroke=1, color="black", aes(fill=Samples)) + scale_fill_manual(values=c("white", "black"))
  } else if (component=='2') {
    p <- ggplot(df, aes(PC2, PC3, colour=Samples)) + geom_point(size=2)
    #p <- ggplot(df, aes(PC3, PC4, colour=Samples)) + geom_point(size=2)
  }
  
  p <- p + theme(legend.position='right', legend.title=element_blank())
  p <- p + labs(title=title)
  
  return(p)
}

#change file name
e <- exprs(resmir5a6vsn.mss)
p <- plotPCA_sc_v2(e, pd, '1', title=paste('', '')) +
  theme_classic()
tiff(paste("Figures/PCAplot.tiff"), width = 4, height = 4, units = 'in', res = 600)

plot(p)
dev.off()


group <- factor(phenoData(resmir5a6vsn.mss)$TreatmentGroup)
design <- model.matrix(~0+group)
colnames(design) <- c('Ctrl', 'Transgn')
fit <- lmFit(e, design)
cm <- makeContrasts(Ctrl-Transgn, levels=design)

fit2 <- contrasts.fit(fit, cm)
fit2 <- eBayes(fit2)




ttUp.df <- topTable(fit2, number=Inf, sort.by ='p', p.value=1)[, c(1, 4, 5)]
#write.table(ttUp.df, "up_fit2ebayes_shinyapp_matrix.csv", sep=",")
ttUp.df$symbol <- unlist(mget(rownames(ttUp.df), uniprotmir5a62sym, ifnotfound=rownames(ttUp.df)))
ttUp.df$FC <- ifelse(ttUp.df$logFC >= 0, inv.glog2(ttUp.df$logFC), -inv.glog2(-ttUp.df$logFC))
ttUp.df$logPval <- -log10(ttUp.df[,c(2)])
ttUp.df <- ttUp.df[which(ttUp.df$logFC  >= 0.58 & ttUp.df$logPval >= 1.3),]
write.table(ttUp.df, file=paste("Figures/StatsUpregulated.csv"), quote=FALSE, sep=',', row.names = FALSE)
rm(ttUp.df)

ttDown.df <- topTable(fit2, number=Inf, sort.by ='p', p.value=1)[, c(1, 4, 5)]
ttDown.df$symbol <- unlist(mget(rownames(ttDown.df), uniprotmir5a62sym, ifnotfound=rownames(ttDown.df)))
ttDown.df$FC <- ifelse(ttDown.df$logFC >= 0, inv.glog2(ttDown.df$logFC), -inv.glog2(-ttDown.df$logFC))
ttDown.df$logPval <- -log10(ttDown.df[,c(2)])
ttDown.df <- ttDown.df[which(ttDown.df$logFC  <= -0.58 & ttDown.df$logPval >= 1.3),]
write.table(ttDown.df, file=paste("Figures/StatsDownregulated.csv"), quote=FALSE, sep=',', row.names = FALSE)
rm(ttDown.df)

tt.df <- topTable(fit2, number=Inf, sort.by ='p', p.value=1)[, c(1, 4, 5)]
tt.df$symbol <- unlist(mget(rownames(tt.df), uniprotmir5a62sym, ifnotfound=rownames(tt.df)))
tt.df$FC <- ifelse(tt.df$logFC >= 0, inv.glog2(tt.df$logFC), -inv.glog2(-tt.df$logFC))
tt.df$logPval <- -log10(tt.df[,c(2)])

write.table(tt.df, file=paste("Figures/StatsTable.csv"), quote=FALSE, sep=',', row.names = FALSE)

 
highlight_df <- dataFilter() %>%
  filter(symbol==input$protint)
highlight_df_down <- dataFilter() %>%
  filter(logFC<=-1)
highlight_df_up <- dataFilter() %>%
  filter(logFC>=1)
ggplot(dataFilter(),aes(x=logFC,y=logPval)) + geom_point(size=2, alpha=1, col='black') +
  labs(title=input$plottitle, x=input$xaxis, y=input$yaxis) +
  theme_update(plot.title=element_text(hjust=0.5), legend.position='none') +
  geom_point(data=dataFilter(), stat='identity', aes(colour=cut(logFC, c(-Inf,-1,1,5))), size=1) + geom_hline(yintercept=-log(0.05,10), linetype="3313", colour="grey") + geom_vline(xintercept=1, linetype="3313", colour="grey") + geom_vline(xintercept=-1, linetype="3313", colour="grey") +
  scale_color_manual(name = "logFC",
                     values = c("(-Inf,-1]" = "blue",
                                "(-1,1]" = "gray",
                                "(1,5]" = "red"),
                     labels = c("decreased", "insignificant", "increased")) +
  geom_point(data=highlight_df, aes(x=logFC,y=logPval), color='green',size=2,alpha=1, col='black') +
  theme_classic()
