library(ggplot2)
library(tidyr)
library(reshape)
library(pca3d)
library(wesanderson)
library(stringr)

pn.df <- read.csv("ProteinNormalization_matrix_box.csv", sep=',', header=FALSE)
pn.df <- melt(pn.df)
pn.df['variable'] <- NULL
columnNames <- c('TreatmentGroup', 'Channel', 'ProteinNorm')
colnames(pn.df) <- columnNames

p <- ggplot(data=pn.df, aes(x=ProteinNorm, y=Channel, colour=TreatmentGroup))+geom_boxplot()
tiff(paste("NormalizationBoxPlot.tiff"), width = 4, height = 4, units = 'in', res=200)
plot(p)
dev.off()



#change file name
pn_matrix.df <- read.csv("ProteinNormalization_matrix_transpose.csv", sep=',', row.names = 'Protein', header=TRUE)
pca <- prcomp(pn_matrix.df[1:length(pn_matrix.df)], scale=FALSE)
df <- as.data.frame(pca$rotation[, 1:4])
df <- namerows(df, col.name='Repl')
df <- df %>% separate(Repl, c('Samples', 'Replicates'))

#MetaEDIT.df <- read.table('PCA_matrix.csv', header=TRUE,quote='\"', sep=',', comment.char='')

p <- ggplot(df, aes(PC1, PC2, colour=Samples)) + geom_point(size=4) + scale_color_manual(values=wes_palette(n=3, name="Darjeeling1")) + ggtitle("Protein Normalization, Test")
tiff("PCAPlot.tiff", width = 6, height = 8, units = 'in', res=200)
plot(p)
dev.off()


pn_norm_matrix.df <- read.csv("ProteinNormalization_matrix.csv", sep=',')
group <- factor(pn_norm_matrix.df$TreatmentGroup)
design <- model.matrix(~0+group)
colnames(design) <- c('Ctrl', 'Transgn')
fit <- lmFit(pn_matrix.df, design)
cm <- makeContrasts(Ctrl-Transgn, levels=design)

fit2 <- contrasts.fit(fit, cm)
fit2 <- eBayes(fit2)

uniprot2genename.df <- read.table('Template_uni_gn.csv', header=FALSE, sep=',', quote='')
uniprotmir5a62sym <- new.env(hash=TRUE)
apply(uniprot2genename.df, 1, function(x) {
  x <- as.character(x)
  uniprotmir5a62sym[[x[1]]] <- x[2]
})


ttUp.df <- topTable(fit2, number=Inf, sort.by ='p', p.value=1)[, c(1, 4, 5)]
#write.table(ttUp.df, "up_fit2ebayes_shinyapp_matrix.csv", sep=",")
ttUp.df$symbol <- unlist(mget(rownames(ttUp.df), uniprotmir5a62sym, ifnotfound=rownames(ttUp.df)))
ttUp.df$FC <- ifelse(ttUp.df$logFC >= 0, inv.glog2(ttUp.df$logFC), -inv.glog2(-ttUp.df$logFC))
ttUp.df$logPval <- -log10(ttUp.df[,c(2)])
ttUp.df <- ttUp.df[which(ttUp.df$logFC  >= 0.58 & ttUp.df$logPval >= 1.3),]
write.table(ttUp.df, file=paste("StatsUpregulated.csv"), quote=FALSE, sep=',', row.names = FALSE)
rm(ttUp.df)

ttDown.df <- topTable(fit2, number=Inf, sort.by ='p', p.value=1)[, c(1, 4, 5)]
ttDown.df$symbol <- unlist(mget(rownames(ttDown.df), uniprotmir5a62sym, ifnotfound=rownames(ttDown.df)))
ttDown.df$FC <- ifelse(ttDown.df$logFC >= 0, inv.glog2(ttDown.df$logFC), -inv.glog2(-ttDown.df$logFC))
ttDown.df$logPval <- -log10(ttDown.df[,c(2)])
ttDown.df <- ttDown.df[which(ttDown.df$logFC  <= -0.58 & ttDown.df$logPval >= 1.3),]
write.table(ttDown.df, file=paste("StatsDownregulated.csv"), quote=FALSE, sep=',', row.names = FALSE)
rm(ttDown.df)

tt.df <- topTable(fit2, number=Inf, sort.by ='p', p.value=1)[, c(1, 4, 5)]
tt.df$symbol <- unlist(mget(rownames(tt.df), uniprotmir5a62sym, ifnotfound=rownames(tt.df)))
tt.df$FC <- ifelse(tt.df$logFC >= 0, inv.glog2(tt.df$logFC), -inv.glog2(-tt.df$logFC))
tt.df$logPval <- -log10(tt.df[,c(2)])

write.table(tt.df, file=paste("StatsTable.csv"), quote=FALSE, sep=',', row.names = FALSE)

 
v <- ggplot(tt.df,aes(x=logFC,y=logPval)) + geom_point(size=2, alpha=1, col='black') +
  theme_update(plot.title=element_text(hjust=0.5), legend.position='none') +
  geom_point(data=tt.df, stat='identity', aes(colour=cut(logFC, c(-Inf,-1,1,100))), size=1) + geom_hline(yintercept=-log(0.05,10), linetype="3313", colour="grey") + geom_vline(xintercept=1, linetype="3313", colour="grey") + geom_vline(xintercept=-1, linetype="3313", colour="grey") +
  scale_color_manual(name = "logFC",
                     values = c("(-Inf,-1]" = "blue",
                                "(-1,1]" = "gray",
                                "(1,100]" = "red"),
                     labels = c("decreased", "insignificant", "increased")) +
  theme_classic()
tiff("VolcanoPlot.tiff", width = 6, height = 8, units = 'in', res=200)
plot(v)
dev.off()

