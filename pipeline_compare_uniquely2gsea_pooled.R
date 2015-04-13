#starting with uniquely sets, tests them for enrichment using binom.test in gsea lists.
#a heatmap of -log10 pvals is written to pdf (gsea lists are rows, uniq groups are columns)
#gsea lists passing a threshold are added to gsea_passing.save
#the corresponding best uniquely group is added to uniq_passing.save

##parameters
#must be true to update input for ngsplots
updatePassing = F
#name of heatmap pdf
pdfName = 'output_pooled/gsea_uniquely_enriched.pdf'
savePlot = T
#starting pvalue - will be reduced until multiple gsea lists pass
pthresh = 9
#bg size, all detected?  all passing FE?  all genes?
bg_size = 20000

##dependencies
#load('data_intermediate///uniquely_membership.save')
load('data_intermediate/fe_uniquely.save')
load('ref//gsea_dataset.save')
source('scripts/heatmap.3-split.R')

if(!file.exists(dirname(pdfName)))
  dir.create(dirname(pdfName))

gsea_passing = list()
uniq_passing = list()

list_num = 1

all_genes = character()
ncols = 0
for(gm in gsea_membs){
  all_genes = union(all_genes, rownames(gm))
  ncols = ncols + ncol(gm)
}
gsea_membership = matrix(F, ncol = ncols, nrow = length(all_genes))
rownames(gsea_membership) = all_genes
colnames(gsea_membership) = 1:ncol(gsea_membership)

ncolsnamed = 0
for(gm in gsea_membs){
  start_col = ncolsnamed + 1
  end_col = ncolsnamed + ncol(gm)
  ncolsnamed = ncolsnamed + ncol(gm)
  colnames(gsea_membership)[start_col:end_col] = colnames(gm)
  gsea_membership[rownames(gm),colnames(gm)] = gm
}

lsize = ncol(gsea_membership)
list_num = list_num + 1

gsea_sizes = colSums(gsea_membership)
layout(matrix(1:9, nrow = 3))
toPlot = matrix(0, nrow = ncol(uniquely_K4_membership), ncol = ncol(gsea_membership))
rownames(toPlot) = colnames(uniquely_K4_membership)
colnames(toPlot) = colnames(gsea_membership)

for(j in 1:ncol(uniquely_K4_membership)){
  keep = rownames(uniquely_K4_membership)[uniquely_K4_membership[,j]]#gene symbols positive in this uniquely group
  uniq_size = length(keep)
  keep = intersect(keep, rownames(gsea_membership))#filter out symbols that don't occur in gsea
  gsea_hits = colSums(gsea_membership[keep,])#number of genes in from uniquely group that occur in each gsea list
  perc_in_gsea = gsea_hits / length(keep)
  #o = order(perc_in_gsea, decreasing = T)
  #plot(perc_in_gsea[o])
  #title(paste(colnames(uniquely_K4_membership)[i], colnames(gsea_membership)[o][1], sep = '\n'))
  toPlot[j,] = gsea_hits
}
#save(uniquely_K4_membership, file = "uniquely_H3K4.save") 
perc_of_bg = colSums(uniquely_K4_membership) / bg_size

pvals = toPlot
for(i in 1:nrow(toPlot)){
  for(j in 1:ncol(toPlot)){
    tmp = binom.test(toPlot[i,j], n = sum(gsea_membership[,j]), p = perc_of_bg[i], alternative = 'greater')
    pvals[i,j] = tmp$p.value
  }
}
dat = -log10(t(pvals))

keep = F
pthresh = pthresh + 1#compensate for first while loop
while(sum(keep) < 5){#find most stringent pval with results
  pthresh = pthresh - 1
  keep = apply(dat, 1, max) > pthresh
  
}

tmp = apply(dat[keep,], 1, function(x)return((1:length(x))[x == max(x)]))#uniq category with best pval in kept
tmp = colnames(dat)[tmp]
uniq_passing = tmp
o = c(4,1,7,5,2,8,6,3,9)
o = 1:9
cr = colorRamp(c('white','blue'))
colors = cr(x = (0:13^2)^.5/10 - .3)
colors = ifelse(is.na(colors), 255, colors)
colors = rgb(colors / 255)
byChance = round(lsize * 10^-pthresh, 2) * 9
title = paste(
  'Of all gsea lists',
  '\n', sum(keep), ' of ', lsize, ' lists pass 10^-', pthresh, 
  '\n', format(byChance, nsmall = 2), ' expected by chance',
  sep = '')


dat = dat[keep,o]
clust = hclust(dist(dat))
dat = dat[rev(clust$order),]
noteText = t(toPlot[o,rownames(dat)])

cr = colorRamp(c('black', 'red'))
noteColors = cr(t(dat) / max(dat))
noteColors = rgb(noteColors / 255)
noteColors = t(matrix(noteColors, nrow = length(o), ncol = sum(keep), byrow = F))
tmp = character()
for(i in nrow(dat):1){
  tmp = c(tmp, noteColors[i,])
}
noteColors = tmp
#noteColors = rgb(cr(1:10/10) / 255)
if(savePlot) pdf(pdfName, width = 9, height = nrow(dat)/6+3)
heatmap.2(
  dat,
  labRow = paste0('(', colSums(gsea_membership[,rownames(dat)]), ') ', rownames(dat)),
  labCol = paste0(colnames(dat), ' (', colSums(uniquely_K4_membership[,colnames(dat)]), ') '),
  cellnote = noteText,
  notecol = noteColors,
  scale = 'n',
  dendrogram = 'n',
  Colv = NA,
  Rowv = NA,
  margins = c(8,20),
  cexRow = 1, cexCol = 1.4,
  trace = 'n',
  col = colors,
  colsep = c(3,6),
  sepcolor = 'black',
  main = title,
  key.title = '', density.info = 'n')
gsea_passing = names(keep)[keep]
if(savePlot) dev.off()

if(updatePassing){
  save(gsea_passing, file = 'data_intermediate/gsea_passing_pooled.save')
  save(uniq_passing, file = 'data_intermediate/uniq_passing_pooled.save')
}
