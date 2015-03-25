#starting with uniquely sets, tests them for enrichment using binom.test in gsea lists.
#a heatmap of -log10 pvals is written to pdf (gsea lists are rows, uniq groups are columns)
#gsea lists passing a threshold are added to gsea_passing.save
#the corresponding best uniquely group is added to uniq_passing.save

##parameters
#name of heatmap pdf
pdfName = 'each_gsea_uniquely_enriched.pdf'
#starting pvalue - will be reduced until multiple gsea lists pass
pthresh = 9
#bg size, all detected?  all passing FE?  all genes?
bg_size = 20000

##dependencies
load('data_intermediate///uniquely_membership.save')
load('ref//gsea_dataset.save')
source('scripts/heatmap.3-split.R')

gsea_passing = list()
uniq_passing = list()
pdf(pdfName, width = 9)
list_num = 1
for(gsea_membership in gsea_membs){
  lname = names(gsea_membs)[list_num]
  lsize = ncol(gsea_membership)
  list_num = list_num + 1
  
  gsea_sizes = colSums(gsea_membership)
  layout(matrix(1:9, nrow = 3))
  toPlot = matrix(0, nrow = ncol(uniquely_K4_membership), ncol = ncol(gsea_membership))
  rownames(toPlot) = colnames(uniquely_K4_membership)
  colnames(toPlot) = colnames(gsea_membership)
  
  for(j in 1:ncol(uniquely_K4_membership)){
    keep = rownames(uniquely_K4_membership)[uniquely_K4_membership[,j]]
    uniq_size = length(keep)
    keep = intersect(keep, rownames(gsea_membership))
    gsea_hits = colSums(gsea_membership[keep,])
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
  pthresh = pthresh + 1#compensate for first while loop
  keep = F
  while(sum(keep) < 5){#find most stringent pval with results
    pthresh = pthresh - 1
    keep = apply(dat, 1, max) > pthresh
  }
  
  tmp = apply(dat[keep,], 1, function(x)return((1:length(x))[x == max(x)]))#uniq category with best pval in kept
  tmp = colnames(dat)[tmp]
  uniq_passing[[length(uniq_passing) + 1]] = tmp
  o = c(4,1,7,5,2,8,6,3,9)
  cr = colorRamp(c('white','blue'))
  colors = cr(x = (0:13^2)^.5/10 - .3)
  colors = ifelse(is.na(colors), 255, colors)
  colors = rgb(colors / 255)
  byChance = round(lsize * 10^-pthresh, 2) * 9
  title = paste(
    lname,
    '\n', sum(keep), ' of ', lsize, ' lists pass 10^-', pthresh, 
    '\n', format(byChance, nsmall = 2), ' expected by chance',
    sep = '')
  heatmap.2(dat[keep,o], scale = 'n', dendrogram = 'n',  Colv = NA, margins = c(8,20), cexRow = .7, cexCol = .6, trace = 'n', col = colors, main = title, key.title = '', density.info = 'n')
  gsea_passing[[length(gsea_passing) + 1]] = names(keep)[keep]
}
dev.off()
names(gsea_passing) = names(gsea_membs)
save(gsea_passing, file = 'data_intermediate/gsea_passing.save')
names(uniq_passing) = names(gsea_membs)
save(uniq_passing, file = 'data_intermediate/uniq_passing.save')
