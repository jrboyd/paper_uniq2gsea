#starting with uniquely sets, tests them for enrichment using binom.test in gsea lists.
#a heatmap of -log10 pvals is written to pdf (gsea lists are rows, uniq groups are columns)
#gsea lists passing a threshold are added to gsea_passing.save
#the corresponding best uniquely group is added to uniq_passing.save

##parameters
#must be true to update input for ngsplots
updatePassing = F
#name of heatmap pdf
pdfName = 'output_pooled/gsea_important_genes.pdf'
savePlot = T
#starting pvalue - will be reduced until multiple gsea lists pass
pthresh = 9
#bg size, all detected?  all passing FE?  all genes?
bg_size = 20000
#gene most be present in at least this many gsea lists
min_shared = 4

##dependencies
load('data_intermediate///uniquely_membership.save')
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
gsea_membership_sig = gsea_membership[,rownames(dat)]
keep = rowSums(gsea_membership_sig) > 0
gsea_membership_sig = gsea_membership_sig[keep,]
uniq_in_gsea_sig = uniquely_K4_membership[intersect(rownames(gsea_membership_sig), rownames(uniquely_K4_membership)),]
uniq_in_gsea_sig = uniq_in_gsea_sig[,c(4,1,7,5,2,8,6,3,9)]
if(savePlot) pdf(pdfName, width = 12, height = 8)
for(i in 1:ncol(uniq_in_gsea_sig)){
  uniq_members = rownames(uniq_in_gsea_sig)[uniq_in_gsea_sig[,i]]
  toPlot = t(ifelse(gsea_membership_sig[uniq_members,], 1, 0))
  r_order = order(rowSums(toPlot), decreasing = T)
  c_order = order(colSums(toPlot), decreasing = T)
  toPlot = toPlot[r_order, c_order]
  
  
  j = 1
  tot = 0
  end = ncol(toPlot)
  while(tot < ncol(toPlot) & j <= nrow(toPlot)){
    start = tot + 1
    o = order(toPlot[j,start:end], decreasing = T)
    n = sum(toPlot[j,start:end])
    #   if(n == 0)
    #     next
    
    o = o + start - 1
    tot = tot + n 
    print(range(o))
    
    
    toPlot[,start:end] = toPlot[,o, drop = F]
    
    j = j + 1
    if(j > nrow(toPlot))
      next
    o = order(toPlot[j,start:tot], decreasing = F)
    o = o + min(start, tot) - 1
    print(range(o))
    toPlot[,start:tot] = toPlot[,o, drop = F]
  }
  #toPlot = rbind(toPlot, colSums(toPlot))
  heatmap.2((toPlot), scale = 'n', margins = c(4,34), cexCol = 20/ncol(toPlot),Colv = F, Rowv = F, col = c('white', 'black'), trace = 'n', dendrogram = 'n')
}

layout(1)
plot(0:1, 0:1, type = 'n', xlim = c(0,1), ylim = c(-0,10))
cr = colorRamp(c('gray','black', 'orange','red'))
for(i in 1:ncol(uniq_in_gsea_sig)){
  uniq_members = rownames(uniq_in_gsea_sig)[uniq_in_gsea_sig[,i]]
  toPlot = t(ifelse(gsea_membership_sig[uniq_members,], 1, 0))
  r_order = order(rowSums(toPlot), decreasing = T)
  c_order = order(colSums(toPlot), decreasing = T)
  toPlot = toPlot[r_order, c_order]
  
  
  j = 1
  tot = 0
  end = ncol(toPlot)
  while(tot < ncol(toPlot) & j <= nrow(toPlot)){
    start = tot + 1
    o = order(toPlot[j,start:end], decreasing = T)
    n = sum(toPlot[j,start:end])
    #   if(n == 0)
    #     next
    
    o = o + start - 1
    tot = tot + n 
    print(range(o))
    
    
    toPlot[,start:end] = toPlot[,o, drop = F]
    
    j = j + 1
    if(j > nrow(toPlot))
      next
    o = order(toPlot[j,start:tot], decreasing = F)
    o = o + min(start, tot) - 1
    print(range(o))
    toPlot[,start:tot] = toPlot[,o, drop = F]
  }
  
  colors = rgb(cr(colSums(toPlot) / max(colSums(toPlot)))/255)
  text(.1*i,9.8, colnames(uniq_in_gsea_sig)[i],adj = c(1,1), cex = .4, col = rgb(0,109/255,44/255))
  for(k in 1:ncol(toPlot)){
    col = colors[k]
    txt = colnames(toPlot)[k]
    text(.1*i,9.5-(k-1)*.1, txt,adj = c(1,1), cex = .4, col = col)
#     if(colSums(toPlot)[k] > 2)
#       text(.3,10-(k-1)*.1, txt,adj = c(1,1), cex = ncol(toPlot)/150, col = col)
  }
  
}


layout(1)
plot(0:1, 0:1, type = 'n', xlim = c(0,1), ylim = c(-0,10))
cr = colorRamp(c('gray','black', 'orange','red'))
for(i in 1:ncol(uniq_in_gsea_sig)){
  uniq_members = rownames(uniq_in_gsea_sig)[uniq_in_gsea_sig[,i]]
  toPlot = t(ifelse(gsea_membership_sig[uniq_members,], 1, 0))
  r_order = order(rowSums(toPlot), decreasing = T)
  c_order = order(colSums(toPlot), decreasing = T)
  toPlot = toPlot[r_order, c_order]
  
  
  j = 1
  tot = 0
  end = ncol(toPlot)
  while(tot < ncol(toPlot) & j <= nrow(toPlot)){
    start = tot + 1
    o = order(toPlot[j,start:end], decreasing = T)
    n = sum(toPlot[j,start:end])
    #   if(n == 0)
    #     next
    
    o = o + start - 1
    tot = tot + n 
    print(range(o))
    
    
    toPlot[,start:end] = toPlot[,o, drop = F]
    
    j = j + 1
    if(j > nrow(toPlot))
      next
    o = order(toPlot[j,start:tot], decreasing = F)
    o = o + min(start, tot) - 1
    print(range(o))
    toPlot[,start:tot] = toPlot[,o, drop = F]
  }
  colors = rgb(cr(colSums(toPlot) / max(colSums(toPlot)))/255)
  text(.1*i,9.8, colnames(uniq_in_gsea_sig)[i],adj = c(1,1), cex = .4, col = rgb(0,109/255,44/255))
  for(k in 1:ncol(toPlot)){
    col = colors[k]
    txt = colnames(toPlot)[k]
    #text(.1*i,10-(k-1)*.1, txt,adj = c(1,1), cex = .4, col = col)
         if(colSums(toPlot)[k] >= min_shared)
           text(.1*i,9.5-(k-1)*.1, txt,adj = c(1,1), cex = .4, col = col)
  }
  
}

layout(1)
plot(0:1, 0:1, type = 'n', xlim = c(0,1), ylim = c(-0,10))
cr = colorRamp(c('gray','black', 'orange','red'))
output = character(length = ncol(uniq_in_gsea_sig))
names(output) = colnames(uniq_in_gsea_sig)
for(i in 1:ncol(uniq_in_gsea_sig)){
  uniq_members = rownames(uniq_in_gsea_sig)[uniq_in_gsea_sig[,i]]
  toPlot = t(ifelse(gsea_membership_sig[uniq_members,], 1, 0))
  r_order = order(rowSums(toPlot), decreasing = T)
  c_order = order(colSums(toPlot), decreasing = T)
  toPlot = toPlot[r_order, c_order]
  
  
  j = 1
  tot = 0
  end = ncol(toPlot)
  while(tot < ncol(toPlot) & j <= nrow(toPlot)){
    start = tot + 1
    o = order(toPlot[j,start:end], decreasing = T)
    n = sum(toPlot[j,start:end])
    #   if(n == 0)
    #     next
    
    o = o + start - 1
    tot = tot + n 
    print(range(o))
    
    
    toPlot[,start:end] = toPlot[,o, drop = F]
    
    j = j + 1
    if(j > nrow(toPlot))
      next
    o = order(toPlot[j,start:tot], decreasing = F)
    o = o + min(start, tot) - 1
    print(range(o))
    toPlot[,start:tot] = toPlot[,o, drop = F]
  }
  colors = rgb(cr(colSums(toPlot) / max(colSums(toPlot)))/255)
  text(.1*i,9.8, colnames(uniq_in_gsea_sig)[i],adj = c(1,1), cex = .4, col = rgb(0,109/255,44/255))
  l = 1
  for(k in 1:ncol(toPlot)){
    col = colors[k]
    txt = colnames(toPlot)[k]
    #text(.1*i,10-(k-1)*.1, txt,adj = c(1,1), cex = .4, col = col)
    if(colSums(toPlot)[k] >= min_shared){
      text(.1*i,9.5-(l-1)*.25, txt,adj = c(1,1), cex = .7, col = col)
      l = l + 1
    }
  }
  keep = colSums(toPlot) >= min_shared
  tmp = colnames(toPlot)[keep]
  output[i] = paste(tmp, collapse = ',')
  
}
output = paste(names(output), output, sep = ',', collapse = '\n')
write.table(output, file = 'important_genes.csv', col.names = F, row.names = F, quote = F)
if(savePlot) dev.off()



