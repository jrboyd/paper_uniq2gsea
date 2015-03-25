#outputs all ngs profiles for gsea lists in gsea_passing.save.
#separately plots groups in uniq_passing.save, highlighting whether mark trend hold for full gsea list.

lines = c('MCF10A', 'MCF7', 'MDA231')
library(RColorBrewer)
l2col = RColorBrewer::brewer.pal(n = 3, 'Dark2')
names(l2col) = lines
l2col_bg = RColorBrewer::brewer.pal(n = 3, 'Set2')
names(l2col_bg) = lines
mods = c('H3K4AC', 'H3K4ME3')
other_mod = function(m){
  if(m == mods[1])
    return(mods[2])
  else if(m == mods[2])
    return(mods[1])
  return('unrecognized mod')
}

load('data_raw//mycounts_data.save')
load('ref//dictionaries.save')
load('data_raw//ngs_k4_data.save')
NGS_SYMBOL_SET = rownames(ac_dat[[1]])
load('ref//gsea_dataset.save')
load('data_intermediate//gsea_passing.save')
load('data_intermediate//uniq_passing.save')

outdir_name = 'output_all'
if(!file.exists(outdir_name))
  dir.create(outdir_name)
load('data_intermediate//uniquely_membership.save')
grp_name = 'All'
grp = gsea_passing
dir_name = paste(outdir_name, '//ngs_', grp_name, sep = '')
#membership = gsea_membs[[grp_name]]
tmp = append(gsea_lists[[1]], gsea_lists[[2]])
for(i in 3:length(gsea_lists)){
  tmp = append(tmp, gsea_lists[[i]])
}
all_gsea_lists = tmp

grp_genes = character()
for(i in 1:length(grp)){
  grp_genes = union(grp_genes, all_gsea_lists[[grp[i]]])
}

membership = matrix(F, nrow = length(grp_genes), ncol = length(grp))
rownames(membership) = grp_genes
colnames(membership) = grp
for(i in 1:length(grp)){
  membership[all_gsea_lists[[grp[i]]],grp[i]] = T
}

best_uniques = uniq_passing

if(!file.exists(dir_name))
  dir.create(dir_name)
for(j in 1:length(grp)){
  list_name = grp[j]
  fname = paste(dir_name, '//', list_name, '.pdf', sep = '')
  pdf(fname, width = 10)
  toPlot = membership[,list_name]
  toPlot = names(toPlot)[toPlot] #character vector of gene symbols in list
  best_uniq = best_uniques[j]
  isSigUniq = uniquely_K4_membership[,best_uniq]
  uniq_name = sub('_uniqAConly', ' uniquely marked by ac, not me3', best_uniq)
  uniq_name = sub('_uniqMEonly', ' uniquely marked by me3, not ac', uniq_name)
  uniq_name = sub('_both', ' uniquely marked by both me3 and ac', uniq_name)
  if(sum(isSigUniq) == 0){
    next
  }
  isUniq = uniquely_K4_membership[intersect(toPlot, rownames(uniquely_K4_membership[isSigUniq,])),]
  fg_toPlot = rownames(isUniq) #those genes that are uniq in a cell line
  isBg = rep(T, length(toPlot))
  names(isBg) = toPlot
  isBg[fg_toPlot] = F
  bg_toPlot = toPlot[isBg]#those genes that are not uniq in a cell line
  
  fg_toPlot = sym2cut(fg_toPlot)
  bg_toPlot = sym2cut(bg_toPlot)
  plotNGS_wBG(fg_toPlot, bg_toPlot, list_name, uniq_name, ymax = 5)
  dev.off()
}
