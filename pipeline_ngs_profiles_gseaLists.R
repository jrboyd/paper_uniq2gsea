#outputs all ngs profiles

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
load('data_raw///ngs_k4_data.save')
NGS_SYMBOL_SET = rownames(ac_dat[[1]])
load('ref//gsea_dataset.save')
#load('data_intermediate/gsea_passing.save')
#load('data_intermediate/uniq_passing.save')

outdir_name = 'output_allngs'
if(!file.exists(outdir_name))
  dir.create(outdir_name)
#load('data_intermediate///uniquely_membership.save')
grps_toPlot = c('C6', 'C2')
for(i in 1:length(grps_toPlot)){#plot profile of genes in all passing lists, genes in uniquely groups masked
  grp_name = grps_toPlot[i]
  grp = gsea_passing[[grp_name]]
  dir_name = paste(outdir_name, '//ngs_', grp_name, sep = '')
  membership = gsea_membs[[grp_name]]
  
  if(!file.exists(dir_name))
    dir.create(dir_name)
  for(j in 1:length(grp)){
    list_name = grp[j]
    fname = paste(dir_name, '//', list_name, '.pdf', sep = '')
    if(file.exists(fname))
      next
    pdf(fname, width = 10)
    toPlot = membership[,list_name]
    toPlot = names(toPlot)[toPlot] #character vector of gene symbols in list
    plotNGS(sym2cut(toPlot), list_name, ymax = 5)
    #plotNGS_wBG(fg_toPlot, bg_toPlot, list_name, uniq_name, ymax = 5)
    dev.off()
  }
}
