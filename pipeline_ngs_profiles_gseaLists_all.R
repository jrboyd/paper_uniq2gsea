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

primeColors = list(#from color brewer
  c(27,158,119),
  c(217,95,2),
  c(117,112,176),
  c(102,166,30),
  c(230,171,2),
  c(231,41,138)
)
names(primeColors) = colnames(markData_4me3_4ac)

secondColors = list(
  c(102,194,165),
  c(252,141,98),
  c(141,160,203),
  c(166,216,84),
  c(255,217,47),
  c(231,138,195)
)
names(secondColors) = colnames(markData_4me3_4ac)
secondColors = primeColors #mask secondary colors

applyWindow = function(dat, win = 10){
  if(win < 2)
    return(dat)
  out = matrix(0, nrow = nrow(dat), ncol = ncol(dat)/win)
  for(i in 1:ncol(out)){
    start = (i-1)*win+1
    end = i*win
    #out[,i] = apply(dat[,start:end],1,median)
    out[,i] = rowMeans(dat[,start:end])
  }
  return(out)
}

plotNGS_wBG = function(fg_ENSGcut_list, bg_ENSGcut_list, list_name, sel_name = 'selected', invert = F, ymax = 4){
  #plot ngs profile style plots
  #ENSGcut_list : is a character vector of cut ensg ids, cut means version number removed
  #list_name : name or description of input list, used in title
  #invert : if T, everything not in list will be plotted
  #ymax : the upper ylim for plots
  layout(matrix(c(1:4), ncol = 1, byrow = T))
  par(mai = c(0,1,0,.2))
  fg_lwd = 3
  bg_lwd = 3
  
  xs = 0:100
  xs = (20 * xs) - 1000
  bg_keep = bg_ENSGcut_list
  fg_keep = fg_ENSGcut_list
  if(length(bg_keep) < 2 && is.na(bg_keep)){
    bg_keep = NGS_SYMBOL_SET
  }
  if(invert){
    tmp = rep(T, length(NGS_SYMBOL_SET))
    names(tmp) = NGS_SYMBOL_SET
    tmp[bg_keep] = F
    bg_keep = NGS_SYMBOL_SET[tmp]
  }
  if(length(fg_keep) < 2 && is.na(fg_keep)){
    fg_keep = NGS_SYMBOL_SET
  }
  if(invert){
    tmp = rep(T, length(NGS_SYMBOL_SET))
    names(tmp) = NGS_SYMBOL_SET
    tmp[fg_keep] = F
    fg_keep = NGS_SYMBOL_SET[tmp]
  }
  bg_keep = intersect(bg_keep, NGS_SYMBOL_SET)
  fg_keep = intersect(fg_keep, NGS_SYMBOL_SET)
  
  plot(c(0,1),c(0,1), type = 'n', axes = F, xlab = '', ylab = '')
  text(.7,.5,paste(list_name, '\n', length(bg_keep) + length(fg_keep),' genes', sep = ''), cex = 1.5)
  legend(x = 0, y = .7,legend = lines, fill = l2col[lines], bty = 'n')
  legend(x = 0, y = .4,legend = c(paste(length(fg_keep), 'genes in', sel_name), paste(length(bg_keep), 'other genes in list')), lty = c(3,1), bty = 'n', lwd = 2)
  plot(c(0,1), type = 'n', xlim = c(-1000,1000), ylim = c(0,ymax), ylab = 'H3K4ac', lwd = 2, xaxt = 'n')
  for(l in lines){
    ac_d = ac_dat[[l]]
    lines(xs, colMeans(ac_d[fg_keep,]), col = l2col_bg[l], lwd = fg_lwd, lty = 3)
  }
  for(l in lines){
    ac_d = ac_dat[[l]]
    lines(xs, colMeans(ac_d[bg_keep,]), col = l2col[l], lwd = bg_lwd, lty = 1)
  }
  
  plot(c(0,1), type = 'n', xlim = c(-1000,1000), ylim = c(0,ymax), ylab = 'H3K4me3', lwd = 2)
  for(l in lines){
    me_d = me_dat[[l]]
    lines(xs, colMeans(me_d[fg_keep,]), col = l2col_bg[l], lwd = fg_lwd, lty = 3)
  }
  for(l in lines){
    me_d = me_dat[[l]]
    lines(xs, colMeans(me_d[bg_keep,]), col = l2col[l], lwd = bg_lwd, lty = 1)
  }
  
  plot(c(0,1),c(0,1), type = 'n', axes = F, xlab = '', ylab = '')
}


load('data_raw//ngs_k4_data.save')
NGS_SYMBOL_SET = rownames(ac_dat[[1]])
load('ref//gsea_dataset.save')
load('data_intermediate//gsea_passing.save')
load('data_intermediate//uniq_passing.save')

sym2cut = function(sym){
  keep = sapply(ensg2sym, function(x)return(any(x == sym)))
  cut = ensg2cut[names(ensg2sym)[keep]]
  return(cut)
}

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
