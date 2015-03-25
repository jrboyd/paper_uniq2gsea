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

plotNGS = function(ENSGcut_list, list_name, invert = F, ymax = 4){
  #plot ngs profile style plots
  #ENSGcut_list : is a character vector of cut ensg ids, cut means version number removed
  #list_name : name or description of input list, used in title
  #invert : if T, everything not in list will be plotted
  #ymax : the upper ylim for plots
  layout(matrix(c(1:4), ncol = 1, byrow = T))
  par(mai = c(0,1,0,.2))
  xs = 0:100
  xs = (20 * xs) - 1000
  keep = ENSGcut_list
  if(length(keep) < 2 && is.na(keep)){
    keep = NGS_SYMBOL_SET
  }
  if(invert){
    tmp = rep(T, length(NGS_SYMBOL_SET))
    names(tmp) = NGS_SYMBOL_SET
    tmp[keep] = F
    keep = NGS_SYMBOL_SET[tmp]
  }
  plot(c(0,1),c(0,1), type = 'n', axes = F, xlab = '', ylab = '')
  text(.5,.5,paste(list_name, '\n', length(keep),' genes', sep = ''), cex = 1.5)
  legend(x = 'left',legend = lines, fill = l2col[lines], bty = 'n')
  plot(c(0,1), type = 'n', xlim = c(-1000,1000), ylim = c(0,ymax), ylab = 'H3K4ac', lwd = 2, xaxt = 'n')
  for(l in lines){
    ac_d = ac_dat[[l]]
    keep = intersect(keep, rownames(ac_d))
    lines(xs, colMeans(ac_d[keep,]), col = l2col[l], lwd = 2)
  }
  plot(c(0,1), type = 'n', xlim = c(-1000,1000), ylim = c(0,ymax), ylab = 'H3K4me3', lwd = 2)
  for(l in lines){
    me_d = me_dat[[l]]
    keep = intersect(keep, rownames(me_d))
    lines(xs, colMeans(me_d[keep,]), col = l2col[l], lwd = 2)
  }
  plot(c(0,1),c(0,1), type = 'n', axes = F, xlab = '', ylab = '')
}

load('data_raw///ngs_k4_data.save')

NGS_SYMBOL_SET = rownames(ac_dat[[1]])
load('ref//gsea_dataset.save')
load('data_intermediate/gsea_passing.save')
load('data_intermediate/uniq_passing.save')

sym2cut = function(sym){
  keep = sapply(ensg2sym, function(x)return(any(x == sym)))
  cut = ensg2cut[names(ensg2sym)[keep]]
  return(cut)
}

outdir_name = 'output3'
if(!file.exists(outdir_name))
  dir.create(outdir_name)
load('data_intermediate///uniquely_membership.save')
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
    pdf(fname, width = 10)
    toPlot = membership[,list_name]
    toPlot = names(toPlot)[toPlot] #character vector of gene symbols in list
    plotNGS(sym2cut(toPlot), list_name, ymax = 5)
    #plotNGS_wBG(fg_toPlot, bg_toPlot, list_name, uniq_name, ymax = 5)
    dev.off()
  }
}
