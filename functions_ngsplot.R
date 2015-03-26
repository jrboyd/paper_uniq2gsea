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
  bg_lty = 2
  bg_pch = 21
  
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
    #points(xs, colMeans(ac_d[fg_keep,]), col = l2col[l], pch = 19, cex = 1)
    #points(xs, colMeans(ac_d[fg_keep,]), col = l2col_bg[l], pch = 16, cex = 1)
    lines(xs, colMeans(ac_d[fg_keep,]), col = l2col[l], lwd = bg_lwd, lty = bg_lty)
  }
  for(l in lines){
    ac_d = ac_dat[[l]]
    lines(xs, colMeans(ac_d[bg_keep,]), col = l2col[l], lwd = fg_lwd, lty = 1)
  }
  
  plot(c(0,1), type = 'n', xlim = c(-1000,1000), ylim = c(0,ymax), ylab = 'H3K4me3', lwd = 2)
  for(l in lines){
    me_d = me_dat[[l]]
    #points(xs, colMeans(me_d[fg_keep,]), col = l2col[l], pch = 19, cex = 1)
    #points(xs, colMeans(me_d[fg_keep,]), col = l2col_bg[l], pch = 16, cex = 1)
    lines(xs, colMeans(me_d[fg_keep,]), col = l2col[l], lwd = bg_lwd, lty = bg_lty)
  }
  for(l in lines){
    me_d = me_dat[[l]]
    lines(xs, colMeans(me_d[bg_keep,]), col = l2col[l], lwd = fg_lwd, lty = 1)
  }
  
  plot(c(0,1),c(0,1), type = 'n', axes = F, xlab = '', ylab = '')
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

sym2cut = function(sym){
  keep = sapply(ensg2sym, function(x)return(any(x == sym)))
  cut = ensg2cut[names(ensg2sym)[keep]]
  return(cut)
}