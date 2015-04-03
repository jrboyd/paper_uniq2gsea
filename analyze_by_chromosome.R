library(ggplot2)
source('jrb_R_scripts/read_countDir.R')
source('jrb_R_scripts/norm_counts.R')
source('jrb_R_scripts/jrb_matrix.R')
source('jrb_R_scripts/heatmap.3-split.R')

savePlot = T

exData = countDir.example()
exData = norm.remove_unaligned(exData)
exData = norm.read_depth(exData)

chrSizes = read.table('jrb_R_scripts/hg38.chrom.sizes', row.names = 1)
exData = norm.feature_length(exData, chrSizes)


tmp = jrb.split_colnames(exData)
cell_lines = tmp[1,]
histone_mods = tmp[2,]
colnames(exData) = paste(cell_lines, histone_mods)





matrix2data.frame = function(mat){
  chrms = factor(sub('chr', '', rep(colnames(mat), nrow(mat))), levels = sub('chr', '', colnames(mat)))
  
  i = unlist(lapply(1:nrow(mat), function(x)rep(x, ncol(mat))))#appends each row into single vector
  
  is_diff = rep('other', length(i))
  is_diff[chrms == '17'] = 'chr17'
  is_diff[chrms == '20'] = 'chr20'
  is_diff[chrms == '8'] = 'chr8'
  is_diff[chrms == '19'] = 'chr19'
  
  cl = jrb.split_colnames(t(mat), split = ' ')[1,]
  cl = cl[i]
  
  hm = jrb.split_colnames(t(mat), split = ' ')[2,]
  hm = hm[i]
  
  plot_data = data.frame(values = t(mat)[1:length(mat)], chrms = chrms, cell_lines = cl, is_input = ifelse(histone_mods[i]=='input', 'input','histone mark'),  histone_mods = hm, is_diff = is_diff)
  return(plot_data)
}

tmp = t(exData)
o = mixedorder(rownames(exData))
tmp = tmp[,o]
tmp = tmp[,colnames(tmp) != 'chrY']
k4me3_dat = tmp[histone_mods == 'H3K4ME3',]
k4ac_dat = tmp[histone_mods == 'H3K4AC',]
tmp = tmp[histone_mods == 'H3K4AC' | histone_mods == 'H3K4ME3',]
plot_data = matrix2data.frame(tmp)



p <- ggplot(data = plot_data, mapping = aes(x = cell_lines, y = values, color = is_diff))
p = p + scale_colour_manual(values=c(RColorBrewer::brewer.pal(length(unique(plot_data$is_diff))-1, 'Set1'), '#595959'))

if(savePlot) pdf('k4 reads per chromosome.pdf', height = 3, width = 6, useDingbats = F)
print(p + geom_point(position = position_jitter(width = .04, height = 0)) + facet_grid(. ~ histone_mods, scales = "free"))
if(savePlot) dev.off()