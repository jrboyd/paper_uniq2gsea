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



if(savePlot) pdf('reads by cell line - all marks.pdf')

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
plot_data = matrix2data.frame(tmp)

p <- ggplot(data = plot_data, mapping = aes(x = cell_lines, y = values, color = is_diff, shape = is_input))
p = p + scale_colour_manual(values=c(RColorBrewer::brewer.pal(length(unique(plot_data$is_diff))-1, 'Set1'), '#ffffff'))
print(p + layer(geom = "point", position = position_jitter(width = .04, height = 0)))

p <- ggplot(data = plot_data, mapping = aes(x = cell_lines, y = values, color = is_diff, shape = is_input))
p = p + scale_colour_manual(values=c(RColorBrewer::brewer.pal(length(unique(plot_data$is_diff))-1, 'Set1'), '#808080'))
print(p + layer(geom = "point", position = position_jitter(width = .04, height = 0)))
if(savePlot) dev.off()

if(savePlot) pdf('reads by chromosome - all marks.pdf')
p <- ggplot(data = plot_data, mapping = aes(x = chrms, y = values, color = cell_lines, shape = is_input))
plot1 = (p + 
  scale_colour_manual(values=RColorBrewer::brewer.pal(3, 'Dark2')) +
  #position_jitter() +
  ggtitle('Distribution of reads among chromosomes for 6 marks and input') +
  xlab('chromosome') + 
  ylab('chromosome RPKM') + 
  layer(geom = "point", position = position_jitter(width = .1, height = 0)))
# +
#   smry
print(plot1)




p <- ggplot(data = plot_data, mapping = aes(x = chrms, y = values, color = cell_lines))
smry = stat_summary(
  fun.data="mean_cl_boot", conf.int=0.95,
  geom="crossbar", width=0.8
)

plot2 = (p + 
  scale_colour_manual(values=RColorBrewer::brewer.pal(3, 'Dark2')) +
  #position_jitter() +
  ggtitle('Distribution of reads among chromosomes for 6 marks and input') +
  xlab('chromosome') + 
  ylab('chromosome RPKM') + 
  geom_boxplot(position = position_identity(), outlier.size = 0))

print(plot2)

plot3 = (p + 
  scale_colour_manual(values=RColorBrewer::brewer.pal(3, 'Dark2')) +
  #position_jitter() +
  ggtitle('Distribution of reads among chromosomes for 6 marks and input') +
  xlab('chromosome') + 
  ylab('chromosome RPKM') + 
  layer(geom = "point")+
  smry)

print(plot3)

if(savePlot) dev.off()

if(savePlot) pdf('reads genome - k4 only.pdf')
tmp = tmp[(histone_mods == 'H3K4ME3' | histone_mods == 'H3K4AC' | histone_mods == 'input'),]
plot_data = matrix2data.frame(tmp)

p <- ggplot(data = plot_data, mapping = aes(x = cell_lines, y = values, color = is_diff, shape = histone_mods))
p = p + scale_colour_manual(values=c(RColorBrewer::brewer.pal(length(unique(plot_data$is_diff))-1, 'Set1'), '#ffffff'))
print(p + layer(geom = "point", position = position_jitter(width = .04, height = 0)))

p <- ggplot(data = plot_data, mapping = aes(x = cell_lines, y = values, color = is_diff, shape = histone_mods))
p = p + scale_colour_manual(values=c(RColorBrewer::brewer.pal(length(unique(plot_data$is_diff))-1, 'Set1'), '#808080'))
print(p + layer(geom = "point", position = position_jitter(width = .04, height = 0)))

p <- ggplot(data = plot_data, mapping = aes(x = chrms, y = values, color = cell_lines, shape = histone_mods))
p = (p + 
           scale_colour_manual(values=RColorBrewer::brewer.pal(3, 'Dark2')) +
           #position_jitter() +
           ggtitle('Distribution of reads among chromosomes for 6 marks and input') +
           xlab('chromosome') + 
           ylab('chromosome RPKM') + 
           layer(geom = "point", position = position_dodge(width = .2, height = .1)))

print(p)

p <- ggplot(data = plot_data, mapping = aes(x = chrms, y = values, color = cell_lines))
smry = stat_summary(
  fun.data="mean_cl_boot", conf.int=0.95,
  geom="crossbar", width=0.8
)

plot2 = (p + 
           scale_colour_manual(values=RColorBrewer::brewer.pal(3, 'Dark2')) +
           #position_jitter() +
           ggtitle('Distribution of reads among chromosomes for 6 marks and input') +
           xlab('chromosome') + 
           ylab('chromosome RPKM') + 
           geom_boxplot(position = position_identity(), outlier.size = 0))

print(plot2)

plot3 = (p + 
           scale_colour_manual(values=RColorBrewer::brewer.pal(3, 'Dark2')) +
           #position_jitter() +
           ggtitle('Distribution of reads among chromosomes for 6 marks and input') +
           xlab('chromosome') + 
           ylab('chromosome RPKM') + 
           layer(geom = "point")+
           smry)

print(plot3)
if(savePlot) dev.off()