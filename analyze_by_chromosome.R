source('read_countDir.R')
source('norm_counts.R')
source('jrb_matrix.R')
source('heatmap.3-split.R')
exData = countDir.example()
exData = norm.remove_unaligned(exData)
exData = norm.read_depth(exData)

chrSizes = read.table('hg38.chrom.sizes', row.names = 1)
exData = norm.feature_length(exData, chrSizes)


tmp = jrb.split_colnames(exData)
cell_lines = tmp[1,]
histone_mods = tmp[2,]
colnames(exData) = paste(cell_lines, histone_mods)

heatmap.2((t(exData)+1), trace = 'n', scale = 'c', Rowv = F, margins = c(5,12), Colv = T, dendrogram = 'c')

feData = norm.fe(exData)
feData = feData[rownames(feData)!='chrY',]

heatmap.2((t(feData)+1), trace = 'n', scale = 'n', Rowv = T, margins = c(5,12))

o = order(histone_mods)
forBoxplot = matrix(0, nrow = nrow(exData), ncol = 0)
for(l in unique(cell_lines)){
  keep = cell_lines == l
}

tmp = t(exData)
o = mixedorder(rownames(exData))
tmp = tmp[,o]
tmp = tmp[,colnames(tmp) != 'chrY']
forBoxplot = cbind(tmp[1:7,], tmp[8:14,], tmp[15:21,])
forBoxplot = forBoxplot[,mixedorder(colnames(forBoxplot))]
heatmap(forBoxplot, scale = 'c', Colv = NA)
boxplot(forBoxplot[,],  range = 0, border = rep(c('red', 'green', 'blue'),24))

library(ggplot2)

chrms = factor(sub('chr', '', rep(colnames(tmp), nrow(tmp))), levels = sub('chr', '', colnames(tmp)))
i = unlist(lapply(1:nrow(tmp), function(x)rep(x, ncol(tmp))))



pdf('reads by cell line.pdf')
is_diff = rep('other', length(i))
is_diff[chrms == '17'] = 'chr17'
is_diff[chrms == '20'] = 'chr20'
is_diff[chrms == '8'] = 'chr8'
is_diff[chrms == '19'] = 'chr19'

plot_data = data.frame(values = t(tmp)[1:length(tmp)], chrms = chrms, cell_lines = cell_lines[i], is_input = ifelse(histone_mods[i]=='input', 'input','histone mark'),  histone_mods = histone_mods[i], is_diff = is_diff)

p <- ggplot(data = plot_data, mapping = aes(x = cell_lines, y = values, color = is_diff, shape = is_input))
p = p + scale_colour_manual(values=c(RColorBrewer::brewer.pal(length(unique(is_diff))-1, 'Dark2'), '#ffffff'))
print(p + layer(geom = "point", position = position_jitter(width = .04, height = 0)))

p <- ggplot(data = plot_data, mapping = aes(x = cell_lines, y = values, color = is_diff, shape = is_input))
p = p + scale_colour_manual(values=c(RColorBrewer::brewer.pal(length(unique(is_diff))-1, 'Dark2'), '#808080'))
print(p + layer(geom = "point", position = position_jitter(width = .04, height = 0)))
dev.off()

pdf('reads by chromosome.pdf')
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

dev.off()