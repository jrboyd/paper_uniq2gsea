#plot membership table of gsea lists by cell line cluster
#members must be in corresponding uniquely list
#plots full table and table requiring membership shared by > 50%

##dependencies
load('ref//gsea_dataset.save')
load('data_intermediate///uniquely_membership.save')
load('data_intermediate/gsea_passing_pooled.save')
load('data_intermediate/uniq_passing_pooled.save')
source('scripts//heatmap.3-split.R')

o = c(4,1,7,5,2,8,6,3,9)#standard line order, ac - both - me
uniquely_K4_membership = uniquely_K4_membership[,o]

clust_mcf7 = c(
  'CHARAFE_BREAST_CANCER_LUMINAL_VS_MESENCHYMAL_UP',
  'SMID_BREAST_CANCER_BASAL_DN',
  'CHARAFE_BREAST_CANCER_LUMINAL_VS_BASAL_UP',
  'GOZGIT_ESR1_TARGETS_DN',
  'CREIGHTON_ENDOCRINE_THERAPY_RESISTANCE_1'
)

clust_10a = c(
  'CHARAFE_BREAST_CANCER_LUMINAL_VS_BASAL_DN',
  'COLDREN_GEFITINIB_RESISTANCE_DN',
  'ONDER_CDH1_TARGETS_2_DN'
)

clust_mda = c(
  'ONDER_CDH1_TARGETS_2_UP',
  'LIM_MAMMARY_STEM_CELL_UP',
  'CHICAS_RB1_TARGETS_CONFLUENT',
  'LEF1_UP.V1_UP',
  'CHARAFE_BREAST_CANCER_LUMINAL_VS_MESENCHYMAL_DN'
)

clusters = list(clust_10a, clust_mcf7, clust_mda)
names(clusters) = c('MCF10A', 'MCF7', 'MDA231')

tmp = append(gsea_lists[[1]], gsea_lists[[2]])
for(i in 3:length(gsea_lists)){
  tmp = append(tmp, gsea_lists[[i]])
}
gsea_lists = tmp



all_intersect = character()

for(i in 1:length(gsea_passing)){
  gsea_name = gsea_passing[i]
  uniq_name = uniq_passing[i]
  gsea_list = gsea_lists[[gsea_name]]
  uniq_list = uniquely_K4_membership[,uniq_name]
  uniq_list = names(uniq_list)[uniq_list]
  keep = intersect(gsea_list, uniq_list)
  all_intersect = union(all_intersect, keep)
}

gsea_toPlot = matrix(F, nrow = length(all_intersect), ncol = length(gsea_passing))
rownames(gsea_toPlot) = all_intersect
colnames(gsea_toPlot) = gsea_passing
for(g in gsea_passing){
  keep = intersect(rownames(gsea_toPlot), gsea_lists[[g]])
  gsea_toPlot[keep,g] = T
}

uniq_toPlot = matrix(F, nrow = length(all_intersect), ncol = ncol(uniquely_K4_membership))
rownames(uniq_toPlot) = all_intersect
colnames(uniq_toPlot) = colnames(uniquely_K4_membership)
for(un in colnames(uniquely_K4_membership)){
  uniqs = uniquely_K4_membership[,un]
  uniqs = names(uniqs)[uniqs]
  keep = intersect(rownames(uniq_toPlot), uniqs)
  uniq_toPlot[keep,un] = T
}

plotMembership = function(memb, title = ''){
  tmp = ifelse(memb, 1, 0)
  tmp = tmp[,order(colSums(tmp), decreasing = T)]
  for(i in ncol(tmp):1){
    tmp = tmp[order(tmp[,i], decreasing = T),]
  }
  tmp = tmp[order(rowSums(tmp), decreasing = T),]
  heatmap.2(tmp, scale = 'n', trace = 'n', dendrogram = 'n', Colv = F, Rowv = F, col = c('white', 'blue'), cexRow = 40/nrow(tmp), cexCol = .6, main = title, margins = c(20,16), key = F)  
}

# plotMembership(gsea_toPlot)
# plotMembership(uniq_toPlot)

gsea2grp = matrix(unlist(strsplit(uniq_passing, split = '_')), nrow = 2)[1,]
names(gsea2grp) = gsea_passing

lines = unique(gsea2grp)
pdf('heatmap_finalshared.pdf', height = 12, width = 10)
for(l in lines){
  keep = gsea2grp == l
  tmp = gsea_toPlot[,keep]
  keep = rowSums(tmp) > 0
  plotMembership(tmp[keep,], l)
  keep = rowSums(tmp) >= (ncol(tmp) / 2)
  plotMembership(tmp[keep,], l)
}
dev.off()
