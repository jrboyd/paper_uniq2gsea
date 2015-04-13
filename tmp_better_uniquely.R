library(limma)
sym2dat = function(sym){
  ensg = names(ensg2sym)[ensg2sym == sym]
  return(markData_4me3_4ac[ensg,])
}
minDiff = log2(3)
load('data_raw/mycounts_data.save')
isMarked = matrix(F, nrow = nrow(markData_4me3_4ac), ncol = 0)
for(i in 1:3){
  isMarked = cbind(isMarked, apply(markData_4me3_4ac[,1:3], 1, function(x){
    sorted_x = sort(x)
    isUniq =  x[i] > (sorted_x[2] + minDiff)
    return(isUniq)
  }))
}
for(i in 1:3){
  isMarked = cbind(isMarked, apply(markData_4me3_4ac[,4:6], 1, function(x){
    sorted_x = sort(x)
    isUniq =  x[i] > (sorted_x[2] + minDiff)
    return(isUniq)
  }))
}
colnames(isMarked) = colnames(markData_4me3_4ac)
layout(1:2)
vennDiagram(isMarked[,1:3], cex = .7)
vennDiagram(isMarked[,4:6], cex = .7)

layout(matrix(1:9, nrow = 3))
for(i in 1:3){
  vennDiagram(isMarked[,c(i,4)])
  vennDiagram(isMarked[,c(i,4+1)])
  vennDiagram(isMarked[,c(i,4+2)])
}

load('ref//dictionaries.save')
crossMarked = apply(isMarked[,c(1,4)], 1, all)
crossMarked = names(crossMarked)[crossMarked]
ensg2sym[crossMarked]

uniquely_K4_membership = matrix(F, nrow = nrow(isMarked), ncol = 9)
uniquely_K4_membership[,1] = isMarked[,1] & !isMarked[,4]
uniquely_K4_membership[,2] = isMarked[,1] & isMarked[,4]
uniquely_K4_membership[,3] = !isMarked[,1] & isMarked[,4]
uniquely_K4_membership[,4] = isMarked[,2] & !isMarked[,5]
uniquely_K4_membership[,5] = isMarked[,2] & isMarked[,5]
uniquely_K4_membership[,6] = !isMarked[,2] & isMarked[,5]
uniquely_K4_membership[,7] = isMarked[,3] & !isMarked[,6]
uniquely_K4_membership[,8] = isMarked[,3] & isMarked[,6]
uniquely_K4_membership[,9] = !isMarked[,3] & isMarked[,6]

rownames(uniquely_K4_membership) = ensg2sym[rownames(isMarked)]
colnames(uniquely_K4_membership) = c('10a ac', '10a both', '10a me3', '7 ac', '7 both', '7 me3', 'mda ac', 'mda both', 'mda me3')
keep = rowSums(uniquely_K4_membership) > 0
uniquely_K4_membership = uniquely_K4_membership[keep,]
save(uniquely_K4_membership, file = 'data_intermediate/fe_uniquely.save')
