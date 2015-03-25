#reads c2 and c6 gmt files from gsea
#save assembled gsea membership tables, gene lists, and union list in ref/gsea_dataset.save
#c2 is split into into subgroups KEGG, BIOCARTA, REACTOME, PID, rest of C2
#each output is divided into subgroup (6 total)

c2 = read.table('ref//c2.all.v5.0.symbols.gmt', sep = '\n', stringsAsFactors = F)[,1]
gsea_lists = list()
gsea_unions = list()
keys = c('KEGG', 'REACTOME', 'BIOCARTA', 'PID', 'C6', 'C2')
for(k in keys){
  gsea_lists[[length(gsea_lists) + 1]] = list()
  names(gsea_lists)[length(gsea_lists)] = k
  gsea_unions[[length(gsea_unions) + 1]] = character()
  names(gsea_unions)[length(gsea_unions)] = k
}
for(l in c2){
  tmp = strsplit(l, '\t')[[1]]
  nam = tmp[1]
  memb = tmp[3:length(tmp)]
  key = strsplit(nam, '_')[[1]][1]
  if(!any(key == keys))#is an uncategorized C2 list
    key = 'C2'
  gsea_unions[[key]] = union(gsea_unions[[key]], memb)
  gsea_lists[[key]][[length(gsea_lists[[key]]) + 1]] = memb
  names(gsea_lists[[key]])[length(gsea_lists[[key]])] = nam
}

c6 = read.table('ref//c6.all.v5.0.symbols.gmt', sep = '\n', stringsAsFactors = F)[,1]
for(l in c6){
  key = 'C6'
  tmp = strsplit(l, '\t')[[1]]
  nam = tmp[1]
  memb = tmp[3:length(tmp)]
  gsea_unions[[key]] = union(gsea_unions[[key]], memb)
  gsea_lists[[key]][[length(gsea_lists[[key]]) + 1]] = memb
  names(gsea_lists[[key]])[length(gsea_lists[[key]])] = nam
}

gsea_membs = list()
for(i in 1:length(gsea_lists)){
  rnames = gsea_unions[[i]]
  grp_lists = gsea_lists[[i]]
  members = matrix(F, nrow = length(rnames), ncol = length(grp_lists))
  rownames(members) = rnames
  colnames(members) = names(grp_lists)
  for(j in 1:length(grp_lists)){
    key = grp_lists[[j]]
    nam = names(grp_lists)[j]
    members[key, nam] = T
  }
  gsea_membs[[i]] = members
}
names(gsea_membs) = names(gsea_lists)

save(gsea_lists, gsea_unions, gsea_membs, file = 'ref//gsea_dataset.save')
