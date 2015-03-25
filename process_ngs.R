lines = c('MCF10A', 'MCF7', 'MDA231')
ac_dat = list()
me_dat = list()
for(l in lines){
  ac_dir = dir('data_raw/ngsplot_output', full.names = T, pattern = paste(l,'.+','AC',sep = ""))
  me_dir = dir('data_raw/ngsplot_output', full.names = T, pattern = paste(l,'.+','ME3',sep = ""))
  ac_prof = paste(ac_dir, '/hm1.txt', sep = "")
  me_prof = paste(me_dir, '/hm1.txt', sep = "")
  
  tmp = read.table(ac_prof, stringsAsFactors = F)
  ensgs = tmp[2:nrow(tmp),1]
  
  strand = tmp[2:nrow(tmp),4]
  ac = as.matrix(tmp[2:nrow(tmp), 5:ncol(tmp)])
  rownames(ac) = ensgs
  ac_dat[[length(ac_dat) + 1]] = ac
  
  tmp = read.table(me_prof, stringsAsFactors = F)
  me = as.matrix(tmp[2:nrow(tmp), 5:ncol(tmp)])
  rownames(me) = ensgs
  me_dat[[length(me_dat) + 1]] = me
} 
names(ac_dat) = lines
names(me_dat) = lines
save(ac_dat, me_dat, file = 'data_raw//ngs_k4_data.save')