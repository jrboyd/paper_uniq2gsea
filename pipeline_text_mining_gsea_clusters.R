#read descriptions for gsea lists by hardcoded clusters
#assemble into corpus and process text
#construct comparison wordcloud

##dependecies
library(wordcloud)
library(tm)

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

all_desc = read.table('ref/all_desc.txt', sep = '\t', quote = '', row.names = 1, stringsAsFactors = F)
all_corpus = character(length = 3)
for(i in 1:length(clusters)){
  nam = names(clusters)[i]
  clust = clusters[[i]]
  corpus = paste(all_desc[clust, 1:2], collapse = ' ')
  corpus = paste(corpus, collapse = ' ')
  
  all_corpus[i] = corpus
}



corp = VCorpus(VectorSource(all_corpus))
corp <- tm_map(corp, removePunctuation)
corp <- tm_map(corp, content_transformer(tolower))
corp <- tm_map(corp, removeNumbers)
corp <- tm_map(corp, function(x)removeWords(x,stopwords()))
corp <- tm_map(corp, function(x)removeWords(x,c('one', 'two', 'three', 'four', 'five')))

find = c('tumours', 'tumors', 'tumour',
         'signatures',
         'genes',
         'subtypes',
         'levels',
         'lines',
         'models',
         'subpopulations',
         'tumorigeneis')
repl = c('tumor', 'tumor', 'tumor',
         'signature',
         'gene',
         'subtype',
         'level',
         'line',
         'model',
         'subpopulation',
         'tumorigenesis')
#corp = tm_map(corp, stemDocument)
synonyms <- function(x){
  for(i in 1:length(find)){
    x = gsub(find[i], repl[i], x)
  }
  return(PlainTextDocument(x))
}
corp <- tm_map(corp, synonyms)

term.matrix <- TermDocumentMatrix(corp)
term.matrix <- as.matrix(term.matrix)
colnames(term.matrix) <- names(clusters)

savePlot = F
if(savePlot) pdf('wordcloud_cell_lines.pdf')
comparison.cloud(term.matrix,max.words=40,random.order=FALSE)

#commonality.cloud(term.matrix,max.words=40,random.order=FALSE)
if(savePlot) dev.off()


