library('STRINGdb')
library('tidyverse')
library('RCurl')
library('stringr')
library(ggplot2)
setwd("~/travail Bioinfo/projet matthieu/idh/projet")



sdb = STRINGdb$new(version='10',species=7955,score_threshold=400,input_directory='data')


#Plusieurs fonctions de regroupements
proteins = sdb$get_proteins()
neighbors = sdb$get_neighbors('7955.ENSDARP00000051548')
interaction=sdb$get_interactions(neighbors)

#regroupement en fonction du text mining
interaction_text_mining = interaction[order(interaction[,14],decreasing=T),]
interaction_text_mining= interaction_text_mining
set_text_mining_information = unique(c(interaction_text_mining$from,interaction_text_mining$to))
prot_text_mining = proteins[proteins$protein_external_id %in% set_text_mining_information,]

sink("set_text_mining_500_information_p53.txt")
cat(paste0(prot_text_mining$preferred_name[1:500],"\t"),sep = "")
sink()

sink("set_text_mining_50_information_p53.txt")
cat(paste0(prot_text_mining$preferred_name[1:50],"\t"),sep = "")
sink()

#regroupement en fonction de la coexpression
interaction_coexpression= interaction[order(interaction[,8],decreasing=T),]
interaction_coexpression= interaction_coexpression
set_coexpression_information = unique(c(interaction_coexpression$from,interaction_coexpression$to))
prot_coexpression = proteins[proteins$protein_external_id %in% set_coexpression_information,]

sink("set_coexpression_500_information_p53.txt")
cat(paste0(prot_coexpression$preferred_name[1:500],"\t"),sep = "")
sink()

sink("set_coexpression_50_information_p53.txt")
cat(paste0(prot_coexpression$preferred_name[1:50],"\t"),sep = "")
sink()


#utilisation de l'information globale pour regrouper
interaction=sdb$get_interactions(neighbors)
interaction_all= interaction[order(interaction[,16],decreasing=T),]
interaction_all= interaction_all
set_all_information = unique(c(interaction_all$from,interaction_all$to))
prot_all = proteins[proteins$protein_external_id %in% set_all_information,]

sink("set_all_500_information_p53.txt")
cat(paste0(prot_all$preferred_name[1:500],"\t"),sep = "")
sink()


sink("set_all_50_information_p53.txt")
cat(paste0(prot_all$preferred_name[1:50],"\t"),sep = "")
sink()

#sdb$plot_network(interaction$from)

set_500 = sample(1:length(proteins$preferred_name),500,replace = F)
set_50 = sample(1:length(proteins$preferred_name),50,replace = F)

for(i in c(1:10)){
  
  set_500 = sample(1:length(proteins$preferred_name),500,replace = F)
  set_50 = sample(1:length(proteins$preferred_name),50,replace = F)
  
  sink(paste0("set_aleatoire_500_information_p53_",i,".txt"))
  cat(paste0(proteins$preferred_name[set_500],"\t"),sep = "")
  sink()
  
  sink(paste0("set_aleatoire_50_information_p53_",i,".txt"))
  cat(paste0(proteins$preferred_name[set_50],"\t"),sep = "")
  sink()
}

for(i in c(1:10)){
  
  sink('files.txt')
  cat(paste0("set_aleatoire_500_information_p53_",i,".txt"),sep = " ")
  cat(paste0("set_aleatoire_50_information_p53_",i,".txt"),sep = " ")
  sink()
}

set = unique(neighbors[sample(1:1825, 1825, replace=FALSE)])
prot = proteins[proteins$protein_external_id %in% set,]

write.csv2(prot$preferred_name,file='set_gene_p53.txt',quote = F)
sink("set_all_information_p53.txt")
cat(paste0(prot$preferred_name,"\t"))
sink()


# Stockage des IDs des pathways
all_pathways = read.table(file = "pathways_uniq.txt",colClasses = "factor")
all_pathways = unlist(all_pathways$V1)

# Creation de la liste des genes de pathways:
pathway.genes = list()
pathway.genes <- vector(mode='list', length=length(all_pathways))
names(pathway.genes) <- all_pathways

# Remplissage de pathway.genes avec tous les genes pour chaque pathway:
for (p in (1:(length(all_pathways)))){
  temp <- sdb$get_term_proteins(all_pathways[[p]])$preferred_name
  pathway.genes[[p]] <- temp
}

# Ecriture dans un fichier
sink("zebrafish_pathway_sets.txt")
cat(paste0("# format: sets\n# version: 1.2\n# strain: Danio rerio\n# date: ",Sys.Date(),"\n# comment: All Danio rerio genes in pathways\n"))
for (i in 1:length(pathway.genes)){
  cat(((names(pathway.genes)[i])),sep = "")
  cat(paste0("\t",(str_trim(getURL(paste("http://togows.dbcls.jp/entry/pathway/dre",names(pathway.genes)[i],"/name", sep="")))),sep = ""))
  cat(paste0("\t",(pathway.genes[[i]])),sep = "")
  cat("\n")
}
sink()


#test

#les regroupements trouvés ici sont les mêmes que ceux que j'ai trouvé sur la GO etc avec les mesures binomiales par exemple
enrichment_Process=sdb$get_enrichment(prot$protein_external_id, category = "Process", methodMT = "fdr", iea = TRUE, minScore=NULL )
enrichment_Component=sdb$get_enrichment(prot$protein_external_id, category = "Component", methodMT = "fdr", iea = TRUE, minScore=NULL )
enrichment_Function=sdb$get_enrichment(prot$protein_external_id, category = "Function", methodMT = "fdr", iea = TRUE, minScore=NULL )
enrichment_KEGG=sdb$get_enrichment(prot$protein_external_id, category = "KEGG", methodMT = "fdr", iea = TRUE, minScore=NULL )
enrichment_Pfam=sdb$get_enrichment(prot$protein_external_id, category = "Pfam", methodMT = "fdr", iea = TRUE, minScore=NULL )
enrichment_InterPro=sdb$get_enrichment(prot$protein_external_id, category = "InterPro", methodMT = "fdr", iea = TRUE, minScore=NULL )



#histogramme

distribution = read.csv(file="distribution.txt")
hist(log(distribution$X0),main = 'distribution du nombre de Gene Product par GoTerm dans la GO')
ggplot(data=distribution, aes(log(distribution$X0))) + 
  geom_histogram(binwidth = 0.1, color = 'black', fill = '#099DD9')+
  scale_x_continuous(breaks = 1:10)+
  xlab(label = 'log(distribution des GoTerm ==> Gene Products)')+
  ylab(label= 'count(GoTerm)') +
  ggtitle('distribution du nombre de Gene Product par GoTerm dans la GO',)
  
