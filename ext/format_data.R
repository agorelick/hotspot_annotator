library(data.table)

setwd('~/luna/repos/hotspot_annotator')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# format raw data from 3d hotspots
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d <- fread('data/hotspots_3d_11k.txt')

parse_string <- function(str) { 
    str_data <- data.table(str=unlist(strsplit(str,';')))
    str_data$str[grepl('[(]',str_data$str)==F] <- paste0(str_data$str[grepl('[(]',str_data$str)==F],'(NA)')

    ## parse mutations
    str_data$str2 <- gsub('[)]','',str_data$str)
    info <- rbindlist(lapply(strsplit(str_data$str2,'[(]'),as.list))
    names(info) <- c('val1','val2')
    str_data <- cbind(str_data,info)
    str_data[,c('val1','val2'),with=F]
}


## expand cluster to get the residues mutated in each
expand_cluster <- function(d) {
    out <- d[1,]
    out[,residues_mutation_numbers:=NULL]
    mutation_info <- parse_string(d$residues_mutation_numbers[1])
    names(mutation_info) <- c('mutation','samples')
    out <- cbind(out,mutation_info)
    out    
}
dd <- d[,expand_cluster(.SD),by=cluster_id]

## format and rename fields
setnames(dd,'gene','Hugo_Symbol')
setnames(dd,'p_val','hotspot_3d.pval')
setnames(dd,'cluster_id','hotspot_3d.cluster_id')
setnames(dd,'pdb_chain_p_val','hotspot_3d.pdb_chain_p_val')
setnames(dd,'total_mutation_number','hotspot_3d.total_mutation_number')
setnames(dd,'number_residues','hotspot_3d.number_residues')
setnames(dd,'mutation','hotspot')
dd$class <- 'hotspot_3d'
dd$Reference_Amino_Acid <- substr(dd$hotspot,1,1)
dd$Amino_Acid_Position <- as.integer(substr(dd$hotspot,2,nchar(dd$hotspot)))
dd$hotspot <- paste0(dd$Hugo_Symbol,' ',dd$hotspot)
dd$start_position <- dd$Amino_Acid_Position; dd$end_position <- dd$Amino_Acid_Position+1
## order fields for output
front.fields <- c('Hugo_Symbol','Reference_Amino_Acid','Amino_Acid_Position','hotspot','class','start_position','end_position') 
hotspots_3d <- dd[,front.fields,with=F]
hotspots_3d <- hotspots_3d[!duplicated(hotspots_3d$hotspot),]
hotspots_3d$type <- 'snp'



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# format raw data from single-codon mutations in hotspots v2
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d <- fread('data/hotspots_v2.txt')
d[,Allele_Freq_Rank:=NULL]

dd <- d[,c('Hugo_Symbol','Amino_Acid_Position','Reference_Amino_Acid'),with=F]
dd$Reference_Amino_Acid <- substr(dd$Reference_Amino_Acid,1,1)
dd$Reference_Amino_Acid[grepl('splice',dd$Amino_Acid_Position)] <- ''
dd$hotspot <- paste0(dd$Hugo_Symbol,' ',dd$Reference_Amino_Acid,dd$Amino_Acid_Position)
dd$class <- 'hotspot_v2'
dd$start_position <- as.integer(gsub('X','',gsub('_splice','',dd$Amino_Acid_Position)))
dd$end_position <- dd$start_position + 1
front.fields <- c('Hugo_Symbol','Reference_Amino_Acid','Amino_Acid_Position','hotspot','class','start_position','end_position') 
hotspots_v2 <- dd[,front.fields,with=F]
hotspots_v2 <- hotspots_v2[!duplicated(hotspots_v2$hotspot),]
hotspots_v2$type <- 'snp'
hotspots_v2$type[grep('splice',hotspots_v2$hotspot)] <- 'splice_site'
hotspots_v2 <- hotspots_v2[!duplicated(hotspots_v2$hotspot),]




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# format raw data from in-frame indels in hotspots v2
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d <- fread('data/hotspots_v2_inframe_indels.txt')
dd <- d[,c('Hugo_Symbol','Amino_Acid_Position'),with=F]
dd$Reference_Amino_Acid <- ''
dd$hotspot <- paste0(dd$Hugo_Symbol,' ',dd$Reference_Amino_Acid,dd$Amino_Acid_Position)

f <- function(s) {
    if(length(s)==1) s <- c(s,s)
    as.list(s)
}
to.merge <- rbindlist(lapply(strsplit(dd$Amino_Acid_Position,'-'),f))
names(to.merge) <- c('start_position','end_position')
dd <- cbind(dd,to.merge)
dd$class <- 'hotspot_v2'
front.fields <- c('Hugo_Symbol','Reference_Amino_Acid','Amino_Acid_Position','hotspot','class','start_position','end_position') 
hotspots_v2_indel <- dd[,front.fields,with=F]
hotspots_v2_indel$type <- 'indel'
hotspots_v2_indel <- hotspots_v2_indel[!duplicated(hotspots_v2_indel$hotspot),]
hotspots_v2_indel$start_position <- as.integer(hotspots_v2_indel$start_position)
hotspots_v2_indel$end_position <- as.integer(hotspots_v2_indel$end_position)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# format raw data from hotspots V1
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d <- fread('data/hotspots_11k.txt')
dd <- d[!is.na(d$qval)]
setnames(dd,'Hugo Symbol','Hugo_Symbol')
setnames(dd,'codon','hotspot')
ss <- dd[grepl('splice',dd$hotspot),]
ss$Reference_Amino_Acid <- ''
ss$start_position <- as.integer(substr(gsub('_splice','',ss$hotspot),2,nchar(ss$hotspot)))
ss$Amino_Acid_Position <- paste0('X',substr(ss$hotspot,2,nchar(ss$hotspot)))
ss <- ss[,c('Hugo_Symbol','Reference_Amino_Acid','Amino_Acid_Position','hotspot','start_position'),with=F]
ss$hotspot <- paste(ss$Hugo_Symbol,ss$Amino_Acid_Position)
dd <- dd[!grepl('splice',dd$hotspot),]
dd$Reference_Amino_Acid <- substr(dd$hotspot,1,1)
dd$Amino_Acid_Position <- as.integer(substr(dd$hotspot,2,nchar(dd$hotspot)))
dd$start_position <- dd$Amino_Acid_Position
dd$hotspot <- paste0(dd$Hugo_Symbol,' ',dd$Reference_Amino_Acid,dd$Amino_Acid_Position)
dd <- dd[,c('Hugo_Symbol','Reference_Amino_Acid','Amino_Acid_Position','hotspot','start_position'),with=F]
hotspots_v1 <- rbind(dd,ss)
hotspots_v1$class <- 'hotspot_v1'
hotspots_v1$type <- 'snp'
hotspots_v1$type[grep('splice',hotspots_v1$hotspot)] <- 'splice_site'
hotspots_v1 <- hotspots_v1[!duplicated(hotspots_v1$hotspot),]
hotspots_v1$end_position <- hotspots_v1$start_position+1
front.fields <- c('Hugo_Symbol','Reference_Amino_Acid','Amino_Acid_Position','hotspot','class','start_position','end_position','type') 
hotspots_v1 <- hotspots_v1[,front.fields,with=F]



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# merge results
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

all_hotspots <- rbind(hotspots_v1,hotspots_v2,hotspots_v2_indel,hotspots_3d)

## collapse SNP/Splice_Site hotspots
snp <- all_hotspots[all_hotspots$type %in% c('snp','splice_site'),]
snp$class <- factor(gsub('hotspot_','',snp$class),levels=c('v1','v2','3d'))
summarize_hotspot <- function(d) {
    out <- d[1,]
    class <- paste(sort(unique(d$class)),collapse=',')
    out$class <- class
    out
}
snp <- snp[,summarize_hotspot(.SD),by=hotspot]
indel <- all_hotspots[all_hotspots$type=='indel',]
indel$class <- 'v2'
indel$class <- factor(indel$class,levels=c('v1','v2','3d'))
result <- rbind(snp,indel)
setnames(result,'start_position','start')
setnames(result,'end_position','end')
result <- result[order(result$Hugo_Symbol,result$start,result$end),]
write.tsv(result,'data/hotspots_merged.tsv')
save(result,file='data/hotspots.rda')



