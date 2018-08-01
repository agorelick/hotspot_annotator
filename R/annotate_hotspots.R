##' annotate_hotspots
##' @export 
annotate_hotspots <- function(maf,min_indel_overlap=0.5) {
    ## this function will annotate hotspot SNP, Splice Site, and Indels in the MAF with Taylor Lab published hotspots
    require(data.table)
    data(hotspots)

    hotspots$class <- as.character(hotspots$class)
    hotspots$type <- as.character(hotspots$type)
    indel_hotspots <- hotspots[hotspots$type=='indel',c('Hugo_Symbol','hotspot','start','end'),with=F]
    setkey(indel_hotspots,'Hugo_Symbol','start','end')
    maf$mutation <- paste(maf$Hugo_Symbol,gsub('p[.]','',maf$HGVSp_Short))
    maf_original <- maf
    maf <- maf[maf$HGVSp_Short!='',]

    ## snp hotspots
    d_snp <- maf[maf$HGVSp_Short!='.' & maf$Variant_Classification %in% c('Missense_Mutation','Nonsense_Mutation'),c('Hugo_Symbol','HGVSp_Short'),with=F]
    if(nrow(d_snp)==0) {
        message('No missense or nonsense mutations found!')
        d_snp <- as.data.table(matrix(nrow=0,ncol=6))
        names(d_snp) <- c('Hugo_Symbol','mutation','mutation.orig','class','type','hotspot')
    } else {
        message('Annotating Missense/Nonsense hotspots ...')
        d_snp$mutation.orig <- paste(d_snp$Hugo_Symbol,gsub('p[.]','',d_snp$HGVSp_Short))
        d_snp$mutation <- substr(d_snp$mutation.orig,1,(nchar(d_snp$mutation.orig)-1))
        d_snp <- merge(d_snp,hotspots[,c('hotspot','class','type'),with=F],by.x='mutation',by.y='hotspot',all.x=T)
        d_snp[,HGVSp_Short:=NULL]
        d_snp$class[is.na(d_snp$class)] <- ''
        d_snp$type[is.na(d_snp$type)] <- ''
        d_snp$hotspot <- d_snp$mutation; d_snp$hotspot[d_snp$class==''] <- ''
    }

    ## splice site hotspots
    d_ss <- maf[maf$HGVSp_Short!='.' & maf$Variant_Classification %in% 'Splice_Site',c('Hugo_Symbol','HGVSp_Short'),with=F]
    if(nrow(d_ss)==0) {
        message('No splice site mutations found!')
        d_ss <- as.data.table(matrix(nrow=0,ncol=6))
        names(d_ss) <- c('Hugo_Symbol','mutation','mutation.orig','class','type','hotspot')
    } else {
        message('Annotating Splice Site hotspots ...')
        d_ss$mutation.orig <- paste(d_ss$Hugo_Symbol,gsub('p[.]','',d_ss$HGVSp_Short))
        d_ss$mutation <- d_ss$mutation.orig
        d_ss <- merge(d_ss,hotspots[,c('hotspot','class','type'),with=F],by.x='mutation',by.y='hotspot',all.x=T)
        d_ss[,HGVSp_Short:=NULL]
        d_ss$class[is.na(d_ss$class)] <- ''
        d_ss$type[is.na(d_ss$type)] <- ''
        d_ss$hotspot <- d_ss$mutation; d_ss$hotspot[d_ss$class==''] <- ''
    }

    ## annotate indel hotspots
    d_indel <- maf[maf$HGVSp_Short!='.' & maf$Variant_Classification %in% c('In_Frame_Ins','In_Frame_Del'),c('Hugo_Symbol','HGVSp_Short'),with=F]
    if(nrow(d_indel)==0) {
        message('No in-frame INDELs found!')
        d_indel <- as.data.table(matrix(nrow=0,ncol=6))
        names(d_indel) <- c('Hugo_Symbol','mutation','mutation.orig','class','type','hotspot')
    } else {
        message('Annotating in-frame INDEL hotspots ...')
        d_indel$mutation.orig <- paste(d_indel$Hugo_Symbol,gsub('p[.]','',d_indel$HGVSp_Short))
        ## change p.T454_Q455>Q -> p.T454_Q455 to extract the start/end amino acid positions
        f <- function(x) as.list(x[[1]])
        d_indel$HGVSp_Short <- rbindlist(lapply(strsplit(d_indel$HGVSp_Short,'>'),f))[[1]]
        indels <- gsub('p[.]','',d_indel$HGVSp_Short)
        prune_delins <- function(s) unlist(strsplit(s,'delins'))[1]
        prune_del <- function(s) unlist(strsplit(s,'del'))[1]
        prune_ins <- function(s) unlist(strsplit(s,'ins'))[1]
        prune_dup <- function(s) unlist(strsplit(s,'dup'))[1]
        f <- function(x) {
            if(length(x)==1) x <- c(x,x)
            as.list(x)
        }
        l <- lapply(strsplit(indels,'_'),f)
        indels <- cbind(Hugo_Symbol=d_indel$Hugo_Symbol,HGVSp_Short=indels,rbindlist(l),mutation.orig=d_indel$mutation.orig)
        indels$mutation <- paste(indels$Hugo_Symbol,indels$HGVSp_Short)
        indels$V1 <- sapply(indels$V1,prune_delins,USE.NAMES=F)
        indels$V1 <- sapply(indels$V1,prune_del,USE.NAMES=F)
        indels$V1 <- sapply(indels$V1,prune_ins,USE.NAMES=F)
        indels$V1 <- sapply(indels$V1,prune_dup,USE.NAMES=F)
        indels$V2 <- sapply(indels$V2,prune_delins,USE.NAMES=F)
        indels$V2 <- sapply(indels$V2,prune_del,USE.NAMES=F)
        indels$V2 <- sapply(indels$V2,prune_ins,USE.NAMES=F)
        indels$V2 <- sapply(indels$V2,prune_dup,USE.NAMES=F)
        indels[,HGVSp_Short:=NULL]
        d_indel <- indels; rm(indels)
        d_indel$start <- as.integer(substr(d_indel$V1,2,nchar(d_indel$V1)))
        d_indel$end <- as.integer(substr(d_indel$V2,2,nchar(d_indel$V2)))
        d_indel[,c('V1','V2'):=NULL]
        d_indel$start <- as.integer(d_indel$start)
        d_indel$end <- as.integer(d_indel$end)
        setkey(d_indel,'Hugo_Symbol','start','end')
        hits <- foverlaps(d_indel,indel_hotspots)
        hits$start <- as.integer(hits$start)
        hits$end <- as.integer(hits$end)
        hits$i.start <- as.integer(hits$i.start)
        hits$i.end <- as.integer(hits$i.end)

        ## for each indel hotspot, get % overlap
        get_pct_overlap <- function(i,d) {
            d <- d[i,]
            indel_region <- seq(d$i.start,d$i.end)
            if(is.na(d$start)) d$start <- -2
            if(is.na(d$end)) d$end <- -1
            hotspot_region <- seq(d$start,d$end)
            overlap <- sum(indel_region %in% hotspot_region) / length(indel_region)
            if(overlap==0) overlap <- as.numeric(NA)
            d$overlap <- overlap
            d
        }
        l <- lapply(1:nrow(hits),get_pct_overlap,hits)
        hits <- rbindlist(l)
        hits <- hits[,c('Hugo_Symbol','mutation','mutation.orig','overlap','hotspot'),with=F]
        hits$class <- ''; hits$type <- '';
        hits$hotspot[hits$overlap < min_indel_overlap | is.na(hits$overlap)] <- ''
        hits$class[hits$overlap >= min_indel_overlap] <- '2d'
        hits$type[hits$overlap >= min_indel_overlap] <- 'indel'
        hits[,overlap:=NULL]
        d_indel <- hits[,c('Hugo_Symbol','mutation','mutation.orig','class','type','hotspot'),with=F]
    }

    ## merge annotated hotspots
    result <- rbind(d_snp,d_ss,d_indel,fill=T)
    result[,c('Hugo_Symbol','mutation'):=NULL]
    result <- result[!duplicated(result$mutation.orig),]

    names(result) <- c('mutation','hotspot_class','hotspot_type','hotspot')
    maf_annotated <- merge(maf_original,result,by='mutation',all.x=T)
    maf_annotated$hotspot_class[is.na(maf_annotated$hotspot_class)] <- ''
    maf_annotated$hotspot_type[is.na(maf_annotated$hotspot_type)] <- ''
    maf_annotated$hotspot[is.na(maf_annotated$hotspot)] <- ''
    maf_annotated[,mutation:=NULL]
    maf_annotated
} 





