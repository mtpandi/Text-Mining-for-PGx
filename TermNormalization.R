termnorm <- function(inpdf){
  
  
  #genes
  genenorm <- function(df){
    #Genes
    gene_anns <- function(genes){
      
      gene.batches <- split(genes, ceiling(seq_along(genes)/100))
      gene.batches <- pblapply(gene.batches, function(y) paste0(y, collapse = ","))
      
      
      recs <- pblapply(gene.batches, function(y) rentrez::entrez_summary(db = "gene", id = y))
      recs <- recs %>% unlist(recursive = F) %>% unname
      
      
      genesres<-list()
      
      for(i in 1:length(recs)){
        x<-recs[[i]]
        
        symbol<-x['name'] %>% as.character()
        uid<-x['uid'] %>% as.character()
        org<-x['organism'] %>% unname %>% unlist(recursive = F)
        orgid<-org['taxid'] %>% as.character()
        orgname<-org['scientificname'] %>% as.character()
        
        ires<-list(data.table(ID=uid, 
                              TERM=symbol,
                              ORGANISM=orgname,
                              GeneOrganismID=orgid,
                              stringsAsFactors = FALSE))
        
        genesres<-append(genesres, ires)
      }
      
      res<-rbindlist(genesres, fill = TRUE) %>% mutate_if(is.list, as.character)
      
      return(res)
      
    }
    
    gg <- df %>% subset(type=='Gene') %>% select(identifier) %>% unlist %>% unique
    genedf <- gene_anns(gg) %>% as.data.frame(stringsAsFactors = FALSE) %>% mutate(type = 'Gene')
    genedf <- genedf %>% subset(GeneOrganismID=='9606') %>% select(ID, TERM, type) %>% unique()
    
    return(genedf)
    
  }
  
  gnormdf <- genenorm(inpdf)
  
  
  #chemicals
  inpdf$identifier <- pblapply(inpdf$identifier, function(x) gsub("MESH:", "", x))
  inpdf <- inpdf %>%  mutate_if(is.list, as.character)
  
  chemanns <- function(df){
    
    # for the mesh ids we have we must get the uids
    
    #get unique mesh ids and create a df with the mismatching ones
    # 2 cols: original(identifier), toquery(New)
    cc <- df %>% subset(type=='Chemical') %>% select(identifier) %>% 
      subset(!is.na(identifier))  %>% unique 
    
    meshs <- read.delim('concept_ID_updates.txt', sep='\t', header=F, stringsAsFactors = F)
    colnames(meshs) <- c('Old', 'New')
    
    meshs$Old <- pblapply(meshs$Old, function(x) gsub("MESH:", "", x)) 
    meshs$New <- pblapply(meshs$New, function(x) gsub("MESH:", "", x)) 
    meshs <- meshs %>%  mutate_if(is.list, as.character)
    
    
    cc <- left_join(cc, meshs, by=c('identifier'='Old'))
    cc <- cc %>% mutate(Query = ifelse(is.na(New), identifier, New)) %>%
      select(-New)
    
    #ctd has up-to-date mesh terms -> filter the Query column
    
    ok <- read.csv("FilterChemicals.csv", header = T, sep=',', stringsAsFactors = F)
    meshdb_ids <- ok %>% unlist %>% unique %>% as.character()
    
    cc <- cc %>% subset(Query %in% meshdb_ids)
    
    uidfun <- function(chem){
      
      res<-reutils::uid(reutils::esearch(chem, db = "mesh", rettype = "uilist"))
      
      if(!is.na(res)){
        ncbiid <- res %>% as.character()
      }else{
        ncbiid <- 'NotFound'
      }
      
      return(ncbiid)
    }
    
    cc$NCBIuid <- pblapply(cc$Query, uidfun) 
    cc <- cc %>% mutate_if(is.list, as.character) %>% subset(NCBIuid!='NotFound')
    
    chemnorm <- function(ncbiid){
      
      ncbiid.batches <- split(ncbiid, ceiling(seq_along(ncbiid)/100))
      ncbiid.batches <- pblapply(ncbiid.batches, function(y) paste0(y, collapse = ","))
      
      
      recs <- pblapply(ncbiid.batches, function(y) rentrez::entrez_summary(db = "mesh", id = y))
      recs <- recs %>% unlist(recursive = F) %>% unname
      
      cheml <- list()
      
      for(i in 1:length(recs)){
        x<-recs[[i]]
        if((length(x)!=0)&('ds_meshterms' %in% names(x))&('ds_meshui' %in% names(x))){
          meshterm<-x['ds_meshterms'][[1]][1] %>% as.character
          meshid<-x['ds_meshui'] %>% as.character
        }else{
          meshterm<- 'None'
          meshid<- 'None'
        }
        
        
        cheml[[i]] <- data.frame(ID=meshid, 
                                 TERM=meshterm,
                                 stringsAsFactors = FALSE)
        
        
        
      }
      
      cheml <- cheml  %>% bind_rows() %>% mutate(type = 'Chemical') %>% subset(!ID=='None')
      cheml <- cheml %>% distinct(ID, .keep_all=T)
      return(cheml)
    }
    
    
    
    cnormdf <- chemnorm(cc[['NCBIuid']])
    cnormdf <- left_join(cnormdf, cc[,c('identifier', 'Query')], by=c('ID'='Query'))
    cnormdf <- cnormdf %>% select(-ID) 
    names(cnormdf)[names(cnormdf)=='identifier'] <-'ID'
    
    return(cnormdf)
  }
  
  chemnormdf <- chemanns(inpdf)

  
  #merge the dfs
  norms <- bind_rows(gnormdf, chemnormdf)
  
  #combine annotations with data
  inpdf <- left_join(inpdf, norms, by=c('identifier'='ID', 'type'))
  inpdf <- inpdf %>% mutate(TERM=ifelse(type=='Mutation', identifier, TERM))

  #filter mutations (keep only with rsIDs)
  inpdf <- inpdf %>% subset(!((!grepl('^rs[0-9]+', inpdf$identifier))&(type=='Mutation')))
  
 

  return(inpdf)
}

