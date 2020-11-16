#This program is designed to scrape genbank for accessions of symbiotic microbial species. It 
# performs a directed search on the host and isolation source fields using host species as input. 
# It will return any records with matches then can be later queried to get sequence data for downstream 
# analysis. 
# If a record does not have a match in the host or isolation source fields, it will the entire record and return 
# any records with matches. Matches there will have to be manually assessed as to whether the sequence data came from 
# a host of interest. 

#Note, you will need an entrez key from NCBI to authorize the frequency of NCBI queries that occur in this 
# script. This key is free and generated instantly at https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/

#The host species list must be a dataframe or matrix with genus in first column and species in the second column

setwd("~/Box/Nematode_macroevolution/Macroevolution_species_interactions/")
library(rentrez)
library(stringdist)
library(stringr)
"(Bradyrhizobium[Organism]) AND (16s rRNA[Title] OR 16s rRNA[Gene] OR 16s ribosomal RNA[All Fields])"
symbiosis_scraper <- function(entrez_key, database, search_term, host_species){
  #set the key
  set_entrez_key(entrez_key)
  #make the search
  r_search <- entrez_search(db=database, term=search_term)
  r_search <- entrez_search(db=database, term=search_term, retmax = r_search$count)
  #Enter the host data that will be quried
  hosts <- host_species
  #set the empty vector 
  nat_host <- vector()
  isolation_source <- vector()
  GI <- vector()
  acc_manual_rev <- vector()
  #Set "state-counters" - these will be used to build the output tables 
  a <- 1
  j <- 1
  k <- 1
  l <- 1
  m <- 1
  n <- 1
  o <- 1
  n_records <- r_search$count
  for(i in 1:n_records){
    if(i %% 25 == 0){
      #This is to make sure that we meet NCBIs 10 querries per second regulation
      Sys.sleep(1)
    }
    #pull down the record
    rec <- entrez_fetch(db = "nucleotide", id = r_search$ids[i], rettype = "native")
    #If record is not genome sized, parse it...
    if(object.size(rec) > 100000){
      genome_rec <- r_search$ids[i]
      if(a > 2){
        genomes <- c(genomes, genome_rec)
      }else{
        genomes <- genome_rec
        a <-2
      }
    }else{
      if(grepl("nat-host", rec) == F & grepl("isolation-source", rec) == F){
        GI <- r_search$ids[i]
        nat_host <- NA
        isolation_source <- NA
        acc_manual_rev <- "Nucleotide"
        #If there is a match in host, we parse the record to put into the host vector, if no source record we put NA in source vector
      }else if(grepl("nat-host", rec) == T & grepl("isolation-source", rec) == F){
        s <- unlist(strsplit(rec, split = "\n"))
        x <- s[grep("nat-host",s)+1]
        nat_host <- gsub('"', '', regmatches(x, gregexpr('"([^"]*)"', x))[[1]])
        isolation_source <- NA
        GI <- r_search$ids[i]
        acc_manual_rev <- "Nucleotide"
        #If there no match in host, we put an NA in the host vector, but if there is a match in isolation source we parse the record and put it in the source vector. 
      }else if(grepl("nat-host", rec) == F & grepl("isolation-source", rec) == T){
        s <- unlist(strsplit(rec, split = "\n"))
        x <- s[grep("isolation-source",s)+1]
        isolation_source <- gsub('"', '', regmatches(x, gregexpr('"([^"]*)"', x))[[1]])
        nat_host <- NA
        GI <- r_search$ids[i]
        acc_manual_rev <- "Nucleotide"
        # If there are matches in both, they both go into their respective vectors
      }else if(grepl("nat-host", rec) == T & grepl("isolation-source", rec) == T){
        s <- unlist(strsplit(rec, split = "\n"))
        x <- s[grep("isolation-source",s)+1]
        isolation_source <- gsub('"', '', regmatches(x, gregexpr('"([^"]*)"', x))[[1]])
        x <- s[grep("nat-host",s)+1]
        nat_host <- gsub('"', '', regmatches(x, gregexpr('"([^"]*)"', x))[[1]])
        GI <- r_search$ids[i]
        acc_manual_rev <- "Nucleotide"
      }
      #Ok, now all those respective vectors go inot a data frame that has the shape:
      #           NCBI-IDs   Host    isolation source    record type 
      # Record...
      search_dat <- data.frame(Ids = GI, host = nat_host, isolation_source = isolation_source, rec_type = acc_manual_rev, stringsAsFactors = F)
      #Coerce any NAs into character 
      search_dat$host <- as.character(search_dat$host)
      search_dat$isolation_source <- as.character(search_dat$isolation_source)
      if(is.na(nat_host) & is.na(isolation_source)){
        #IF there is nothing in these expected fields 
        deep_rec_search <- !is.na(str_locate(pattern = Legumes[,1], rec)[,1])
        deep_rec_search <- data.frame(Ids = rep(GI, sum(deep_rec_search)), leg_genus = Legumes[deep_rec_search,1], leg_species = Legumes[deep_rec_search,2])
        if(nrow(deep_rec_search)> 0){
          #search the species name 
          deep_rec_search_spec <- !is.na(str_locate(pattern = deep_rec_search[,3], rec)[,1])
          if(deep_rec_search_spec > 0){
            if(m > 1){
              species_match_deep <- rbind(species_match_deep, deep_rec_search[deep_rec_search_spec,])
            }else{
              species_match_deep <- deep_rec_search[deep_rec_search_spec,]
              m <- 2
            }
            #if no species match... put in the genus match
          }else{
            if(n>1){
              genus_match_deep <- rbind(genus_match_deep, deep_rec_search)
            }else{
              genus_match_deep <- deep_rec_search
              n <- 2
            }
          }
        }else{
          if(o >1){
            no_match_deep <- rbind(no_match_deep, search_dat)
          }else{
            no_match_deep <- rbind(search_dat)
            o <- 2
          }
        }
      }else{
        #split up the host field and search it against the legumes 
        split_host <- unlist(strsplit(search_dat[,2], split = " "))
        match_ind <- amatch(split_host, Legumes[,1], maxDist = 1, nomatch = 0)
        #If there's a match, we'll search deeper to pull out the legume. 
        if(sum(match_ind) > 0){
          leg_match <- Legumes[match_ind,]
          leg_match <- subset(Legumes, Legumes[,1] == leg_match[1])
          match_ind_2 <- amatch(tolower(split_host), leg_match[,2], maxDist = 1, nomatch = 0)
          #If there is a species match, we will put it in its own table.
          if(sum(match_ind_2) > 0){
            #Put species matched into their own table 
            if(j > 1){
              species_match <- rbind(species_match, search_dat)
            }else{
              species_match <- search_dat
              j <- 2
            }
          }else{
            #If no species level match, then it goes in the "genus" table for further consideration
            if(k > 1){
              genus_match <- rbind(genus_match, search_dat)
            }else{
              genus_match <- search_dat
              k <- 2
            }
          }
        }else{
          #We search the isolation source field for matches, if there aren't any in host
          split_source <- unlist(strsplit(search_dat[,3], split = " "))
          match_ind <- amatch(split_source, Legumes[,1], maxDist = 1, nomatch = 0)
          if(sum(match_ind) > 0){
            leg_match <- Legumes[match_ind,]
            leg_match <- subset(Legumes, Legumes[,1] == leg_match[1])
            match_ind_2 <- amatch(tolower(split_source), leg_match[,2], maxDist = 1, nomatch = 0)
            if(sum(match_ind_2) > 0){
              if(j > 1){
                species_match <- rbind(species_match, search_dat)
              }else{
                species_match <- search_dat
                j <- 2
              }
            }else{
              if(k > 1){
                genus_match <- rbind(genus_match, search_dat)
              }else{
                genus_match <- search_dat
                k <- 2
              }
            }
          }else{
            #If there is no match in the isolation source, we will do a deep search as above. 
            deep_rec_search <- !is.na(str_locate(pattern = Legumes[,1], rec)[,1])
            deep_rec_search <- data.frame(Ids = rep(GI, sum(deep_rec_search)), leg_genus = Legumes[deep_rec_search,1], leg_species = Legumes[deep_rec_search,2])
            if(nrow(deep_rec_search)> 0){
              #search the species name 
              deep_rec_search_spec <- !is.na(str_locate(pattern = deep_rec_search[,3], rec)[,1])
              if(deep_rec_search_spec > 0){
                if(m > 1){
                  species_match_deep <- rbind(species_match_deep, deep_rec_search[deep_rec_search_spec,])
                }else{
                  species_match_deep <- deep_rec_search[deep_rec_search_spec,]
                  m <- 2
                }
                #if no species match... put in the genus match
              }else{
                if(n>1){
                  genus_match_deep <- rbind(genus_match_deep, deep_rec_search)
                }else{
                  genus_match_deep <- deep_rec_search
                  n <- 2
                }
              }
            }else{
              if(o >1){
                no_match_deep <- rbind(no_match_deep, search_dat)
              }else{
                no_match_deep <- rbind(search_dat)
                o <- 2
              }
            }
          }
        }
      }
    }
    print(i)
  }
  if(a == 1){
    genomes <- NA
  }
  if(j == 1){
    species_match <- NA
  }
  if(k == 1){
    genus_match <- NA
  }
  if(m == 1){
    species_match_deep <- NA 
  }
  if(n == 1){
    genus_match_deep <- NA
  }
  if(o == 1){
    no_match_deep <- NA
  }
  data_lists <- list(species_match, genus_match, species_match_deep, genus_match_deep, no_match_deep, genomes)
}

Legumes <- read.table("../Nema_macro_old/Nemabase/legume_in_tree.txt", stringsAsFactors = F)
Legumes <- matrix(unlist(strsplit(Legumes$V1, split = "_")), nrow = nrow(Legumes), ncol = 2, byrow = T)
Legumes <- unique(Legumes)

tim <- Sys.time()
test <- symbiosis_scraper(entrez_key = "56361f6f407d6e7d2386b0f5f06fdeaebb09", 
                          database = "nucleotide",
                          search_term = "(Bradyrhizobium[Organism]) AND (16s rRNA[Title] OR 16s rRNA[Gene] OR 16s ribosomal RNA[All Fields])", 
                          host_species = Legumes)
Sys.time() - tim

#get an API key from NCBI to make 10 querries per second.
set_entrez_key("56361f6f407d6e7d2386b0f5f06fdeaebb09")

#Check querries against the otus gathered from Harrison et al. 


#Make a search and get all of the Ids - these records will then be queried individually below. 
#Can also get total records returned by a search with r_search$count
r_search <- entrez_search(db="nucleotide", term="(Bradyrhizobium[Organism]) AND (16s rRNA[Title] OR 16s rRNA[Gene] OR 16s ribosomal RNA[All Fields])")
r_search <- entrez_search(db="nucleotide", 
                          term="(Bradyrhizobium[Organism]) AND (16s rRNA[Title] OR 16s rRNA[Gene] or 16s ribosomal RNA[All Fields])", 
                          retmax = r_search$count)
#r_search <- entrez_search(db="nucleotide", term="16s rRNA Rhizobium", retmax = r_search$count)

#Create empty vectors to feed data into as search is conducted
nat_host <- vector()
isolation_source <- vector()
GI <- vector()
acc_manual_rev <- vector()
a <- Sys.time()

#Bring in the list of species that are present in the megaphylogeny
Legumes <- read.table("legume_in_tree.txt", stringsAsFactors = F)
Legumes <- matrix(unlist(strsplit(Legumes$V1, split = "_")), nrow = nrow(Legumes), ncol = 2, byrow = T)
Legumes <- unique(Legumes)
a <- 1
j <- 1
k <- 1
m <- 1
n <- 1
o <- 1
p <- 1
#Iterate through all records in the search 
#n_records <- r_search$count
n_records <- 20
for(i in 1:n_records){
  if(i %% 25 == 0){
    #This is to make sure that we meet NCBIs 10 querries per second regulation
    Sys.sleep(1)
  }
  #pull down the record
  rec <- entrez_fetch(db = "nucleotide", id = r_search$ids[i], rettype = "native")
  #If record is not genome sized, parse it...
  if(object.size(rec) > 100000){
    genome_rec <- r_search$ids[i]
    if(a > 2){
      genomes <- c(genomes, genome_rec)
    }else{
      genomes <- genome_rec
      a <-2
    }
  }else{
    if(grepl("nat-host", rec) == F & grepl("isolation-source", rec) == F){
      GI <- r_search$ids[i]
      nat_host <- NA
      isolation_source <- NA
      acc_manual_rev <- "Nucleotide"
      #If there is a match in host, we parse the record to put into the host vector, if no source record we put NA in source vector
    }else if(grepl("nat-host", rec) == T & grepl("isolation-source", rec) == F){
      s <- unlist(strsplit(rec, split = "\n"))
      x <- s[grep("nat-host",s)+1]
      nat_host <- gsub('"', '', regmatches(x, gregexpr('"([^"]*)"', x))[[1]])
      if(identical(nat_host,character(0))){
        nat_host <- x
      }
      isolation_source <- NA
      GI <- r_search$ids[i]
      acc_manual_rev <- "Nucleotide"
      #If there no match in host, we put an NA in the host vector, but if there is a match in isolation source we parse the record and put it in the source vector. 
    }else if(grepl("nat-host", rec) == F & grepl("isolation-source", rec) == T){
      s <- unlist(strsplit(rec, split = "\n"))
      x <- s[grep("isolation-source",s)+1]
      isolation_source <- gsub('"', '', regmatches(x, gregexpr('"([^"]*)"', x))[[1]])
      if(identical(isolation_source,character(0))){
        isolation_source <- x
      }
      nat_host <- NA
      GI <- r_search$ids[i]
      acc_manual_rev <- "Nucleotide"
      # If there are matches in both, they both go into their respective vectors
    }else if(grepl("nat-host", rec) == T & grepl("isolation-source", rec) == T){
      s <- unlist(strsplit(rec, split = "\n"))
      x <- s[grep("isolation-source",s)+1]
      isolation_source <- gsub('"', '', regmatches(x, gregexpr('"([^"]*)"', x))[[1]])
      if(identical(isolation_source,character(0))){
        isolation_source <- x
      }
      x <- s[grep("nat-host",s)+1]
      nat_host <- gsub('"', '', regmatches(x, gregexpr('"([^"]*)"', x))[[1]])
      if(identical(nat_host,character(0))){
        nat_host <- x
      }
      GI <- r_search$ids[i]
      acc_manual_rev <- "Nucleotide"
  }
    #Ok, now all those respective vectors go inot a data frame that has the shape:
  #           NCBI-IDs   Host    isolation source    record type 
  # Record...
  search_dat <- data.frame(Ids = GI, host = nat_host, isolation_source = isolation_source, rec_type = acc_manual_rev, stringsAsFactors = F)
  #Coerce any NAs into character 
  search_dat$host <- as.character(search_dat$host)
  search_dat$isolation_source <- as.character(search_dat$isolation_source)
  if(is.na(nat_host) & is.na(isolation_source)){
    #IF there is nothing in these expected fields 
    deep_rec_search <- !is.na(str_locate(pattern = Legumes[,1], rec)[,1])
    deep_rec_search <- data.frame(Ids = rep(GI, sum(deep_rec_search)), leg_genus = Legumes[deep_rec_search,1], leg_species = Legumes[deep_rec_search,2])
    if(nrow(deep_rec_search)> 0){
      #search the species name 
      deep_rec_search_spec <- !is.na(str_locate(pattern = deep_rec_search[,3], rec)[,1])
      if(deep_rec_search_spec > 0){
        if(m > 1){
          species_match_deep <- rbind(species_match_deep, deep_rec_search[deep_rec_search_spec,])
        }else{
          species_match_deep <- deep_rec_search[deep_rec_search_spec,]
          m <- 2
        }
        #if no species match... put in the genus match
      }else{
        if(n>1){
          genus_match_deep <- rbind(genus_match_deep, deep_rec_search)
        }else{
          genus_match_deep <- deep_rec_search
          n <- 2
        }
      }
    }else{
      if(o >1){
        no_match_deep <- rbind(no_match_deep, search_dat)
      }else{
        no_match_deep <- rbind(search_dat)
        o <- 2
      }
    }
  }else{
    #split up the host field and search it against the legumes 
    split_host <- unlist(strsplit(search_dat[,2], split = " "))
    split_host <- gsub("[[:punct:]]", "", split_host)
    match_ind <- amatch(split_host, Legumes[,1], maxDist = 1, nomatch = 0)
    #If there's a match, we'll search deeper to pull out the legume. 
    if(sum(match_ind) > 0){
      leg_match <- Legumes[match_ind,]
      leg_match <- subset(Legumes, Legumes[,1] == leg_match[1])
      match_ind_2 <- amatch(tolower(split_host), leg_match[,2], maxDist = 1, nomatch = 0)
      #If there is a species match, we will put it in its own table.
      if(sum(match_ind_2) > 0){
        #Put species matched into their own table 
        if(j > 1){
          species_match <- rbind(species_match, search_dat)
        }else{
          species_match <- search_dat
          j <- 2
        }
      }else{
        #If no species level match, then it goes in the "genus" table for further consideration
        if(k > 1){
          genus_match <- rbind(genus_match, search_dat)
        }else{
          genus_match <- search_dat
          k <- 2
        }
      }
    }else{
      #We search the isolation source field for matches, if there aren't any in host
      split_source <- unlist(strsplit(search_dat[,3], split = " "))
      match_ind <- amatch(split_source, Legumes[,1], maxDist = 1, nomatch = 0)
      if(sum(match_ind) > 0){
        leg_match <- Legumes[match_ind,]
        leg_match <- subset(Legumes, Legumes[,1] == leg_match[1])
        match_ind_2 <- amatch(tolower(split_source), leg_match[,2], maxDist = 1, nomatch = 0)
        if(sum(match_ind_2) > 0){
          if(j > 1){
            species_match <- rbind(species_match, search_dat)
          }else{
            species_match <- search_dat
            j <- 2
          }
        }else{
          if(k > 1){
            genus_match <- rbind(genus_match, search_dat)
          }else{
            genus_match <- search_dat
            k <- 2
          }
        }
      }else{
        #If there is no match in the isolation source, we will do a deep search as above. 
        deep_rec_search <- !is.na(str_locate(pattern = Legumes[,1], rec)[,1])
        deep_rec_search <- data.frame(Ids = rep(GI, sum(deep_rec_search)), leg_genus = Legumes[deep_rec_search,1], leg_species = Legumes[deep_rec_search,2])
        if(nrow(deep_rec_search)> 0){
          #search the species name 
          deep_rec_search_spec <- !is.na(str_locate(pattern = deep_rec_search[,3], rec)[,1])
          if(deep_rec_search_spec > 0){
            if(m > 1){
              species_match_deep <- rbind(species_match_deep, deep_rec_search[deep_rec_search_spec,])
            }else{
              species_match_deep <- deep_rec_search[deep_rec_search_spec,]
              m <- 2
            }
            #if no species match... put in the genus match
          }else{
            if(n>1){
              genus_match_deep <- rbind(genus_match_deep, deep_rec_search)
            }else{
              genus_match_deep <- deep_rec_search
              n <- 2
            }
          }
        }else{
          if(o >1){
            no_match_deep <- rbind(no_match_deep, search_dat)
          }else{
            no_match_deep <- rbind(search_dat)
            o <- 2
          }
        }
      }
    }
  }
  }
}
