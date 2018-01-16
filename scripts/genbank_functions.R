#!/usr/bin/env Rscript

# functions to make the GenBank submission files (fasta and feature table) for CYTB and RAG1

# make genbank fasta for Cytb
gb_format_fasta_cytb <- function(df, gene){
    organism.name <- ifelse(test=df$taxonRank=="species", yes=paste(df$genus, df$specificEpithet), no=paste(df$genus, "sp."))
    id.by <- str_replace_all(string=iconv(df$identifiedBy, to='ASCII//TRANSLIT'), pattern=" \\| ", replacement="; ")
    specimen.code <- ifelse(test=df$basisOfRecord=="MaterialSample", 
        yes=ifelse(test=is.na(df$collectionCode), yes=paste0(df$institutionCode, ":", df$otherCatalogNumbers), no=paste0(df$institutionCode, ":", df$collectionCode, ":", df$otherCatalogNumbers)), #
        no=ifelse(test=is.na(df$collectionCode), yes=paste0(df$institutionCode, ":", df$catalogNumber), no=paste0(df$institutionCode, ":", df$collectionCode, ":", df$catalogNumber)))
    date <- format(as.Date(df$eventDate, format="%Y-%m-%d"), format="%d-%b-%Y")
    country <- paste0(df$country, ": ", iconv(df$stateProvince, to='ASCII//TRANSLIT'), " State, ", iconv(df$waterBody, to='ASCII//TRANSLIT'), " River drainage")
    lat <- ifelse(test=df$decimalLatitude < 0, yes=paste(str_replace_all(df$decimalLatitude, "-", ""), "S"), no=paste(str_replace_all(df$decimalLatitude, "-", ""), "N"))
    lon <- ifelse(test=df$decimalLongitude < 0, yes=paste(str_replace_all(df$decimalLongitude, "-", ""), "W"), no=paste(str_replace_all(df$decimalLongitude, "-", ""), "E"))
    lat.lon <- paste(lat, lon)
    fas <- ifelse(test=df$basisOfRecord=="MaterialSample", #
        yes=paste0(">", df$otherCatalogNumbers, "_", gene, " ", "[organism=", organism.name, "]", " ", "[Bio_material=", specimen.code, "]", " ", "[location=mitochondrion] [mgcode=2]", " ", "[Collection_date=", date, "]", " ", "[Country=", country, "]", " ", "[Lat_Lon=", lat.lon, "]",  " ", "[Identified_by=", id.by, "]"), #
        no=paste0(">", df$otherCatalogNumbers, "_", gene, " ", "[organism=", organism.name, "]", " ", "[Specimen-voucher=", specimen.code, "]", " ", "[location=mitochondrion] [mgcode=2]", " ", "[Collection_date=", date, "]", " ", "[Country=", country, "]", " ", "[Lat_Lon=", lat.lon, "]", " ", "[Identified_by=", id.by, "]"))
    fas <- str_replace_all(string=fas, pattern=" \\[Lat_Lon=NA NA\\]| \\[Country=NA\\]| \\[Collection_date=NA\\]", replacement="")
    fas <- paste(fas, df$cytb, sep="\n")
        return(fas)
}


# make genbank fasta for Rag1
gb_format_fasta_rag1 <- function(df, gene){
    organism.name <- ifelse(test=df$taxonRank=="species", yes=paste(df$genus, df$specificEpithet), no=paste(df$genus, "sp."))
    id.by <- str_replace_all(string=iconv(df$identifiedBy, to='ASCII//TRANSLIT'), pattern=" \\| ", replacement="; ")
    specimen.code <- ifelse(test=df$basisOfRecord=="MaterialSample", 
        yes=ifelse(test=is.na(df$collectionCode), yes=paste0(df$institutionCode, ":", df$otherCatalogNumbers), no=paste0(df$institutionCode, ":", df$collectionCode, ":", df$otherCatalogNumbers)), #
        no=ifelse(test=is.na(df$collectionCode), yes=paste0(df$institutionCode, ":", df$catalogNumber), no=paste0(df$institutionCode, ":", df$collectionCode, ":", df$catalogNumber)))
    date <- format(as.Date(df$eventDate, format="%Y-%m-%d"), format="%d-%b-%Y")
    country <- paste0(df$country, ": ", iconv(df$stateProvince, to='ASCII//TRANSLIT'), " State, ", iconv(df$waterBody, to='ASCII//TRANSLIT'), " River drainage")
    lat <- ifelse(test=df$decimalLatitude < 0, yes=paste(str_replace_all(df$decimalLatitude, "-", ""), "S"), no=paste(str_replace_all(df$decimalLatitude, "-", ""), "N"))
    lon <- ifelse(test=df$decimalLongitude < 0, yes=paste(str_replace_all(df$decimalLongitude, "-", ""), "W"), no=paste(str_replace_all(df$decimalLongitude, "-", ""), "E"))
    lat.lon <- paste(lat, lon)
    fas <- ifelse(test=df$basisOfRecord=="MaterialSample", #
        yes=paste0(">", df$otherCatalogNumbers, "_", gene, " ", "[organism=", organism.name, "]", " ", "[Bio_material=", specimen.code, "]", " ", "[Collection_date=", date, "]", " ", "[Country=", country, "]", " ", "[Lat_Lon=", lat.lon, "]",  " ", "[Identified_by=", id.by, "]"), #
        no=paste0(">", df$otherCatalogNumbers, "_", gene, " ", "[organism=", organism.name, "]", " ", "[Specimen-voucher=", specimen.code, "]", " ", "[Collection_date=", date, "]", " ", "[Country=", country, "]", " ", "[Lat_Lon=", lat.lon, "]", " ", "[Identified_by=", id.by, "]"))
    fas <- str_replace_all(string=fas, pattern=" \\[Lat_Lon=NA NA\\]| \\[Country=NA\\]| \\[Collection_date=NA\\]", replacement="")
    fas <- paste(fas, df$rag1, sep="\n")
        return(fas)
}



# feature table containing the locations of attributes of the sequence. CYTB
gb_features_cytb <- function(df, gene, product){
    specimen.code <- ifelse(test=df$basisOfRecord=="MaterialSample", 
        yes=ifelse(test=is.na(df$collectionCode), yes=paste0(df$institutionCode, "-", df$otherCatalogNumbers), no=paste0(df$institutionCode, "-", df$collectionCode, "-", df$otherCatalogNumbers)), #
        no=ifelse(test=is.na(df$collectionCode), yes=paste0(df$institutionCode, "-", df$catalogNumber), no=paste0(df$institutionCode, "-", df$collectionCode, "-", df$catalogNumber)))
    feature.tab <- paste0(paste0(">Feature", " ", df$otherCatalogNumbers, "_", gene),"\n", # 
        "<1", "\t", ">", nchar(df$cytb), "\t", "gene", "\n", #
        "\t", "\t", "\t", "gene", "\t", gene, "\n", #
        "<1", "\t", ">", nchar(df$cytb), "\t", "CDS", "\t", "\t", "\n", #
        "\t", "\t", "\t", "product", "\t", product, "\n", #
        "\t", "\t", "\t", "codon_start", "\t", "1")#
    return(feature.tab)
}

# feature table containing the locations of attributes of the sequence. RAG1
gb_features_rag1 <- function(df, gene, product){
    specimen.code <- ifelse(test=df$basisOfRecord=="MaterialSample", 
        yes=ifelse(test=is.na(df$collectionCode), yes=paste0(df$institutionCode, "-", df$otherCatalogNumbers), no=paste0(df$institutionCode, "-", df$collectionCode, "-", df$otherCatalogNumbers)), #
        no=ifelse(test=is.na(df$collectionCode), yes=paste0(df$institutionCode, "-", df$catalogNumber), no=paste0(df$institutionCode, "-", df$collectionCode, "-", df$catalogNumber)))
    feature.tab <- paste0(paste0(">Feature", " ", df$otherCatalogNumbers, "_", gene),"\n", # 
        "<1", "\t", ">", nchar(df$rag1), "\t", "gene", "\n", #
        "\t", "\t", "\t", "gene", "\t", gene, "\n", #
        "<1", "\t", ">", nchar(df$rag1), "\t", "mRNA", "\n", #
        "\t", "\t", "\t", "gene", "\t", gene, "\n", #
        "\t", "\t", "\t", "product", "\t", product, "\n", #
        "<1", "\t", ">", nchar(df$rag1), "\t", "CDS", "\t", "\t", "\n", #
        "\t", "\t", "\t", "gene", "\t", gene, "\n", #
        "\t", "\t", "\t", "product", "\t", product, "\n", #
        "\t", "\t", "\t", "codon_start", "\t", "1")#
    return(feature.tab)
}
