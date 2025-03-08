################INSTALAR Y CARGAR PAQUETES NECESARIOS#########################
install.packages("dplyr")
install.packages("tidyverse")
install.packages("tidyr")
install.packages("rentrez")
install.packages("purrr")
install.packages("parallel")
install.packages("foreach")
install.packages("doParallel")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("sangerseqR")
remotes::install_github("ropensci/bold")
library(bold)
library(rentrez)
library(purrr)
library(parallel)
library(foreach)
library(doParallel)
library(dplyr)
library(tidyr)
library(tidyverse)

##############################################################################
####ESCOGER CARPETA DE TRABAJO################################################
setwd("C:/Users/GERARDO/Desktop/bold")


##############################################################################
#IMPORTAR DATOS BOLD EN FORMATO TSV###########################################
#01 SOLO GRUPO
datos_pilum <- bold_seqspec(taxon = "Pilumnus", format = "tsv")
datos_micro <- bold_seqspec(taxon = "Microphrys", format = "tsv")
datos_acan <- bold_seqspec(taxon = "Acanthonyx petiverii", format = "tsv")
datos_conopea <- bold_seqspec(taxon = "Conopea", format = "tsv")
datos_bal <- bold_seqspec(taxon = "Balanus", format = "tsv")
datos_ste <- bold_seqspec(taxon = "Stenothoe", format = "tsv")
datos_syna <- bold_seqspec(taxon = "Synalpheus", format = "tsv")

#VARIOS GRUPOS
datos<- bold_seqspec(taxon = c("Pilumnus", "Microphrys", "Acanthonyx petiverii", 
                               "Conopea", "Balanus", "Stenothoe", 
                               "Synalpheus"), format = "tsv")


###############################################################################
#ESCOGER LA DATA DE MI INTERES#################################################
datos_coi5p <- datos %>%
  filter(markercode == "COI-5P") %>% #INDICAR EL MARCADOR DE INTRES
  mutate(valid_nucleotide_count = str_count(nucleotides, "[^-]")) %>%  # Contar solo caracteres distintos de "-"
  filter(valid_nucleotide_count > 500) %>%  # Filtrar por cantidad de caracteres sin contar guiones
  select(processid,
         genbank_accession, 
         genus_name, 
         species_name, 
         country, 
         province_state, 
         trace_names,
         nucleotides)


##################################################################################
#IMPORTAR EL GenBank ID DEL NCBI USANDO EL "genbank_accession" ###########################
# Función para obtener el ID de GenBank desde el número de acceso
get_genbank_id <- function(accession) {
  if (is.na(accession) | accession == "") {
    return(NA)
  }
  
  search_result <- tryCatch(
    entrez_search(db = "nucleotide", term = accession),
    error = function(e) return(NA)
  )
  
  if (length(search_result$ids) > 0) {
    Sys.sleep(0.5)  # Evita bloqueos en la API
    return(search_result$ids[1])
  } else {
    return(NA)
  }
}

#Obtener GenBank ID para cada número de acceso en la columna "genbank_accession"
datos_coi5p$genbank_id <- sapply(datos_coi5p$genbank_accession, get_genbank_id)





################################################################################
#OBTENER DATA DEL FORMATO GENBANK DE LAS SECUENCIAS DE MI INTERES###############
# Filtrar solo IDs válidos (sin NA ni valores vacíos)
valid_ids <- na.omit(datos_coi5p$genbank_id)
valid_ids <- valid_ids[valid_ids != ""]  # Quitar valores vacíos


# Inicializar las nuevas columnas en datos_coi5p con NA
datos_coi5p$Title <- NA
datos_coi5p$Authors <- NA
datos_coi5p$Journal <- NA

# Recorrer solo los IDs válidos
for (i in seq_along(valid_ids)) {
  
  # Obtener ID GenBank
  id <- valid_ids[i]
  
  # Validar que el ID sea un número válido antes de hacer la consulta
  if (!is.na(id) & id != "") {
    
    # Intentar obtener datos de GenBank
    df <- tryCatch(
      entrez_fetch(db="nucleotide", id=id, rettype="gb", retmode="text"),
      error = function(e) return(NA)
    )
    
    # Si la consulta falla, continuar con el siguiente ID
    if (is.na(df)) next
    
    # Dividir el texto en líneas
    lines <- unlist(strsplit(df, "\n"))
    
    # Extraer el primer TITLE
    title_line <- lines[grep("^  TITLE", lines)]
    title <- ifelse(length(title_line) > 0, sub("^  TITLE +", "", title_line[1]), NA)
    
    # Extraer el primer AUTHORS
    authors_line <- lines[grep("^  AUTHORS", lines)]
    authors <- ifelse(length(authors_line) > 0, sub("^  AUTHORS +", "", authors_line[1]), NA)
    
    # Extraer el primer JOURNAL
    journal_line <- lines[grep("^  JOURNAL", lines)]
    journal <- ifelse(length(journal_line) > 0, sub("^  JOURNAL +", "", journal_line[1]), NA)
    
    
    # Guardar resultados en las columnas de datos_coi5p (asignar en la fila correcta)
    datos_coi5p$Title[datos_coi5p$genbank_id == id] <- title
    datos_coi5p$Authors[datos_coi5p$genbank_id == id] <- authors
    datos_coi5p$Journal[datos_coi5p$genbank_id == id] <- journal
  }
}


##########################################################################
####FILTRAR SECUENCIAS ASOCIADAS A ARTICULOS CIENTIFICOS##################
datos_coi5p_fil <- datos_coi5p %>%
  filter(!(trace_names == "" & Journal == "Unpublished"))


##########################################################################
###########EXPORTAR EN FORMATO EXCEL######################################
write.csv(datos_coi5p_fil, "datos_coi5p_fil.csv", row.names = FALSE) 



