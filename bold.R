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
library(stringr)
library(Biostrings)

##############################################################################
####ESCOGER CARPETA DE TRABAJO################################################
setwd("C:/Users/GERARDO/Desktop/bold")


##############################################################################
#IMPORTAR DATOS BOLD EN FORMATO TSV###########################################
#VARIOS GRUPOS
datosc<- bold_seqspec(taxon = c("Caprella scaura",
                               "Caprella scauroides",
                               "Caprella californica"), format = "tsv")




###############################################################################
#ESCOGER LA DATA DE MI INTERES#################################################
caprella_1<- caprella %>%
  filter(markercode == "COI-5P") %>% #INDICAR EL MARCADOR DE INTRES
  mutate(valid_nucleotide_count = str_count(nucleotides, "[^-]")) %>%  # Contar solo caracteres distintos de "-"
  filter(valid_nucleotide_count > 599) %>%  # Filtrar por cantidad de caracteres sin contar guiones
  select(processid,
         genbank_accession, 
         genus_name, 
         species_name, 
         country, 
         province_state, 
         trace_names,
         trace_links,
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
caprella_1$genbank_id <- sapply(caprella_1$genbank_accession, get_genbank_id)





################################################################################
#OBTENER DATA DEL FORMATO GENBANK DE LAS SECUENCIAS DE MI INTERES###############
# Filtrar solo IDs válidos (sin NA ni valores vacíos)
valid_ids <- na.omit(caprella_1$genbank_id)
valid_ids <- valid_ids[valid_ids != ""]  # Quitar valores vacíos


# Inicializar las nuevas columnas en datos_coi5p con NA
caprella_1$Title <- NA
caprella_1$Authors <- NA
caprella_1$Journal <- NA

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
    caprella_1$Title[caprella_1$genbank_id == id] <- title
    caprella_1$Authors[caprella_1$genbank_id == id] <- authors
    caprella_1$Journal[caprella_1$genbank_id == id] <- journal
  }
}


##########################################################################
####FILTRAR SECUENCIAS ASOCIADAS A ARTICULOS CIENTIFICOS##################
caprella_2 <- caprella_1 %>%
  filter(!(trace_names == "" & Journal == "Unpublished"))




##########################################################################
#########################LLENAR ESPACION VACIOS EN ESPECIE################
datos_coi5p_1 <- datos_coi5p_1 %>%
  mutate(species_name = if_else(species_name == "" | is.na(species_name), 
                                genus_name, species_name))


##SIJORGECREEQUEDEBERIAGREGARMASGRUPOSTAXONOMICO-DEBERIAEMPEZARAQUI
#########################################################################
#########################################################################
#########################################################################
#########################################################################
###DESCARGAR LOS ELECTROFEROGRAMAS#######################################
#########################################################################
#########################################################################
#DESCARGAR TODOS MIS ELECTROFEROGRAMAS###################################
#########################################################################
for (i in 1:nrow(datos_coi5p_1)) {
  if (!is.na(datos_coi5p_1$trace_links[i]) && datos_coi5p_1$trace_links[i] != "") {
    links <- strsplit(datos_coi5p_1$trace_links[i], "\\|")[[1]]  # Separar links
    processid <- datos_coi5p_1$processid[i]  # Obtener el processid
    
    for (j in seq_along(links)) {
      url <- links[j]
      file_name <- paste0("bold_traces/", processid, "_", j, ".ab1")  # Nombre del archivo
      
      download.file(url, file_name, mode = "wb")
      cat("Descargado:", file_name, "\n")
    }
  }
}



# Definir el directorio con los archivos .ab1
setwd("C:/Users/GERARDO/Desktop/bold/bold_traces")


# Listar todos los archivos .ab1 en la carpeta
archivos_ab1 <- list.files(pattern = "\\.ab1$", full.names = TRUE)


# Función para calcular la calidad de cada archivo
calcular_calidad <- function(archivo) {
  tryCatch({
    trace <- read.abif(archivo)  # Leer el archivo .ab1
    mean_quality <- mean(trace@data$PCON.1, na.rm = TRUE)  # Calcular calidad promedio
    return(mean_quality)
  }, error = function(e) {
    cat("Error con:", archivo, "\n")  # Si hay un error, mostrar mensaje
    return(NA)
  })
}


# Aplicar la función a todos los archivos y guardar en un data.frame
resultados_calidad <- data.frame(
  archivo = archivos_ab1,
  mean_quality = sapply(archivos_ab1, calcular_calidad)
)

print(resultados_calidad)




# Filtrar solo los archivos que terminan en "_1.ab1" o "_2.ab1"
resultados_filtrados <- resultados_calidad %>%
  filter(grepl("_1\\.ab1$|_2\\.ab1$", archivo)) %>%
  mutate(processid = gsub("_\\d+\\.ab1", "", archivo))  # Extraer el identificador sin _1 o _2


# Reestructurar los datos para separar calidadF (Forward) y calidadR (Reverse)
resultados_wide <- resultados_filtrados %>%
  mutate(tipo = ifelse(grepl("_1\\.ab1$", archivo), "calidadF", "calidadR")) %>%
  select(processid, tipo, mean_quality) %>%
  pivot_wider(names_from = tipo, values_from = mean_quality)

resultados_wide <- resultados_wide %>%
  mutate(processid = gsub("^\\./", "", processid))  # Eliminar el prefijo "./"


# Unir los datos con `datos_coi5p_fil_3` según `processid`
datos_coi5p_2 <- datos_coi5p_1 %>%
  left_join(resultados_wide, by = "processid")

#usar el critero de filtro
datos_coi5p_3 <- datos_coi5p_2%>%
  filter(is.na(calidadF) | calidadF > 30)
View(datos_coi5p_2)

# Guardar el data.frame en un archivo CSV (compatible con Excel)
setwd("C:/Users/GERARDO/Desktop/bold")
write.csv(datos_coi5p_3, "datos_coi5p_3.csv", row.names = FALSE)




#########################################################################
##################extraer los fastas#####################################
# Generar líneas en formato FASTA ordenadas por species_name
fasta_lines <- datos_coi5p_3 %>%
  arrange(species_name) %>%  # Ordenar alfabéticamente por species_name
  mutate(fasta_header = paste0(">", genus_name, "_", species_name, "_", country, "_", processid)) %>%
  select(fasta_header, nucleotides) %>%
  apply(1, paste, collapse = "\n") %>%
  paste(collapse = "\n\n")  # Agregar un espacio entre secuencias


#Aqui se guardan todas las secuencias de mi interes. 
writeLines(fasta_lines, "secuencias.fasta")

# Generar líneas en formato FASTA ordenadas por cirripedos
fasta_lines_cirripedos <- datos_coi5p_3 %>%
  filter(genus_name == "Membranobalanus" | genus_name == "Conopea" | 
           genus_name == "Striatobalanus" |
           genus_name == "Solidobalanus"  |
           genus_name == "Chirona")%>%
  arrange(species_name) %>%  # Ordenar alfabéticamente por species_name
  mutate(fasta_header = paste0(">", genus_name, "_", species_name, "_", country, "_", processid)) %>%
  select(fasta_header, nucleotides) %>%
  apply(1, paste, collapse = "\n") %>%
  paste(collapse = "\n\n")  # Agregar un espacio entre secuencias
writeLines(fasta_lines_cirripedos, "secuencias_cirripedos.fasta")




#Generar líneas en formato FASTA ordenadas por microphrys
fasta_lines_microphrys <- datos_coi5p_3 %>%
  filter(genus_name == "Microphrys")%>%
  arrange(species_name) %>%  # Ordenar alfabéticamente por species_name
  mutate(fasta_header = paste0(">", genus_name, "_", species_name, "_", country, "_", processid)) %>%
  select(fasta_header, nucleotides) %>%
  apply(1, paste, collapse = "\n") %>%
  paste(collapse = "\n\n")  # Agregar un espacio entre secuencias
writeLines(fasta_lines_microphrys, "secuencias_microphrys.fasta")





#Geenerar líneas en formato FASTA ordenadas por stenothoe
fasta_lines_stenothoe <- datos_coi5p_3 %>%
  filter(genus_name == "Stenothoe")%>%
  arrange(species_name) %>%  # Ordenar alfabéticamente por species_name
  mutate(fasta_header = paste0(">", genus_name, "_", species_name, "_", country, "_", processid)) %>%
  select(fasta_header, nucleotides) %>%
  apply(1, paste, collapse = "\n") %>%
  paste(collapse = "\n\n")  # Agregar un espacio entre secuencias
writeLines(fasta_lines_stenothoe, "secuencias_stenothoe.fasta")



#Geenerar líneas en formato FASTA ordenadas por stenothoe
fasta_lines_pilumnus <- datos_coi5p_3 %>%
  filter(genus_name == "Pilumnus")%>%
  arrange(species_name) %>%  # Ordenar alfabéticamente por species_name
  mutate(fasta_header = paste0(">", genus_name, "_", species_name, "_", country, "_", processid)) %>%
  select(fasta_header, nucleotides) %>%
  apply(1, paste, collapse = "\n") %>%
  paste(collapse = "\n\n")  # Agregar un espacio entre secuencias
writeLines(fasta_lines_pilumnus, "fasta_lines_pilumnus.fasta")



#Geenerar líneas en formato FASTA ordenadas por stenothoe
fasta_lines_caprella <- caprella_1 %>%
  filter(genus_name == "Caprella")%>%
  arrange(species_name) %>%  # Ordenar alfabéticamente por species_name
  mutate(fasta_header = paste0(">", genus_name, "_", species_name, "_", country, "_", processid)) %>%
  select(fasta_header, nucleotides) %>%
  apply(1, paste, collapse = "\n") %>%
  paste(collapse = "\n\n")  # Agregar un espacio entre secuencias
writeLines(fasta_lines_caprella, "fasta_lines_caprella.fasta")

setwd("C:/Users/GERARDO/Desktop/bold")












