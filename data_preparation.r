# -----------------------------------------------------------------------------
# Data Preparation for sRNA Homolog Heatmap
#
# Description:
# This script processes multiple FASTA files containing sRNA homolog data from 
# GLASSgo to prepare data for a Shiny-based heatmap visualization. It extracts 
# metadata from FASTA headers, retrieves taxonomic information using taxonkit,
# and compiles the data into an RData file for the Shiny app.
#

# --- 1. Define Helper Functions ---

#' Parse GLASSgo FASTA Headers
#'
#' Extracts relevant information from a single GLASSgo FASTA header.
#'
#' @param header A character string representing a FASTA header.
#' @return A character vector with accession number, strand, start, end,
#'         strain name, unique ID, and sequence identity.
parse_glassgo_header <- function(header) {
  # Remove the leading ">"
  header <- sub(">", "", header)
  
  # Split header to extract identity and main content
  parts <- strsplit(header, "p.c.VAL:")[[1]]
  identity <- gsub("%.*", "", parts[2])
  
  main_content <- parts[1]
  
  # Extract accession, coordinates, and name
  content_parts <- strsplit(main_content, ":")[[1]]
  accession <- content_parts[1]
  
  name_section <- strsplit(content_parts[2], " ")[[1]]
  full_name <- paste(name_section[2:(length(name_section) - 1)], collapse = " ")
  strain_name <- strsplit(full_name, ",")[[1]][1]
  
  # Handle coordinate extraction and strand determination
  coords <- strsplit(name_section[1], "-")[[1]]
  coord1 <- as.numeric(gsub("c", "", coords[1]))
  coord2 <- as.numeric(coords[2])
  
  if (coord2 < coord1) {
    start_pos <- coord2
    end_pos <- coord1
    strand <- "-"
  } else {
    start_pos <- coord1
    end_pos <- coord2
    strand <- "+"
  }
  
  unique_id <- paste(accession, start_pos, sep = "_")
  
  return(c(accession, strand, start_pos, end_pos, strain_name, unique_id, identity))
}

#' Process a FASTA File
#'
#' Reads a FASTA file, extracts headers and sequences, and processes the headers.
#'
#' @param file_content A character vector representing the lines of a FASTA file.
#' @return A data frame with parsed information for each sequence.
process_fasta_file <- function(file_content) {
  header_lines <- grep(">", file_content)
  headers <- as.character(file_content[header_lines])
  sequences <- as.character(file_content[header_lines + 1])
  
  # Handle potential empty first lines
  if (length(grep(":", headers[1])) == 0) {
    headers <- headers[-1]
    sequences <- sequences[-1]
  }
  
  parsed_data <- do.call(rbind, lapply(headers, parse_glassgo_header))
  
  # Check if parsed_data has the expected number of columns
  if (ncol(parsed_data) == 7) {
    colnames(parsed_data) <- c("Accession_number", "Strand", "start", "end", "name", "ID", "identity")
    return(data.frame(parsed_data, Full_header = headers, sequence = sequences, stringsAsFactors = FALSE))
  } else {
    warning("Header parsing resulted in an incorrect number of columns. Returning NULL.")
    return(NULL)
  }
}


# --- 2. Fetch Taxonomic Data ---

# Get a list of FASTA files in the working directory
files<-dir()
files_to_process<-files[grepl("\\.fasta$|\\.fa$", files)]

taxid_lookup <- new.env(hash = TRUE)
processed_taxids <- c()

for (fasta_path in files_to_process) {
  fasta_file_path <- fasta_path
  
  if (!file.exists(fasta_file_path)) {
    next
  }
  
  fasta_content <- readLines(fasta_file_path)
  
  if (length(fasta_content) < 4) {
    next
  }
  
  coordinates_data <- process_fasta_file(fasta_content[3:length(fasta_content)])
  
  # Extract unique taxIDs
  taxids_from_fasta <- gsub(".*taxID:", "", coordinates_data$Full_header)
  taxids_from_fasta <- gsub(" synten.*|\\;.*", "", taxids_from_fasta)
  unique_taxids <- setdiff(unique(taxids_from_fasta), processed_taxids)
  
  if (length(unique_taxids) == 0) {
    next
  }
  
  # Write taxIDs to a temporary file for taxonkit
  temp_taxid_file <- tempfile()
  writeLines(unique_taxids, con = temp_taxid_file)
  
  # Use taxonkit to get lineage information
  taxonkit_output <- system(paste0("taxonkit lineage --show-lineage-ranks -t ", temp_taxid_file), intern = TRUE)
  
  # Process taxonkit output
  for (line in taxonkit_output) {

    parts <- strsplit(line, "\t")[[1]]
    taxid <- parts[1]
    
    # Get all accessions for the current taxid
    accessions <- coordinates_data$Accession_number[grep(paste0("taxID:", taxid, ""), coordinates_data$Full_header)]
    accessions <- gsub("\\..*", "", accessions)
    
    for (acc in accessions) {
      if (is.null(taxid_lookup[[acc]])) {
        lineage <- strsplit(parts[2], ";")[[1]]
        ranks <- strsplit(parts[4], ";")[[1]]
        names(lineage) <- ranks
        taxid_lookup[[acc]] <- lineage
      }
    }
  }
  
  processed_taxids <- unique(c(processed_taxids, unique_taxids))
}

# Save the lookup table
save(taxid_lookup, file = "taxid_lookup.Rdata")

# --- 3. Consolidate Taxonomic Information ---

load("taxid_lookup.Rdata")

ranks <- c("kingdom", "phylum", "class", "order", "family", "genus", "species", "subspecies", "serotype", "strain")
accession_names <- names(taxid_lookup)

taxonomy_matrix <- matrix(NA, length(accession_names), length(ranks))
colnames(taxonomy_matrix) <- ranks
rownames(taxonomy_matrix) <- accession_names

for (acc in accession_names) {
  tax_data <- taxid_lookup[[acc]]
  
  # Handle "no rank" and serovar special cases
  no_rank_indices <- grep("no rank", names(tax_data))
  if (length(no_rank_indices) > 0) {
    serovar_indices <- grep("serovar", tax_data[no_rank_indices])
    if (length(serovar_indices) > 0) {
      names(tax_data)[no_rank_indices[serovar_indices]] <- "serotype"
    }
  }
  
  matched_indices <- match(ranks, names(tax_data))
  taxonomy_matrix[acc, ] <- tax_data[matched_indices]
}

# --- 4. Combine Sequence Data with Taxonomy ---

all_coordinates_list <- list()

for (fasta_path in files_to_process) {
  fasta_file_path <- fasta_path
  
  if (!file.exists(fasta_file_path)) {
    next
  }
  
  fasta_content <- readLines(fasta_file_path)
  
  if (length(fasta_content) < 4) {
    next
  }
  
  coordinates_data <- process_fasta_file(fasta_content[3:length(fasta_content)])
  
  accessions <- gsub("\\..*", "", coordinates_data$Accession_number)
  taxonomy_matches <- match(accessions, rownames(taxonomy_matrix))
  
  combined_data <- cbind(coordinates_data[,c("Accession_number","identity")], taxonomy_matrix[taxonomy_matches, ])
  colnames(combined_data)[c(1,2)]<-c("Accession_number","identity")
  dir_name <- gsub(".fasta|fa", "", fasta_file_path)
  all_coordinates_list[[dir_name]] <- combined_data
}

# --- 5. Prepare and Save Final Data for Shiny App ---

# Extract a consolidated taxonomy table from all entries
tax_table_list <- list()
for (i in 1:length(all_coordinates_list)) {
  for (j in 1:nrow(all_coordinates_list[[i]])) {
    accession <- all_coordinates_list[[i]][j, "Accession_number"]
    if (is.null(tax_table_list[[accession]])) {
      tax_table_list[[accession]] <- all_coordinates_list[[i]][j, 2:11]
    }
  }
}
consolidated_tax_table <- do.call("rbind", tax_table_list)

# reduce data size
for (i in 1:length(all_coordinates_list)) {
	all_coordinates_list[[i]]<-all_coordinates_list[[i]][,1:2]
}

# Save the final data object for the Shiny app
heatdata <- list(consolidated_tax_table, all_coordinates_list)
save(heatdata, file = "heatdata.Rdata")

