# Expect correct filter for default settings.

    Code
      suppressMessages(SeedMatchR::build_annotation_filter(standard.chroms = TRUE,
        remove.na.rows = TRUE, protein.coding = TRUE))
    Output
      AnnotationFilterList of length 3 
      ((seq_name %in% c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y')) & tx_id != '') & gene_biotype == 'protein_coding' & tx_biotype == 'protein_coding'

# Expect correct filter for selecting symbol and entrez ID with TxBiotypeFilter enabled.

    Code
      suppressMessages(SeedMatchR::build_annotation_filter(standard.chroms = TRUE,
        remove.na.rows = TRUE, protein.coding = TRUE, canonical = FALSE, symbol = c(
          "Ttr", "Hao1"), entrez_id = c("1", "300"), add.filter = AnnotationFilter::TxBiotypeFilter(
          "Nonsense Mediated Decay")))
    Output
      AnnotationFilterList of length 2 
      (((((seq_name %in% c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y')) & tx_id != '') & gene_biotype == 'protein_coding' & tx_biotype == 'protein_coding') & symbol %in% c('Ttr', 'Hao1')) & entrez %in% c('1', '300')) & tx_biotype == 'Nonsense Mediated Decay'

# Expect correct type for added AnnotationFilter.

    Code
      suppressMessages(SeedMatchR::build_annotation_filter(standard.chroms = TRUE,
        remove.na.rows = TRUE, protein.coding = TRUE, canonical = TRUE, symbol = c(
          "Ttr", "Hao1"), entrez_id = c("1", "300"), add.filter = AnnotationFilter::TxBiotypeFilter(
          "Nonsense Mediated Decay")))
    Output
      AnnotationFilterList of length 2 
      ((((((seq_name %in% c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y')) & tx_id != '') & gene_biotype == 'protein_coding' & tx_biotype == 'protein_coding') & tx_is_canonical == '1') & symbol %in% c('Ttr', 'Hao1')) & entrez %in% c('1', '300')) & tx_biotype == 'Nonsense Mediated Decay'

