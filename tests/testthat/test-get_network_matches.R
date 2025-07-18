test_that("get_network_matches returns expected columns and types", {
  genes <- c("GENE1", "GENE2")
  mock_network <- data.frame(
    source = c("GENE1", "GENE2"),
    target = c("GENE2", "GENE3"),
    edge_type = c("Repression", "Activation"),
    evidence = c("PMID:1", "PMID:2"),
    stringsAsFactors = FALSE
  )
  res <- get_network_matches(genes, organism = "human", aggregator = "mean", network_data = mock_network)
  expect_true(is.data.frame(res))
  expect_true(all(c("ENSG", "repression_matches", "activation_matches") %in% colnames(res)))
  expect_true(all(genes %in% res$ENSG))
})

test_that("get_network_matches errors on unsupported organism", {
  expect_error(
    get_network_matches(c("GENE1"), organism = "yeast"),
    "Organism not supported"
  )
})

test_that("message_passing returns graph with correct attributes", {
  # Create a small igraph object
  g <- igraph::graph_from_data_frame(
    data.frame(source = "A", target = "B"),
    directed = TRUE
  )
  igraph::V(g)$feature <- c(1, 0)
  transcription_network <- data.frame(
    source = "A", target = "B", edge_type = "Repression", evidence = "PMID:1"
  )
  g2 <- message_passing(g, transcription_network, agg = "mean")
  expect_true("feature_rep" %in% igraph::vertex_attr_names(g2))
  expect_true("features_act" %in% igraph::vertex_attr_names(g2))
})

test_that("as_dataframe_mp_matches returns correct data frame", {
  g <- igraph::graph_from_data_frame(
    data.frame(source = "A", target = "B"),
    directed = TRUE
  )
  igraph::V(g)$feature_rep <- c(1, 2)
  igraph::V(g)$features_act <- c(3, 4)
  res <- as_dataframe_mp_matches(g)
  expect_true(is.data.frame(res))
  expect_true(all(c("ENSG", "repression_matches", "activation_matches") %in% colnames(res)))
  expect_equal(nrow(res), 2)
})
