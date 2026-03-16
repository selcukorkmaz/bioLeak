\name{build_dependency_graph}
\alias{build_dependency_graph}
\alias{build_depgraph}
\alias{validate_graph}
\alias{validate_depgraph}
\alias{as_igraph}
\title{Assemble and Validate Dependency Graphs}
\usage{
build_dependency_graph(nodes, edges, graph_name = NULL, dataset_name = NULL, validate = TRUE,
  validation_overrides = list())

build_depgraph(nodes, edges, graph_name = NULL, dataset_name = NULL, validate = TRUE,
  validation_overrides = list())

validate_graph(graph, checks = c("ids", "references", "cardinality", "schema", "time"),
  error_on_fail = FALSE, levels = NULL, severities = NULL)

validate_depgraph(graph, checks = c("ids", "references", "cardinality", "schema", "time"),
  error_on_fail = FALSE, levels = NULL, severities = NULL)

as_igraph(x)
}
\arguments{
  \item{nodes, edges}{Lists of \code{graph_node_set} and \code{graph_edge_set} objects.}
  \item{graph_name, dataset_name}{Optional metadata labels.}
  \item{validate}{If \code{TRUE}, run \code{validate_graph()} before returning.}
  \item{validation_overrides}{Optional named list of explicit validation exceptions.}
  \item{graph}{A \code{dependency_graph}.}
  \item{checks}{Validation checks to run.}
  \item{error_on_fail}{If \code{TRUE}, stop when validation errors are found.}
  \item{levels}{Optional validation layers to run.}
  \item{severities}{Optional severities to retain in the returned report.}
  \item{x}{A \code{dependency_graph}.}
}
\value{
For \code{build_dependency_graph()}, a \code{dependency_graph}. For
\code{validate_graph()} and \code{validate_depgraph()}, a
\code{depgraph_validation_report}. For \code{as_igraph()}, the underlying
\code{igraph} object.
}
\description{
Combine canonical node and edge tables into a typed dependency graph and
perform structural, semantic, and graph-local leakage-aware validation.
}
\examples{
meta <- data.frame(
  sample_id = c("S1", "S2"),
  subject_id = c("P1", "P2")
)

samples <- create_nodes(meta, type = "Sample", id_col = "sample_id")
subjects <- create_nodes(meta, type = "Subject", id_col = "subject_id")
edges <- create_edges(
  meta,
  "sample_id",
  "subject_id",
  "Sample",
  "Subject",
  "sample_belongs_to_subject"
)

g <- build_dependency_graph(list(samples, subjects), list(edges))
validate_graph(g)
}
