\name{query_neighbors}
\alias{query_neighbors}
\alias{query_node_type}
\alias{query_edge_type}
\alias{query_paths}
\alias{query_shortest_paths}
\alias{detect_dependency_components}
\alias{detect_shared_dependencies}
\title{Query Dependency Graph Structure}
\usage{
query_neighbors(graph, node_ids, edge_types = NULL, node_types = NULL,
  direction = c("out", "in", "all"))

query_node_type(graph, node_types, ids = NULL)

query_edge_type(graph, edge_types, node_ids = NULL)

query_paths(graph, from, to, edge_types = NULL, node_types = NULL,
  mode = c("out", "in", "all"), max_length = NULL)

query_shortest_paths(graph, from, to, edge_types = NULL, node_types = NULL,
  mode = c("out", "in", "all"))

detect_dependency_components(graph,
  via = c("Subject", "Batch", "Study", "Timepoint", "Assay", "FeatureSet", "Outcome"),
  edge_types = NULL, min_size = 1)

detect_shared_dependencies(graph, via = c("Subject", "Batch", "Study",
  "Timepoint"), samples = NULL)
}
\arguments{
  \item{graph}{A \code{dependency_graph}.}
  \item{node_ids, from, to}{Node identifiers to use as query seeds or endpoints.}
  \item{edge_types}{Optional edge types used to filter the traversal graph or edge table.}
  \item{node_types}{Optional node types used to filter node results or allowed path members.}
  \item{direction, mode}{Traversal direction.}
  \item{ids}{Optional node identifiers used to further restrict \code{query_node_type()}.}
  \item{max_length}{Maximum path length for \code{query_paths()}.}
  \item{via}{Dependency node types used for sample-level dependency detection.}
  \item{min_size}{Minimum component size retained by \code{detect_dependency_components()}.}
  \item{samples}{Optional sample identifiers or sample node IDs used to restrict direct shared-dependency detection.}
}
\value{
Each function returns a \code{graph_query_result}. Use
\code{as.data.frame()} to obtain the tidy result table.
}
\description{
Query graph neighborhoods, typed nodes and edges, path structure, projected
sample dependency components, and direct shared dependencies within a
\code{dependency_graph}.
}
