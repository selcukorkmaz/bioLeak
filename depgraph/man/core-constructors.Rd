\name{graph_node_set}
\alias{graph_node_set}
\alias{graph_edge_set}
\alias{dependency_graph}
\alias{new_depgraph_nodes}
\alias{new_depgraph_edges}
\alias{new_depgraph}
\alias{graph_query_result}
\alias{dependency_constraint}
\alias{split_constraint}
\alias{leakage_constraint}
\title{Construct Core depgraph S3 Objects}
\usage{
graph_node_set(data = NULL, schema_version = .depgraph_schema_version, source = list())

graph_edge_set(data = NULL, schema_version = .depgraph_schema_version, source = list())

dependency_graph(nodes, edges, graph, metadata = list(), caches = list())

new_depgraph_nodes(data = NULL, schema_version = .depgraph_schema_version, source = list())

new_depgraph_edges(data = NULL, schema_version = .depgraph_schema_version, source = list())

new_depgraph(nodes, edges, graph = NULL, metadata = list(), caches = list())

graph_query_result(query = "", params = list(), nodes = NULL, edges = NULL, table = NULL, metadata = list())

dependency_constraint(constraint_id, relation_types, sample_map, transitive = TRUE, metadata = list())

split_constraint(strategy, sample_map, recommended_bioleak_args = list(), metadata = list())

leakage_constraint(issue_type, severity, affected_samples, evidence = NULL, recommendation = "", metadata = list())
}
\arguments{
  \item{data}{A data frame matching the canonical schema for nodes or edges.}
  \item{schema_version}{Schema version string stored on the object.}
  \item{source}{Optional source metadata.}
  \item{nodes, edges}{A \code{graph_node_set} and \code{graph_edge_set}.}
  \item{graph}{An internal \code{igraph} object.}
  \item{metadata, caches, params, recommended_bioleak_args}{Named lists with auxiliary metadata.}
  \item{query}{Query label stored on a \code{graph_query_result}.}
  \item{table}{Tabular query result payload.}
  \item{constraint_id, relation_types, transitive}{Fields describing a dependency constraint.}
  \item{sample_map}{Sample-level mapping table for constraints.}
  \item{strategy}{Split strategy identifier.}
  \item{issue_type, severity, affected_samples, evidence, recommendation}{Fields describing a leakage warning.}
}
\value{
An S3 object corresponding to the constructor that was called.
}
\description{
Low-level constructors for the core S3 classes used throughout \pkg{depgraph}.
}
