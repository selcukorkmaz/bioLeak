# S4 classes for results and audits

setClass("LeakSplits",
         slots = c(mode = "character",
                   indices = "list",    # list of folds: each has train, test integer vectors
                   info = "list"))      # metadata (group/batch/study/time columns, seeds)

setClass("LeakFit",
         slots = c(
           splits = "LeakSplits",
           metrics = "data.frame",
           metric_summary = "data.frame",  # <── yeni
           audit = "data.frame",           # <── yeni
           predictions = "list",
           preprocess = "list",
           learners = "list",
           outcome = "character",
           task = "character",
           feature_names = "character",
           info = "list"
         ))
# seeds, hashes

setClass("LeakAudit",
         slots = c(
           fit = "LeakFit",
           permutation_gap = "data.frame",
           perm_distribution = "numeric",   # <── eklendi
           batch_assoc = "data.frame",
           duplicates = "data.frame",
           trail = "list",
           info = "list"))
