# S4 classes for results and audits

setClass("LeakSplits",
         slots = c(mode = "character",
                   indices = "list",    # list of folds: each has train, test integer vectors
                   info = "list"))      # metadata (group/batch/study/time columns, seeds)

setClass("LeakFit",
         slots = c(splits = "LeakSplits",
                   metrics = "data.frame",   # per-fold metrics
                   predictions = "list",     # per-fold data.frames: id, truth, pred, (surv if any)
                   preprocess = "list",      # per-fold fitted preprocessors
                   learners = "list",        # per-fold fitted learners
                   outcome = "character",    # outcome column name
                   task = "character",       # "binomial","gaussian","survival"
                   feature_names = "character",
                   info = "list"))           # seeds, hashes

setClass("LeakAudit",
         slots = c(fit = "LeakFit",
                   permutation_gap = "data.frame",
                   batch_assoc = "data.frame",
                   duplicates = "data.frame",
                   trail = "list",
                   info = "list"))
