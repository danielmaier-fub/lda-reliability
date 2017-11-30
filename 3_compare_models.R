
require(lda)
require(topicmodels)
require(htmltools)
require(xtable)
require(igraph)
require(cluster)
require(Matrix)
require(slam)
options(stringsAsFactors = F)

#####################################################
## functions to measure significant co-occurrences ##
#####################################################

computeSignificancePoisson = function(ttmatrix, wordsInDocs, pLevel = 0.05) {
  
  if (FALSE) {
    a <- 10
    b <- 200
    wordsInDocs <- round(runif(1000, 10, 100))
    Lc <- sum(wordsInDocs) ^ 2
    Ld <- sum(wordsInDocs ^ 2)
    lambda <- (a * b / Lc) * Ld
    plot(dpois(0:9, lambda), type="o")
    abline(h = 0.05)
    z <- qpois(pLevel, lambda, lower.tail = F)
    sum(dpois(z:100000, lambda))
    z
  }
  
  Lc <- sum(wordsInDocs) ^ 2
  Ld <- sum(wordsInDocs ^ 2)
  
  m <- nrow(ttmatrix)
  n <- ncol(ttmatrix)
  
  sigMatrix <- matrix(0,nrow=m,ncol=n)
  rownames(sigMatrix) <- rownames(ttmatrix)
  colnames(sigMatrix) <- colnames(ttmatrix)
  print("Computing poisson weights for co-occurrences")
  pb <- txtProgressBar(max=m)
  for (i_a in 1:m) {
    setTxtProgressBar(pb, value=i_a)
    for (i_b in 1:n) {
      a = ttmatrix[i_a, i_a]
      
      b = ttmatrix[i_b, i_b]
      
      k = ttmatrix[i_a, i_b]
      
      if (a == 0 || b == 0 || k == 0 || i_a == i_b) {
        sigMatrix[i_a, i_b] = 0;
      } else {
        lambda <- (a * b / Lc) * Ld
        sigMatrix[i_a, i_b] = k - qpois(pLevel, lambda, lower.tail = F)
      }
    }
  }
  close(pb)
  sigMatrix[sigMatrix < 0] <- 0
  return(sigMatrix)
}


computeSignificanceLL = function(ttmatrix, sizeOfCorpus) {
  
  c <- sizeOfCorpus
  
  m <- nrow(ttmatrix)
  n <- ncol(ttmatrix)
  
  sigMatrix <- matrix(0,nrow=m,ncol=n)
  rownames(sigMatrix) <- rownames(ttmatrix)
  colnames(sigMatrix) <- colnames(ttmatrix)
  print("Computing log-likelihood statistic for co-occurrences")
  pb <- txtProgressBar(max=m)
  for (i_a in seq(1,m)) {
    setTxtProgressBar(pb, value=i_a)
    for (i_b in seq(1,n)) {
      a = ttmatrix[i_a, i_a]
      
      b = ttmatrix[i_b,i_b]
      
      k = ttmatrix[i_a,i_b]
      
      if (a == 0 || b == 0 || k == 0 || i_a == i_b) {
        sigMatrix[i_a,i_b] = 0;
      } else {
        if (k==a) a = a + 0.001;
        if (k==b) b = b + 0.001;
        sigMatrix[i_a,i_b] = 2 * ((c * log(c)) - (a * log(a)) - (b * log(b)) + (k * log(k)) # fehlt hier ein MINUSZEICHEN?
                                  + (c - a - b + k) * log(c - a - b + k) 
                                  + (a - k) * log(a - k) + (b - k) * log(b - k) 
                                  - (c - a) * log(c - a) - (c - b) * log(c - b))
      }      
    }
  }
  close(pb)
  return(sigMatrix)
}


computeSignificanceDice = function(ttmatrix, listOfTerms = NULL) {
  
  m <- nrow(ttmatrix)
  n <- ncol(ttmatrix)
  
  if (!is.null(listOfTerms)) {
    ttmatrix <- ttmatrix[listOfTerms,]
    m <- length(listOfTerms)
  }
  
  sigMatrix <- matrix(0,nrow=m,ncol=n)
  rownames(sigMatrix) <- rownames(ttmatrix)
  colnames(sigMatrix) <- colnames(ttmatrix)
  print("Computing dice statistic for co-occurrences")
  pb <- txtProgressBar(max=m)
  for (i_a in seq(1,m)) {
    setTxtProgressBar(pb, value=i_a)
    for (i_b in seq(1,n)) {
      a = ttmatrix[i_a, i_a]
      
      b = ttmatrix[i_b,i_b]
      
      k = ttmatrix[i_a,i_b]
      
      if ((a + b) == 0) {
        sigMatrix[i_a,i_b] = 0;
      } else {
        if (i_a == i_b) {
          sigMatrix[i_a,i_b] = 0;
        } else {
          sigMatrix[i_a,i_b] = (2 * k) / (a + b);
        }
      }
    }
  }
  close(pb)
  return(sigMatrix)
}


getTopCoocs = function(coocSentenceSignificance, decr = TRUE, dtm2 = NULL, N = 300) {
  dictTermsDTM <- colnames(coocSentenceSignificance)
  topCoocsIdx <- arrayInd(order(coocSentenceSignificance, decreasing=decr)[1:(N*2)], dim(coocSentenceSignificance))
  
  if (is.null(dtm2)) {
    topCoocs <- t(apply(topCoocsIdx, 1, FUN=function(x) c(dictTermsDTM[x],coocSentenceSignificance[x[1],x[2]])))
  } else {
    topCoocs <- t(apply(topCoocsIdx, 1, FUN=function(x) c(dictTermsDTM[x],coocSentenceSignificance[x[1],x[2]],dtm2[x[1],x[2]],   dtm2[x[1],x[2]] * 100 / coocSentenceSignificance[x[1],x[2]])))
  }
  
  topCoocs <- topCoocs[seq(1,(N*2),2),]
  return(topCoocs)  
}


#####################################################
## functions to compute initial assigment-clusters ##
#####################################################

computeInitialAssignments <- function(corpusData, corpusVocab, K, nCoocs = 1000, clusterMethod = "PAM", cooccurrenceMeasure = "LL") {
  
  binDTM <- ldaformat2Matrix(corpusData, corpusVocab)
  binDTM[binDTM > 1] <- 1
  topNTerms <- colnames(binDTM)[order(colSums(binDTM), decreasing = T)[1:nCoocs]]
  binDTM <- binDTM[, topNTerms]
  corpusTTM <- t(binDTM) %*% binDTM
  corpusTTM[corpusTTM < 3] <- 0
  
  if (cooccurrenceMeasure == "DICE") {
    corpusTTMsig <- computeSignificanceDice(as.matrix(corpusTTM))
    corpusTTMsig[corpusTTMsig < 0.005] <- 0
  } else if (cooccurrenceMeasure == "LL") {
    corpusTTMsig <- computeSignificanceLL(as.matrix(corpusTTM), nrow(binDTM))
    corpusTTMsig[corpusTTMsig < 6.63] <- 0
  } else if (cooccurrenceMeasure == "POISSON") {
    DTM <- ldaformat2Matrix(corpusData, corpusVocab)
    topNTerms <- colnames(DTM)[order(colSums(DTM), decreasing = T)[1:nCoocs]]
    DTM <- DTM[, topNTerms]
    corpusTTM <- t(DTM) %*% DTM
    diag(corpusTTM) <- colSums(DTM)
    corpusTTM[corpusTTM < 3] <- 0
    corpusTTMsig <- computeSignificancePoisson(as.matrix(corpusTTM), rowSums(DTM), pLevel = 0.05)
  } else {
    stop(paste0("cooccurrence significance measure ", cooccurrenceMeasure, " is undefined"))
  }
  
  clusteredTerms <- clusterTerms(corpusTTMsig, K, method = clusterMethod)
  set.seed(1000)
  vocabClusterAssignments <- sample.int(K, size = length(corpusVocab), replace=TRUE) - 1
  as.numeric(Sys.time())-> t; set.seed((t - floor(t)) * 1e8 -> seed)
  
  for (k in 1:K) {
    terms <- topNTerms[clusteredTerms == k]
    vocabClusterAssignments[corpusVocab %in% terms] <- k - 1
  }
  initialAssignments <- lapply(corpusData, FUN = function(x) {vocabClusterAssignments[(x[1, ] + 1)]})
  initialAssignments <- lapply(initialAssignments, as.integer) 
  return(initialAssignments)

}


clusterTerms <- function(ttm, K, method = "PAM"){
  print("Clustering terms (may take a while)")
  
  topNTerms <- colnames(ttm)
  
  if (method == "PAM") {
    
    require(cluster)
    pamResult <- pam(ttm, K)
    t <- table(pamResult$clustering)
    t[t > 1]
    return(pamResult$clustering)
  
  } else if (method == "KMST") {
    
    require("igraph")
    g <- as.undirected(graph.adjacency(ttm * -1, mode="undirected", weighted=TRUE))
    g_mst <- mst(g, method="prim")
    g_mst <- delete_edges(g_mst, order(E(g_mst)$weight)[1:(K-1)])
    return(components(g_mst)$membership)
  
  } else if (method == "INFOMAP") {
    
    require("igraph")
    g <- as.undirected(graph.adjacency(ttm, mode="undirected", weighted = TRUE))
    g_imap <- infomap.community(g)
    t <- table(g_imap$membership)
    t[t > 1]
    sapply(1:K, FUN=function(x) topNTerms[g_imap$membership == x])
    return(g_imap$membership)
    
  } else {
    stop(paste0("clustering method ", method, " is undefined"))
  }
  
}


#####################################################
## functions to run model evaluation ################
#####################################################

runLDAModelEvalution <- function(sampleDirectory, K, iterations, alphaPriorsToTest, etaPriorsToTest, blacklist = NULL, initTopicAssignments = T, clusterMethod = "PAM", cooccurrenceMeasure = "LL", nCoocs = 2500) {
  corpusData <- read.documents(filename = paste0(sampleDirectory, "/corpus/corpus"))
  corpusVocab <- readLines(paste0(sampleDirectory, "/corpus/vocabulary.dat"), encoding = "UTF-8")
  
  if (!is.null(blacklist)) {
    corpusDTM <- ldaformat2dtm(corpusData, corpusVocab)
    corpusDTM <- corpusDTM[, !(colnames(corpusDTM) %in% blacklist)]
    corpusLDA <- dtm2ldaformat(corpusDTM)
    corpusData <- corpusLDA$documents
    corpusVocab <- corpusLDA$vocab
  }
  
  if (initTopicAssignments) {
    initAssignments <- computeInitialAssignments(corpusData, corpusVocab, K, nCoocs = nCoocs, clusterMethod = clusterMethod, cooccurrenceMeasure = cooccurrenceMeasure)
    initAssignments <- list(assignments = initAssignments)
    initStrategy <- "fix"
  } else {
    initAssignments <- NULL
    initStrategy <- "rnd"
  }

  modelEvaluationData <- data.frame(modelID = integer(), K = integer(), alpha = double(), eta = double(), modelLikelihood = double(), modelCoherence = integer(), topicsHTML = character(), modelFileName = character())
  modelEvaluationHTML <- data.frame()
  modelID <- 0
  for (alpha in alphaPriorsToTest) {
    for (eta in etaPriorsToTest) {
      modelID <- modelID + 1
      set.seed(123)
	  topicModel <- lda.collapsed.gibbs.sampler(corpusData, K, vocab = corpusVocab, num.iterations = iterations, alpha = alpha, eta = eta, initial = initAssignments, trace = 1L, compute.log.likelihood = T)
      topicTerms <- top.topic.words(topicModel$topics, 25, by.score=TRUE)
      topicProportions <- t(topicModel$document_sums) / colSums(topicModel$document_sums)
      modelLikelihood <- tail(as.vector(topicModel$log.likelihoods[2, ]), 1)
      topicCoherenceForAllTopics <- topicCoherence(ldaformat2Matrix(corpusData, corpusVocab), topicModel)
      modelCoherence <- mean(topicCoherenceForAllTopics)
      
      tProportions <- colSums(topicProportions) / nrow(topicProportions)
      oProportions <- order(tProportions, decreasing = T)
      firstDocTopics <- apply(topicProportions, 1, FUN=function(x) order(x, decreasing=TRUE)[1])
      primaryDocTopics <- factor(firstDocTopics, 1:K)
      nRanks1 <- table(primaryDocTopics)
      tRanks1 <- as.integer(nRanks1)
      oRanks1 <- order(tRanks1, decreasing = T)
      tCoherence <- topicCoherenceForAllTopics
      oCoherence <- order(topicCoherenceForAllTopics, decreasing = T)
      topics <- data.frame(
        TopicID = sprintf("%02d", 1:K),
        Rank1 = tRanks1,
        Prob = tProportions,
        Coherence = tCoherence,
        Terms = paste0("<pre>", apply(topicTerms, 2, paste0, collapse=" "), "</pre>")
      )
      topics <- topics[oRanks1, ]
      topicsHTML <- print(xtable(topics, digits=5), print.results = F, type="html", sanitize.text.function=function(x){x}, include.rownames = F, html.table.attributes = paste0('id="T', modelID, '" class="sortable"'))

      modelFileName <- paste0(storageDirectory, "/", sampleName, "_model-", modelID, 
                              "_K", K, "_a", alpha, "_e", eta, "_i", iterations, initStrategy,
                              "_d", format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), ".RData")
      save(topicModel, file = modelFileName)
      
      modelEvaluationHTML <- rbind(modelEvaluationHTML, data.frame(modelID, K, alpha, eta, modelLikelihood, modelCoherence, topicsHTML, modelFileName))
      print(paste0("LDA Model ", modelID, " has been saved to ", modelFileName))
    }
  }
  evalFileName <- paste0(storageDirectory, "/", sampleName,"_modelEvaluationData_K", K, "_i", iterations, initStrategy, ".html")
  file.copy("html_template.html", evalFileName, overwrite = T)
  print(xtable(modelEvaluationHTML, digits=4), type="html", file=evalFileName, append=T, sanitize.text.function=function(x){x}, include.rownames = F)
  
  print(paste0("Final LDA model evaluation CSV data has been has been saved to ", modelFileName))
  
  return(modelEvaluationHTML)
}


topicCoherence <- function(DTM, ldaModel, N = 25) {
  
  if (is.simple_triplet_matrix(DTM)) {
    DTM <- sparseMatrix(i=DTM$i, j=DTM$j, x=DTM$v, dims=c(DTM$nrow, DTM$ncol), dimnames = dimnames(DTM))
  }
  
  DTMBIN <- DTM
  DTMBIN[DTMBIN > 0] <- 1
  
  documentFrequency <- colSums(DTMBIN)
  names(documentFrequency) <- colnames(DTMBIN)
  
  K <- nrow(ldaModel$topics)
  
  topNtermsPerTopic <- top.topic.words(ldaModel$topics, N, by.score=TRUE)
  allTopicModelTerms <- unique(as.vector(topNtermsPerTopic))
  
  DTMpreprocessed <- DTMBIN[, allTopicModelTerms]
  DTMpreprocessedCooc <- t(DTMpreprocessed) %*% DTMpreprocessed
  DTMpreprocessedCooc <- t((DTMpreprocessedCooc + 1) / colSums(DTMpreprocessed))
  DTMpreprocessedCooc <- log(DTMpreprocessedCooc)
  DTMpreprocessedCooc <- as.matrix(DTMpreprocessedCooc)
  
  coherence <- rep(0, K)
  pb <- txtProgressBar(max = K)
  for (topicIdx in 1:K) {
    setTxtProgressBar(pb, topicIdx)
    topWordsOfTopic <- topNtermsPerTopic[, topicIdx]
    coherence[topicIdx] <- 0
    for (m in 2:length(topWordsOfTopic)) {
      for (l in 1:(m-1)) {
        mTerm <- topWordsOfTopic[m]
        lTerm <- topWordsOfTopic[l]
        coherence[topicIdx] <- coherence[topicIdx] + DTMpreprocessedCooc[mTerm, lTerm]
      }
    }
  }
  close(pb)
  return(coherence)
}


ldaformat2Matrix <- function (documents, vocab) {
  stm <- simple_triplet_matrix(
    i = rep(seq_along(documents), sapply(documents, ncol)), 
    j = as.integer(unlist(lapply(documents,"[", 1, )) + 1L), 
    v = as.integer(unlist(lapply(documents,"[", 2, ))), 
    nrow = length(documents), 
    ncol = length(vocab), 
    dimnames = list(names(documents), vocab))
  dtm <- sparseMatrix(i=stm$i, j=stm$j, x=stm$v, dims=c(stm$nrow, stm$ncol), dimnames = dimnames(stm))
}


peekIntoModelCorpus <- function(sampleCorpusFulltext, topicModel, topicToInvestigate, topicThresholdInDocument = NULL, n = 1) {
  topicProportions <- t(topicModel$document_sums) / colSums(topicModel$document_sums)
  if (!is.null(topicThresholdInDocument)) {
    idx <- which(topicProportions[, topicToInvestigate] > topicThresholdInDocument)
    docSampledIds <- sample(idx, min(n, length(idx)), n)
  } else {
    docSampledIds <- order(topicProportions[, topicToInvestigate], decreasing = T)[1:n]
  }
  sampleTexts <- sampleCorpusFulltext[docSampledIds]
  html_print(HTML(paste0(sampleTexts, collapse = "<hr><br><br/><br/>")))
}


getModelAlignments <- function(modelEvaluationData, topWordsToMatch = 100, similarityThreshold = 0.2, verbose = F) {
  numModels <- nrow(modelEvaluationData)
  if (numModels < 2) stop("Nothing to compare, got just than one model!")
  cat(c("Parameters: nWords", topWordsToMatch, "| threshold", similarityThreshold, "\n"))
  pairs <- combn(as.character(modelEvaluationData$modelFileName), 2, simplify = F)
  allReliabilities <- rep(0, length(pairs))
  i <- 0
  for (pair in pairs) {
    i <- i + 1
    cat(c("------", "\n"))
    cat(c(pair[1], "\n"))
    cat(c(pair[2], "\n"))
    tm1 <- get(load(file = pair[1]))
    tm2 <- get(load(file = pair[2]))
    alignment <- alignTopicModels(tm1, tm2, topWordsToMatch, similarityThreshold)
    printAlignedTopics(alignment, verbose = verbose)
    allReliabilities[i] <- alignment$reliability
  }
  return(allReliabilities)
}


#####################################################
## functions to compare topic models ################
#####################################################

topicProbability <- function(topics) {
  token <- sum(topics)
  prob <- c()
  for(i in 1:nrow(topics)) {
    ts <- sum(topics[i, ]) / token
    prob <- c(prob, ts)
  }
  return(prob)
}


toComparable <- function(topics, topWordsToMatch) {
  for(i in 1:nrow(topics)) {
    ts <- topics[i, ] / sum(topics[i, ])
    left <- names(ts[order(-ts)][1:topWordsToMatch])
    topics[i, ] <- 0
    topics[i, left] <- ts[left]
  }
  return(topics)
}


TM_Aligner <- function(topics1, topics2, thres, probs1, probs2, topWordsToMatch) {
  K <- nrow(topics1)
  reliability <- 0
  cosineDists <- as.matrix(1 - topics1 %*% t(topics2) / (sqrt(rowSums(topics1 ^ 2) %*% t(rowSums(topics2 ^ 2)))))
  minIndexes <- apply(cosineDists, 2, which.min)
  mins <- apply(cosineDists, 2, min)
  
  alignedTopics <-list()
  alignedTopics$ids <- matrix(0, nrow = K, ncol = 2)
  alignedTopics$probabilities <- matrix(0, nrow = K, ncol = 2)
  alignedTopics$sharedTerms <- vector("list", K)
  alignedTopics$distance <- rep(0, K)
  
  for (i in 1:K) {
    
    index <- arrayInd(which.min(cosineDists), .dim = c(K, K))
    value <- min(cosineDists)
    if (value > thres)
      break
    
    cosineDists[index[1], ] <- 1
    cosineDists[, index[2]] <- 1
    
    alignedTopics$ids[i, ] <- c(index[1], index[2])
    alignedTopics$probabilities[i, ] <- c( probs1[index[1]],  probs1[index[2]])
    alignedTopics$sharedTerms[[i]] <- intersect(
      names(topics1[index[1], ][order(-topics1[index[1], ])][1:topWordsToMatch]), 
      names(topics2[index[2], ][order(-topics2[index[2], ])][1:topWordsToMatch]))
    
    reliability = reliability + 1
    
    alignedTopics$distance[i] <- value
    
  }
  alignedTopics$reliability <- reliability / K
  return(alignedTopics)
}


alignTopicModels <- function(tm1, tm2, topWordsToMatch = 50, similarityThreshold = 0.2) {
  c_topics1 <- toComparable(tm1$topics, topWordsToMatch)
  c_topics_p1 <- topicProbability(tm1$topics)
  c_topics2 <- toComparable(tm2$topics, topWordsToMatch)
  c_topics_p2 <- topicProbability(tm2$topics)
  alignedTopics <- TM_Aligner(c_topics1, c_topics2, similarityThreshold, c_topics_p1, c_topics_p2, topWordsToMatch)
  return(alignedTopics)
}


printAlignedTopics <- function(alignedTopics, verbose = F) {
  if (verbose) {
    for (i in 1:length(alignedTopics$sharedTerms)) {
      if (length(alignedTopics$sharedTerms[[i]]) > 0) {
        cat(c("___________________________________________________________","\n"))
        cat(c("Shared terms:", alignedTopics$sharedTerms[[i]], "\n"))
        cat(c("Distance:", sprintf("%.4f", alignedTopics$distance[i]),"\n"))
        cat(c("Alignment:", alignedTopics$ids[i, 1], "TO", alignedTopics$ids[i, 2], "\n"))
        cat(c("Probabilities:", sprintf("%.4f", alignedTopics$probabilities[i, 1]), "TO:", sprintf("%.4f", alignedTopics$probabilities[i, 2]), "\n"))
      }
    }
    cat(c("===========================================================", "\n"))
  }
  cat(c("RELIABILITY:", sprintf("%.4f", alignedTopics$reliability), "\n"))
}

RUN_AS_MAIN <- FALSE
if (RUN_AS_MAIN) {
  tm1 <- get(load(file = "../online_climate_uk/model-1_K40_a0.002_e0.025_i10000rnd_d2015-12-18_16-28-51.RData"))
  tm2 <- get(load(file = "../online_climate_uk/model-1_K40_a0.002_e0.025_i10000fix_d2015-12-16_01-55-11.RData"))
  alignment <- alignTopicModels(tm1, tm2, topWordsToMatch = 100, similarityThreshold = 0.2)
  printAlignedTopics(alignment)
}

#####################################################
## run topic models #################################
#####################################################


#--------------------------------------------------
# new model
#--------------------------------------------------

blacklistBadTokens <- readLines("/data/blacklist_badtokens.txt", encoding = "UTF-8")
tokensToIgnore <- blacklistBadTokens

sampleName <- "online_food_us"
storageDirectory <- getwd() # or different set sorage directory
sampleDirectory <- paste0("data/", sampleName)

K <- 50
iterations <- 1000
alphaPriorsToTest <- c(0.1, 0.05, 0.1, 0.2, 0.5, 1)
etaPriorsToTest <- c(1 / K)

# ------------------------------------------------

modelEvaluationData <- runLDAModelEvalution(sampleDirectory, K, iterations, alphaPriorsToTest, etaPriorsToTest, 
                                            blacklist = tokensToIgnore, 
											initTopicAssignments = T, 
                                            clusterMethod = "PAM", cooccurrenceMeasure = "LL", nCoocs = 3000)

											
K <- 30
modelEvaluationData <- runLDAModelEvalution(sampleDirectory, K, iterations, alphaPriorsToTest, etaPriorsToTest, 
                                            blacklist = tokensToIgnore, 
											initTopicAssignments = T, 
                                            clusterMethod = "PAM", cooccurrenceMeasure = "LL", nCoocs = 3000)

K <- 70
modelEvaluationData <- runLDAModelEvalution(sampleDirectory, K, iterations, alphaPriorsToTest, etaPriorsToTest, 
                                            blacklist = tokensToIgnore, 
											initTopicAssignments = T, 
                                            clusterMethod = "PAM", cooccurrenceMeasure = "LL", nCoocs = 3000)											
