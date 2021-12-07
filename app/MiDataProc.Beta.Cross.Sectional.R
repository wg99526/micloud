######################################
# Quality control and transformation #
######################################

library(phangorn)
library(phyloseq)
library(zCompositions)
library(plotly)
library(dplyr)
library(forestplot)
library(quantreg)
library(fossil)
library(picante)
library(entropart)
library(lme4)
library(lmerTest)
library(dirmult) 
library(robustbase)
library(robCompositions) 
library(BiasedUrn)
library(CompQuadForm)
library(GUniFrac) 
library(ecodist) 
library(MiRKAT)
library(GLMMMiRKAT)
library(proxy)


rem.tax.d <- c("", "gut metagenome", "mouse gut metagenome")
rem.tax.str.d <- c("uncultured", "Incertae", "unidentified")

num.tax.rank <- function(tax.tab, na.code = "NANANA") {
  tax.tab.cleaned <- tax.tab.clean(tax.tab, na.code = na.code)
  num.taxa <- c()
  for (i in 1:6) {
    taxa <- unique(tax.tab.cleaned[,i+1])
    uni.taxa <- sum(taxa == na.code)
    num.taxa[i] <- nrow(taxa) - uni.taxa
  }
  return(num.taxa)
}

lib.size.func <- function(biom) {
  otu.tab <- otu_table(biom)
  lib.size <- colSums(otu.tab)
  lib.size.sum <- c(mean(lib.size), quantile(lib.size))
  names(lib.size.sum) <- c("Mean", "Minimum", "1st quartile", "Median", "3rd quartile", "Maximum")
  return(list(lib.size = lib.size, lib.size.sum = lib.size.sum, num.sams = ncol(otu.tab), num.otus = nrow(otu.tab)))
}

mean.prop.func <- function(biom) {
  otu.tab <- otu_table(biom)
  lib.size <- colSums(otu.tab)
  prop.otu.tab <- otu.tab
  for (i in 1:length(lib.size)) {
    prop.otu.tab[,i] <- otu.tab[,i]/lib.size[i]
  }
  mean.prop <- rowMeans(prop.otu.tab)
  mean.prop.sum <- c(mean(mean.prop), quantile(mean.prop))
  names(mean.prop.sum) <- c("Mean", "Minimum", "1st quartile", "Median", "3rd quartile", "Maximum")
  return(list(mean.prop = mean.prop, mean.prop.sum = mean.prop.sum, num.sams = ncol(otu.tab), num.otus = nrow(otu.tab)))
}

rarefy.func <- function(biom, cut.off, multi.rarefy = FALSE) {
  
  if (!multi.rarefy | multi.rarefy == 1) {
    biom <- rarefy_even_depth(biom, cut.off, rngseed = 487)
  } else {
    otu.tab <- otu_table(biom)
    tax.tab <- tax_table(biom)
    tree <- phy_tree(biom)
    sam.dat <- sample_data(biom)
    
    otu.tab.list <- list()
    for (i in 1:multi.rarefy) {
      otu.tab.list[[i]] <- otu_table(rarefy_even_depth(biom, cut.off, rngseed = i), taxa_are_rows = TRUE)
    }
    
    sum.otu.tab <- otu.tab.list[[1]]
    for (i in 2:multi.rarefy) {
      sum.otu.tab <- sum.otu.tab + otu.tab.list[[i]]
    }
    otu.tab <- otu_table(round(sum.otu.tab/multi.rarefy), taxa_are_rows = TRUE)
    biom <- merge_phyloseq(otu.tab, tax.tab, tree, sam.dat) 
  }
  
  return(biom)
}

tax.trans <- function(otu.tab, tax.tab, rare.otu.tab, rare.tax.tab, sub.com = TRUE, na.code = "NANANA") {
  
  n <- ncol(otu.tab)
  lib.size <- colSums(otu.tab)
  rare.n <- ncol(rare.otu.tab)
  
  tax.count.out <- list()
  tax.rare.count.out <- list()
  tax.prop.out <- list()
  tax.imp.prop.out <- list()
  tax.clr.out <- list()
  tax.sub.clr.out <- list()
  
  for (j in 1:6) {
    tax <- as.vector(unique(tax.tab[,j+1]))
    tax.count <- matrix(NA, n, length(tax))
    for (i in 1:length(tax)) {
      ind.tax <- which(tax.tab[,j+1] == tax[i])
      tax.count[,i] <- colSums(otu.tab[ind.tax,])
    }
    rownames(tax.count) <- colnames(otu.tab)
    colnames(tax.count) <- tax
    
    tax.prop <- matrix(NA, n, length(tax))
    for (i in 1:length(lib.size)) {
      tax.prop[i,] <- tax.count[i,]/lib.size[i]
    }
    rownames(tax.prop) <- colnames(otu.tab)
    colnames(tax.prop) <- tax
    
    tax.imp.prop <- zCompositions::cmultRepl(tax.count)
    tax.clr <- compositions::clr(tax.imp.prop)
    
    rare.tax <- as.vector(unique(rare.tax.tab[,j+1]))
    tax.rare.count <- matrix(NA, rare.n, length(rare.tax))
    for (i in 1:length(rare.tax)) {
      ind.tax <- which(rare.tax.tab[,j+1] == rare.tax[i])
      tax.rare.count[,i] <- colSums(rare.otu.tab[ind.tax,])
    }
    rownames(tax.rare.count) <- colnames(rare.otu.tab)
    colnames(tax.rare.count) <- rare.tax
    
    ind <- which(tax == na.code)
    tax.count.out[[j]] <- as.data.frame(tax.count[,-ind])
    tax.rare.count.out[[j]] <- as.data.frame(tax.rare.count[,-ind])
    tax.prop.out[[j]] <- as.data.frame(tax.prop[,-ind])
    tax.imp.prop.out[[j]] <- tax.sub.imp.prop <- as.data.frame(tax.imp.prop[,-ind])
    tax.clr.out[[j]] <- as.data.frame(tax.clr[,-ind])
    tax.sub.clr.out[[j]] <- as.data.frame(compositions::clr(tax.sub.imp.prop))
  }
  
  names(tax.count.out) <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
  names(tax.rare.count.out) <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
  names(tax.prop.out) <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
  names(tax.imp.prop.out) <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
  names(tax.clr.out) <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
  names(tax.sub.clr.out) <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
  
  if (sub.com) {
    return(list(count = tax.count.out, rare.count = tax.rare.count.out, prop = tax.prop.out, clr = tax.sub.clr.out))
  }
  if (!sub.com) {
    return(list(count = tax.count.out, rare.count = tax.rare.count.out, prop = tax.prop.out, clr = tax.clr.out))
  }
  
}

#####################
# Data manipulation #
#####################

# is.mon.rev.bin.con <- function(sam.dat) {
#   
#   n.var <- ncol(sam.dat)
#   n.sam <- nrow(sam.dat)
#   is.mon <- logical()
#   is.rev <- logical()
#   is.bin <- logical()
#   is.con <- logical()
#   
#   for (i in 1:n.var) {
#     sam.var <- as.matrix(sam.dat[,i])
#     if (length(table(sam.var)) == 1) {
#       is.mon[i] <- TRUE
#     }
#     if (length(table(sam.var)) != 1) {
#       is.mon[i] <- FALSE
#     }
#     if (length(table(sam.var)) == n.sam & sum(is.na(as.numeric(sam.var))) == n.sam) {
#       is.rev[i] <- TRUE
#     }
#     if (length(table(sam.var)) != n.sam | sum(is.na(as.numeric(sam.var))) != n.sam) {
#       is.rev[i] <- FALSE
#     }
#     if (length(table(sam.var)) == 2) {
#       is.bin[i] <- TRUE
#     }
#     if (length(table(sam.var)) != 2) {
#       is.bin[i] <- FALSE
#     }
#     if (length(table(sam.var)) != 2 & sum(is.na(as.numeric(sam.var))) != n.sam) {
#       is.con[i] <- TRUE
#     }
#     if (length(table(sam.var)) == 2 & sum(is.na(as.numeric(sam.var))) != n.sam) {
#       is.con[i] <- FALSE
#     }
#     if (sum(is.na(as.numeric(sam.var))) == n.sam) {
#       is.con[i] <- FALSE
#     }
#     
#   }
#   return(list(is.mon = is.mon, is.rev = is.rev, is.bin = is.bin, is.con = is.con))
# }
# 
# pri.func <- function(sam.dat, mon.rev.bin.con) {
#   colnames(sam.dat)[(mon.rev.bin.con$is.bin | mon.rev.bin.con$is.con) & !mon.rev.bin.con$is.mon]
# }
# 
# is.bin.con.pri <- function(sam.dat, mon.rev.bin.con, sel.pri.var) {
#   ind <- which(colnames(sam.dat) == sel.pri.var)
#   if(length(ind) != 0){
#     if (mon.rev.bin.con$is.bin[ind]) {
#       out <- "Binary"
#     } else {
#       out <- "Continuous"
#     }
#   }else {
#     out = "Neither"
#   }
#   return(out)
# }
# 
# cov.func <- function(sam.dat, mon.rev.bin.con, sel.pri.var) {
#   ind.pri <- colnames(sam.dat) == sel.pri.var
#   ind.mon.rev <- mon.rev.bin.con$is.mon | mon.rev.bin.con$is.rev
#   return(colnames(sam.dat)[!(ind.pri | ind.mon.rev)])
# }

beta.bin.cat.func <- function(sam.dat, sel.bin.var) {
  bin.var <- unlist(sam.dat[,sel.bin.var])
  bin.var.no.na <- bin.var[!is.na(bin.var)]
  bin.cat <- unique(bin.var.no.na)
  return(bin.cat)
}

beta.bin.cat.ref.ori.func <- function(sam.dat, sel.bin.var = "ecig_status") {
  return(levels(as.factor(as.data.frame(as.matrix(sam.dat))[,sel.bin.var])))
}

beta.bin.cat.ref.func <- function(sel.bin.var, sel.ref, sel.com, sam.dat, Ds.Ks) {
  bin.var <- unlist(sam.dat[,sel.bin.var])
  ind.ref <- which(bin.var == sel.ref)
  ind.com <- which(bin.var == sel.com)
  bin.var <- factor(c(bin.var[ind.ref], bin.var[ind.com]), levels = c(sel.ref, sel.com))
  Ds <- Ds.Ks$Ds
  Ks <- Ds.Ks$Ks
  for (i in 1:length(Ds)) {
    Ds[[i]] <- Ds[[i]][c(ind.ref, ind.com), c(ind.ref, ind.com)]
  }
  for (i in 1:length(Ks)) {
    Ks[[i]] <- Ks[[i]][c(ind.ref, ind.com), c(ind.ref, ind.com)]
  }
  return(list(bin.var = bin.var, Ds = Ds, Ks = Ks))
}

beta.bin.cat.recode.func <- function(sam.dat, sel.bin.var = "ecig_status", ori.cat, rename.ref, rename.com) {
  ind.ref <- which(sam.dat[,sel.bin.var] == ori.cat[1])
  ind.com <- which(sam.dat[,sel.bin.var] == ori.cat[2])
  sam.dat[ind.ref,sel.bin.var] <- rename.ref
  sam.dat[ind.com,sel.bin.var] <- rename.com
  return(sam.dat)
}

beta.bin.cat.func <- function(sam.dat, sel.bin.var) {
  bin.var <- unlist(sam.dat[,sel.bin.var])
  bin.var.no.na <- bin.var[!is.na(bin.var)]
  bin.cat <- unique(bin.var.no.na)
  return(bin.cat)
}

##################
# Beta diversity #
##################

Ds.Ks.func <- function(rare.biom, biom.after.qc) {
  rare.otu.tab <- otu_table(rare.biom)
  no.rare.otu.tab <- otu_table(biom.after.qc)
  no.rare.tree <- phy_tree(biom.after.qc)
  
  jac <- as.matrix(proxy::dist(t(rare.otu.tab), method = "Jaccard"))
  bc <- as.matrix(bcdist(t(rare.otu.tab)))
  unifs <- GUniFrac(t(no.rare.otu.tab ), no.rare.tree, alpha = c(0.5, 1))$unifracs
  u.unif <- unifs[, , "d_UW"]
  g.unif <- unifs[, , "d_0.5"]
  w.unif <- unifs[, , "d_1"]
  
  jac.k <- D2K(jac)
  bc.k <- D2K(bc)
  u.unif.k <- D2K(u.unif)
  g.unif.k <- D2K(g.unif)
  w.unif.k <- D2K(w.unif)
  
  rownames(jac.k) <- colnames(jac.k) <- colnames(rare.otu.tab)
  rownames(bc.k) <- colnames(bc.k) <- colnames(rare.otu.tab)
  rownames(u.unif.k) <- colnames(u.unif.k) <- colnames(rare.otu.tab)
  rownames(g.unif.k) <- colnames(g.unif.k) <- colnames(rare.otu.tab)
  rownames(w.unif.k) <- colnames(w.unif.k) <- colnames(rare.otu.tab)
  
  return(
    list(Ds = list(Jaccard = jac, Bray.Curtis = bc, U.UniFrac = u.unif, G.UniFrac = g.unif, W.UniFrac = w.unif),
         Ks = list(Jaccard = jac.k, Bray.Curtis = bc.k, U.UniFrac = u.unif.k, G.UniFrac = g.unif.k, W.UniFrac = w.unif.k))
  )
}

beta.bin.cat.func <- function(sam.dat, sel.bin.var) {
  bin.var <- unlist(sam.dat[,sel.bin.var])
  bin.var.no.na <- bin.var[!is.na(bin.var)]
  bin.cat <- unique(bin.var.no.na)
  return(bin.cat)
}

beta.bin.cov.cat.ref.func <- function(sel.bin.var, sel.ref, sel.com, sel.cov.var, sam.dat, Ds.Ks) {  
  bin.var <- unlist(sam.dat[,sel.bin.var])
  ind.ref <- which(bin.var == sel.ref)
  ind.com <- which(bin.var == sel.com)
  bin.var <- factor(c(bin.var[ind.ref], bin.var[ind.com]), levels = c(sel.ref, sel.com))
  cov.var <- sam.dat[,sel.cov.var]
  cov.var <- as.data.frame(rbind(cov.var[ind.ref,], cov.var[ind.com,]))
  Ds <- Ds.Ks$Ds
  Ks <- Ds.Ks$Ks
  for (i in 1:length(Ds)) {
    Ds[[i]] <- Ds[[i]][c(ind.ref, ind.com), c(ind.ref, ind.com)]
  }
  for (i in 1:length(Ks)) {
    Ks[[i]] <- Ks[[i]][c(ind.ref, ind.com), c(ind.ref, ind.com)]
  }
  
  return(list(bin.var = bin.var, cov.var = cov.var, Ds = Ds, Ks = Ks))
}

beta.con.recode.func <- function(sam.dat, sel.con.var, rename.con.var, Ds.Ks) {
  ind.sel <- which(colnames(sam.dat) == sel.con.var)
  colnames(sam.dat)[ind.sel] <- rename.con.var
  ind.nona <- !is.na(sam.dat[,ind.sel])
  con.var <- as.data.frame(as.matrix(sam.dat[ind.nona,ind.sel]))
  Ds <- Ds.Ks$Ds
  Ks <- Ds.Ks$Ks
  for (i in 1:length(Ds)) {
    Ds[[i]] <- Ds[[i]][ind.nona, ind.nona]
  }
  for (i in 1:length(Ks)) {
    Ks[[i]] <- Ks[[i]][ind.nona, ind.nona]
  }
  return(list(con.var = con.var, Ds = Ds, Ks = Ks))
}

beta.con.cov.recode.func <- function(sam.dat, sel.con.var, sel.cov.var, rename.con.var, Ds.Ks) {
  ind.sel <- which(colnames(sam.dat) == sel.con.var)
  colnames(sam.dat)[ind.sel] <- rename.con.var
  ind.nona <- !is.na(sam.dat[,ind.sel])
  con.var <- as.data.frame(as.matrix(sam.dat[ind.nona,ind.sel]))
  cov.var <- as.data.frame(sam.dat[ind.nona, sel.cov.var])
  Ds <- Ds.Ks$Ds
  Ks <- Ds.Ks$Ks
  for (i in 1:length(Ds)) {
    Ds[[i]] <- Ds[[i]][ind.nona, ind.nona]
  }
  for (i in 1:length(Ks)) {
    Ks[[i]] <- Ks[[i]][ind.nona, ind.nona]
  }
  return(list(con.var = con.var, cov.var = cov.var, Ds = Ds, Ks = Ks))
}

#################
# Data analysis #
#################

### MiRKAT

mirkat.bin <- function(beta.bin.out) {
  set.seed(487)
  out <- MiRKAT(y = as.numeric(beta.bin.out$bin.var)-1, X = NULL, Ks = beta.bin.out$Ks, out_type = "D", nperm = 1000)
  
  par(mfrow = c(3, 2))
  for (i in 1:length(beta.bin.out$Ds)) {
    if (out$p_values[i] < 0.05) {
      sub.tit <- paste("*p:", p.value.0.1(out$p_values[i]), sep="")
    }
    if (out$p_values[i] >= 0.05) {
      sub.tit <- paste("p:", p.value.0.1(out$p_values[i]), sep="")
    }
    mod <- betadisper(as.dist(beta.bin.out$Ds[[i]]), beta.bin.out$bin.var)
    plot(mod, ellipse = TRUE, hull = FALSE, main = names(beta.bin.out$Ds)[i], xlab="PC 1", ylab="PC 2",
         sub = sub.tit, col = c("blue2", "red2"), cex=1.5)
  }
  
  plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
  legend("center", title = NULL, legend = levels(beta.bin.out$bin.var), fil = c("blue2", "red2", cex=2.5, box.lty=0), bty = "n", cex=1.5)
  legend("bottom", paste("Omnibus MiRKAT: ", p.value.0.1(out$omnibus_p), sep=""), bty = "n", cex=1.5)
  return(out)
}

mirkat.bin.cov <- function(beta.bin.cov.out) {
  set.seed(487)
  out <- MiRKAT(y = as.numeric(beta.bin.cov.out$bin.var)-1, X = as.matrix(beta.bin.cov.out$cov.var), Ks = beta.bin.cov.out$Ks, out_type = "D", nperm = 1000)
  
  par(mfrow = c(3, 2))
  for (i in 1:length(beta.bin.cov.out$Ds)) {
    if (out$p_values[i] < 0.05) {
      sub.tit <- paste("*p:", p.value.0.1(out$p_values[i]), sep="")
    }
    if (out$p_values[i] >= 0.05) {
      sub.tit <- paste("p:", p.value.0.1(out$p_values[i]), sep="")
    }
    mod <- betadisper(as.dist(beta.bin.cov.out$Ds[[i]]), beta.bin.cov.out$bin.var)
    plot(mod, ellipse = TRUE, hull = FALSE, main = names(beta.bin.cov.out$Ds)[i], xlab="PC 1", ylab="PC 2",
         sub = sub.tit, col = c("blue2", "red2"), cex=1.5)
  }
  
  plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
  legend("center", title = NULL, legend = levels(beta.bin.cov.out$bin.var), fil = c("blue2", "red2", cex=2.5, box.lty=0), bty = "n", cex=1.5)
  legend("bottom", paste("Omnibus MiRKAT: ", p.value.0.1(out$omnibus_p), sep=""), bty = "n", cex=1.5)
  return(out)
}

mirkat.con <- function(beta.con.out) {
  set.seed(487)
  out <- MiRKAT(y = as.numeric(unlist(beta.con.out$con.var)), X = NULL,  Ks = beta.con.out$Ks, out_type = "C", nperm = 1000)
  #print(out)
  par(mfrow = c(3, 2))
  for (i in 1:length(beta.con.out$Ds)) {
    if (out$p_values[i] < 0.05) {
      sub.tit <- paste("*p:", p.value.0.1(out$p_values[i]), sep="")
    }
    if (out$p_values[i] >= 0.05) {
      sub.tit <- paste("p:", p.value.0.1(out$p_values[i]), sep="")
    }
    con.var <- unlist(beta.con.out$con.var)
    con.var.med <- median(con.var)
    bin.var <- rep(NA, length(con.var))
    ind.gr <- which(con.var >= con.var.med)
    ind.sm <- which(con.var < con.var.med)
    bin.var[ind.gr] <- paste(names(beta.con.out$con.var), ">=", con.var.med)
    bin.var[ind.sm] <- paste(names(beta.con.out$con.var), "<", con.var.med)
    bin.var <- factor(bin.var)
    
    mod <- betadisper(as.dist(beta.con.out$Ds[[i]]), bin.var)
    plot(mod, ellipse = TRUE, hull = FALSE, main = names(beta.con.out$Ds)[i], xlab="PC 1", ylab="PC 2",
         sub = sub.tit, col = c("blue2", "red2"), cex=1.5)
  }
  
  plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
  legend("center", title = NULL, legend = levels(bin.var), fil = c("blue2", "red2", cex=2.5, box.lty=0), bty = "n", cex=1.5)
  legend("bottom", paste("Omnibus MiRKAT: ", p.value.0.1(out$omnibus_p), sep=""), bty = "n", cex=1.5)
  return(out)
}

mirkat.con.cov <- function(beta.con.cov.out) {
  set.seed(487)
  out <- MiRKAT(y = as.numeric(unlist(beta.con.cov.out$con.var)), X = as.matrix(beta.con.cov.out$cov.var),  Ks = beta.con.cov.out$Ks, out_type = "C", nperm = 1000)
  
  par(mfrow = c(3, 2))
  for (i in 1:length(beta.con.cov.out$Ds)) {
    if (out$p_values[i] < 0.05) {
      sub.tit <- paste("*p:", p.value.0.1(out$p_values[i]), sep="")
    }
    if (out$p_values[i] >= 0.05) {
      sub.tit <- paste("p:", p.value.0.1(out$p_values[i]), sep="")
    }
    con.var <- unlist(beta.con.cov.out$con.var)
    con.var.med <- median(con.var)
    bin.var <- rep(NA, length(con.var))
    ind.gr <- which(con.var >= con.var.med)
    ind.sm <- which(con.var < con.var.med)
    bin.var[ind.gr] <- paste(names(beta.con.cov.out$con.var), ">=", con.var.med)
    bin.var[ind.sm] <- paste(names(beta.con.cov.out$con.var), "<", con.var.med)
    bin.var <- factor(bin.var)
    
    mod <- betadisper(as.dist(beta.con.cov.out$Ds[[i]]), bin.var)
    plot(mod, ellipse = TRUE, hull = FALSE, main = names(beta.con.cov.out$Ds)[i], xlab="PC 1", ylab="PC 2",
         sub = sub.tit, col = c("blue2", "red2"), cex=1.5)
  }
  
  plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
  legend("center", title = NULL, legend = levels(bin.var), fil = c("blue2", "red2", cex=2.5, box.lty=0), bty = "n", cex=1.5)
  legend("bottom", paste("Omnibus MiRKAT: ", p.value.0.1(out$omnibus_p), sep=""), bty = "n", cex=1.5)
  return(out)
}

###################
# Other functions #
###################

q.func <- function(out, method = c("BH", "BY")) {
  Q.value <- p.adjust(out$P.value, method = method)
  return(cbind(out, Q.value))
}

p.value.0.1 <- function(x, round.x = 3) {
  x <- format(round(x, 3), digits = 3)
  ind.0 <- which(x == "0.000" | x == 0)
  x[ind.0] <- "<.001"
  ind.1 <- which(x == "1.000" | x == 1)
  x[ind.1] <- ">.999"
  return(x)
}