# Build a NJ tree
dnbin <- dist.dna(seqs, model = "K80")
tree <- nj(dnbin)


# HAPLOTYPE FREQUENCIES PER POPULATION -----------------------------------------  
kn <- as.dna(seqs)
kh <- haplotypes::haplotype(kn)
kh

# Population labels
ncb <- as.matrix(labels(seqs))
n2 <- NULL
for (i in 1:nrow(ncb)) n2[i] <- strsplit(ncb[i], '_')[[1]][1]
n2

# Haplotype frequency matrix per population
hf <- grouping(kh, factors = n2)
hf[["hapvec"]] <- NULL
dhf <- as.data.frame(hf$hapmat)
rownames(dhf) <- paste("H", 1:nrow(mat7), sep = "")
dhf


# EXTRACT SEQUENCE AND HAPLOTYPE INFORMATION -----------------------------------
# Format conversion
an <- as.alignment(seqs)           # converting DNAbin to alignment format
nm <- as.matrix(an)                # converting alignment to matrix
nbinmat <- as.matrix(labels(seqs)) # extraction of the sample names

nrow(nm) # 54 samples
ncol(nm) # 419 bp

sat2 <- NULL
for (i in 1:nrow(nm)) sat2[i] <- paste(nm[i, ], collapse="")
sat2 <- toupper(sat2)
sat3 <- unique(sat2)
sat3

hfreq <- NULL
for (i in 1:length(sat3)) {
  hcount = 0
  s3 <- sat3[i]
  for (j in 1:length(sat2)) {
    s2 <- sat2[j]
    if (s3 == s2) {
      hcount <- (hcount + 1) #counts the number of individuals with the same haplotype sequence. 
    }
  }
  hname<-(paste("H",i, sep =""))
  hfreq[i] <- hcount
}

len <- nchar(sat3[1]) #assume all have same length!!!
cnt <- 1
sat4 = list()
for (j in 1:len) {
  same <- TRUE
  first <- substr(sat3[1], j, j)
  for (i in 2:length(sat3)) {
    ch1 <- substr(sat3[i], j, j)
    if (first != ch1) {
      str <- paste(j, first, ch1)
      print(str)
      same <- FALSE
      break
    }
  }
  if (!same) {
    ss <- NULL
    for (i in 1:length(sat3)) {
      ss <- paste(ss, substr(sat3[i], j, j), sep="")
    }
    sat4[cnt] <- ss
    cnt <- cnt + 1
  }
}
len <- nchar(sat3[1]) #assume all have same length!!!
cnt <- 1
sat5 = list() 
for (j in 1:len) { #scan all columnns and if all elements are the same do not copy
  same <- TRUE
  first <- substr(sat3[1], j, j)
  scol <- first
  for (i in 2:length(sat3)) {
    ch1 <- substr(sat3[i], j, j)
    scol <- paste(scol, ch1, sep="")
    if (first != ch1) {
      str <- paste(j, first, ch1)
      same <- FALSE
    }
  }
  if (!same) {
    scol <- paste("V_", cnt, " ", scol, sep="")
    ss <- NULL
    for (i in 1:length(sat3)) {
      ss <- paste(ss, substr(sat3[i], j, j), sep="")
    } 
    sat5[cnt] <- ss
    cnt <- cnt + 1
  }
}

sat6 <- as.matrix(sat5)
mat6 <- matrix(nrow=nrow(sat6), ncol=nchar(sat6[1]))
for (i in 1:nrow(mat6)) {
  s <- as.vector(strsplit(as.character(sat5[i]), ""))
  for (j in 1:ncol(mat6)) {
    mat6[i, j] <- as.character(s[[1]][j])
  }
}
mat7 <- t(mat6)
hname <- paste("H", 1:nrow(mat7), sep = "")
rownames(mat7) <- hname

str4 <- NULL
str4[1] <- paste(mat7[1, ], collapse="")
for (i in 2:nrow(mat7)) {
  tmp <- NULL
  for (j in 1:ncol(mat7)) {
    chr = "."
    if(mat7[i, j] != mat7[1, j]) chr = mat7[i, j]
    tmp <- paste(tmp, chr, sep="")
  }
  str4[i] <- paste(tmp, collapse="")
}
nchar(str4[1]) #confirmation of number of variable sites
mstr4 <- as.matrix(str4)
rownames(mstr4) <- hname
colnames(mstr4) <- paste("sequences length","(", ncol(mat7), "base pairs", ")")
pct <- round((as.matrix(hfreq) * 100 / colSums(as.matrix(hfreq))), 2)
colnames(pct) <- c("pct")
cmstr4 <- as.data.frame(cbind(mstr4, hfreq, pct))
cmstr4


# DISTANCE MATRIX AND HEATMAP OF HAPLOTYPES ------------------------------------
# Distance matrix
dh <- dist.hamming(mat7)
dhm <- as.matrix(dh)
dhm

# Distance matrix for heatmap [SAME AS ABOVE?]
sat2 <- NULL
for (i in 1:nrow(nm)) sat2[i] <- paste(nm[i, ], collapse= "")
sat2 <- toupper(sat2)
sat3 <- unique(sat2)
comat <- matrix(nrow=length(sat3), ncol=length(sat3))
for (i in 1:length(sat3)) { 
  si <- sat3[i]
  for (j in 1:length(sat3)) { 
    sj <- sat3[j]
    difcnt = 0
    s1 = as.vector(strsplit(as.character(si), ""))
    s2 = as.vector(strsplit(as.character(sj), ""))
    for (k in 1:length(s1[[1]])) {
      if (s1[[1]][k] != s2[[1]][k]) {
        difcnt = difcnt + 1
      }
      comat[i, j] = difcnt
    }
  }
}
comat	#is Hamming distance matrix
colnames(comat) <- paste("H", 1:nrow(comat), sep = "")
rownames(comat) <- paste("H", 1:nrow(comat), sep = "")

# Plot the heatmap
heatmap(
  comat,
  scale="none",
  col = heat.colors(100),
  keep.dendro = TRUE,
  symm = TRUE
)


# INDIVIDUAL-LEVEL HAPLOTYPE NETWORKS ------------------------------------------
h <- pegas::haplotype(seqs, strict = FALSE, trailingGapsAsN = TRUE)
hname <- paste("H", 1:nrow(h), sep = "")
rownames(h) <- paste(hname)
net <- haploNet(h, d = NULL, getProb = TRUE)
net
ind.hap <- with(
  utils::stack(setNames(attr(h, "index"), rownames(h))),
  table(hap = ind, individuals = rownames(seqs))
)

par(mar = c(0.01, 0.01, 0.01, 15))

plot(
  net,
  size = attr(net, "freq"),
  scale.ratio = 2,
  cex = 0.6,
  labels = TRUE,
  pie = ind.hap,
  show.mutation = 1,
  font = 2,
  fast = TRUE
)

legend(
  x = 10,
  y = 15,
  colnames(ind.hap),
  fill = rainbow(ncol(ind.hap)),
  cex = 0.52,
  ncol=6,
  x.intersp=0.2,
  text.width = 11
)



# HAPLOTYPE TREE ---------------------------------------------------------------
# Build the tree
D <- dist.hamming(mat7)
class(D)
htre <- nj(D)
bp <- boot.phylo(htre, mat7, B = 100, function(x) nj(dist.hamming(x)))
bp2 <- data.frame(node = 1:Nnode(htre) + Ntip(htre), bootstrap = bp)
htree <- full_join(htre, bp2, by = "node")

# Plot the three
ggtree(htree, size = 1) +
  geom_tiplab(size = 4) +
  geom_nodepoint(
    aes(fill = cut(bootstrap, c(0, 50, 70, 85, 100))),
    shape = 21,
    size = 4
  ) +
  scale_fill_manual(
    values = c("black", "red", "pink1", "white"),
    guide = "legend",
    name = "Bootstrap Percentage (BP)",
    breaks = c("(85,100]", "(70,85]", "(50,70]", "(0,50]"),
    labels = expression(BP>=85, 70<=BP*"<85", 50<=BP*"<70", BP<50)
  ) +
  theme_tree(legend.position = c(0.75, 0.2))


# HAPLOTYPE TREE ---------------------------------------------------------------
#njmsaplot <- msaplot(
#  p = p1,
#  seqs,
#  offset = 0.02,
#  width = 1,
#  height = 0.7,
#  color = nuc_cols
#  )
#njmsaplot


#heatmap_mat_haplotype <- meta |>
#  select(sample, Haplotype) |>
#  column_to_rownames("sample")
#p2 <- gheatmap(
#  p = p1,
#  data = heatmap_mat_haplotype,
#  offset = 0.01,
#  width = 0.15,
#  colnames_position = "top",
#  colnames_offset_y = 2.8,
#  font.size = 2.5
#)
#p2
