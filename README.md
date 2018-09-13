# SimPhe
Bug report git repository for SimPhe

*SimPhe* is an R package dedicated to simulate multiple phenotyes based on genotyping data. 
The main feature of this package is the possibility to take genetic epistatic effects, 
not only additive-additive interaction but also additive-dominance and dominance-dominance, 
heritability and multiple phenotypes into account. Moreover, we provide a variety of convenient functions, 
including suggestion for set the variation of random variable according to user specified genetic effects 
and flexible support for input formats, as well as output formats. It also supports the input of plink formats.

## Installation
```{r, eval = F, error = F, results='hide'}
install.packages(SimPhe)
```

## Sample implementation
```{r, echo = T, results='hide'}
# load package
library(SimPhe)
```
First, we show how easy to get phenotype(s) by using *sim.phe*. Before running *sim.phe*, we need specify the parameter file and genotype file for simulation. After installing *SimPhe*, these two toy files exist in package folder. To get the path, type
```{r, echo = T, results='hide'}
# get file path of simulation parameters
# (two shared SNP pairs and one independent SNP pair for each phenotype)
fpar.path <- system.file("extdata", "simupars.txt", package="SimPhe")

# get file path of genotype file: rows are individuals and columns are SNPs
fgeno.path <- system.file("extdata", "10SNP.txt", package="SimPhe")
```

Then simulate the phenotypes as designed in the parameter file after loading package:
```{r, echo = T, results='hide'}
# simulate phenotype(s)
phe <- sim.phe(sim.pars = fpar.path,
               fgeno = fgeno.path,
               ftype = "snp.head",
               seed = 123,
               fwrite = FALSE)
```

In the parameter file, we conduct two phenotypes contrbituted by two common SNP pairs with epistatic effects and one independent SNP pair with epistatic effects, the example of simulation parameters, same as the example file, is also included in *SimPhe*. Users could touch it by typing *genepars*.
```{r, echo = T, results='markup'}
genepars
```

Phenotype 1 has been set with certain heritability but phenotype 2 has not. We will show whether the heritability of the simulated phenotype 1 is same as the set. *SimPhe* includes the coefficients and the allele frequencies for simulating phenotype 1: *gene.coef* and *allele.freq*, which is extracted from the simulation parameters.
```{r, echo = T, results='markup'}
# regression coefficients for simulating phenotype 1
gene.coefficients
# allele frequencies of SNPs for simulating phenotype 1
allele.freq

# calculate genetic variance
genevar <- calc.gene.var(gene.coefficients, allele.freq)

# the variance of simulated phenotype 1
phe1var <- var(phe[, "p1"])

# heritability of simulated phenotype 1
simuht <- genevar / phe1var
simuht
```

The result is not the exactly 0.6 due to the (pseudo)random numbers generated in R.
To get phenotype 2 with a specific heritability, for example, 0.45, we
could proceed as:
```{r, echo = T, results='markup'}
# get regression coefficients settings for phenotype 2
genecoef <- get.gene.coef(main.pars = specify.pars(genetic.pars = genepars,
                                                   effect.type  = "main",
                                                   phe.index    = 2),
                          epi.pars  = specify.pars(genetic.pars = genepars,
                                                   effect.type  = "epistasis",
                                                   phe.index    = 2))
# get genotype data
genotype <- read.geno(fname = fgeno.path, ftype = "snp.head")

# get allele frequencies of the SNPs set for phenotype 2
freq2 <- get.freq(geno = genotype,
                  epi.pars  = specify.pars(genetic.pars = genepars,
                                           effect.type  = "epistasis",
                                           phe.index    = 2))

# get noise variation
exp.noise.var <- get.noise.var(gene.coef = genecoef,
                               freq = freq2,
                               heritability = 0.45)
```

Then when simulating a phenotype, just give this value as argument
*noise.var* to function *sim.phe*. It will generate a
phenotype which has a heritability close to 0.45.

As we mentioned earlier in this article, buiding the correlation by
setting the shared interactive SNP pairs cannot be controlled. We can take a look at the correlation between the simulated phenotype 1 and phenotype 2:
```{r, echo = T, results='markup'}
# correlation test
cor.test(phe[, "p1"], phe[, "p2"])
```

According to the result of correlation test, the two simulated
phenotypes are significantly correlated but the correlation
coefficient is small and the amount cannot easily be predicted. To get
a certain amount or higher correlation, we can conduct correlation by
applying the correlation matrix to two independent variables. For two
phenotypes, if we set different SNP pairs for each, we assume these
two phenotypes are independent. Here we give another parameter file
to but same genotype file to the arguments in *sim.phe*,
*fgenetic.pars* and *fgeno*:
```{r, echo = T, results='markup'}
# get file path of simulation parameters
# (different SNP pairs for each phenotype)
fpar.path <- system.file("extdata", "sep_simupars.txt", package="SimPhe")

# simulate phenotype(s)
indphe <- sim.phe(sim.pars = fpar.path,
                  fgeno = fgeno.path,
                  ftype = "snp.head",
                  seed = 123,
                  fwrite = FALSE)
```

We can test the correlation between initial phenotypes with seperated SNP pairs settings:
```{r, echo = T, results='markup'}
# test original correlation
cor.test(indphe[, "p1"], indphe[, "p2"])
```

Apparently, these two phenotypes are not related. We could continue our work on coverting them to be correlated. First, conduct a correlation matrix:
```{r, echo = T, results='markup'}
# correlation matrix
corm <- matrix(c(1, 0.6, 0.6, 1), ncol = 2)
corm
```

Before applying correlation matrix to simulated phenotypes, we would like to know what the data looks like:
```{r, echo = T, results='markup'}
apply(indphe, 2, mean)
apply(indphe, 2, sd)
```

Then we can build correlation between the two initial phenotypes:
```{r, echo = T, results='markup'}
# build correlation
corphe <- build.cor.phe(indphe, corMtr = corm)
```

Now let's test the correlation between the two new phenotypes and see if there is anything difference:
```{r, echo = T, results='markup'}
# check mean and standard deviation of new data set
apply(corphe, 2, mean)
apply(corphe, 2, sd)
```

```{r, echo = T, results='markup'}
# test correlation
cor.test(corphe[, "p1"], corphe[, "p2"])
```
Obviously, there is no significant difference on means and standard deviations between the initial phenotypes and the new phenotypes.


# References







---
- id: cockerham1977quadratic
  title: Quadratic analyses of reciprocal crosses
  author:
  - family: Cockerham
    given: C. Clark
  - family: Weir
    given: Bruce Spencer
  container-title: Biometrics
  volume: 33
  URL: 'http://www.jstor.org/stable/2529312'
  DOI: 10.2307/2529312
  issue: 1
  publisher: JSTOR
  page: 187-203
  type: article-journal
  issued:
    year: 1977
    month: 3
- id: kao2002modeling
  title: Modeling epistasis of quantitative trait loci using Cockerham's model
  author:
  - family: Kao
    given: Chen-Hung
  - family: Zeng
    given: Zhao-Bang
  container-title: Genetics
  volume: 160
  URL: 'http://www.genetics.org/content/160/3/1243'
  DOI: 10.1534/genetics.104.035857
  issue: 3
  publisher: Genetics Society of America
  page: 1243-1261
  type: article-journal
  issued:
    year: 2002
    month: 3
nocite: |
  @cockerham1977quadratic, @kao2002modeling
---
