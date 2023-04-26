---
title: 'iTensor: An R package for independent component analysis-based matrix/tensor decomposition'
tags:
  - R
  - independent component analysis
  - multimodal independent component analysis
  - group independent component analysis
  - multilinear independent component analysis
  - dimension reduction
authors:
  - name: Koki Tsuyuzaki^[first author]
    orcid: 0000-0003-3797-2148
    affiliation: "1, 2"
affiliations:
 - name: Department of Artificial Intelligence Medicine, Graduate School of Medicine, Chiba University, Japan
   index: 1
 - name: Laboratory for Bioinformatics Research, RIKEN Center for Biosystems Dynamics Research, Japan
   index: 2
date: 1 May 2023
bibliography: paper.bib
---

# Summary

Independent Component Analysis (ICA) is a widely used algorithm to extract a small number of mutually independent source signals in high-dimensional data. There are many applications of ICA in signal processing [@ica1; @ica2], neuroscience [@ica1; @ica2], bioinformatics [@bio], and causal discovery [@causal]. ICA has been applied to matrix data but there is a growing demand to apply ICA to more heterogeneous data such as multiple matrices and tensors (high-dimensional arrays), which are higher-order data structures than matrices [@mica; @groupica1; @groupica2; @multilinearica]. To meet these requirements, I originally developed \texttt{iTensor}, which is an R/CRAN package to perform some ICA-based matrix/tensor decomposition algorithms (\url{https://cran.r-project.org/web/packages/iTensor/index.html}).

# Statement of need

Currently, the most comprehensive implementation for ICA-related algorithms is the Group ICA Of fMRI Toolbox (GIFT, http://mialab.mrn.org/software/gift), but it is not freely available because it is implemented in MATLAB. Also, some open-source software is implemented in R and Python but those only focus on fewer algorithms. To fill this gap, I originally implemented some ICA-based matrix/tensor decomposition algorithms in R.

\texttt{iTensor} provides the ICA-based matrix/tensor decomposition functions as follows:

- \texttt{ICA}: ICA (3 classic models including InfoMax [@infomax1; @infomax2], ExtInfoMax [@extinfomax], and FastICA [@fastica])
- \texttt{ICA2}: ICA (9 modern models including JADE [@jade], AuxICA1/2 [@auxica], SIMBEC [@simbec], AMUSE [@amuse], SOBI [@sobi], FOBI [@fobi], ProDenICA [@prodenica], and RICA [@rica])
- \texttt{MICA}: Multimodal ICA [@mica]
- \texttt{GroupICA}: Group ICA [@groupica1; @groupica2]
- \texttt{MultilinearICA}: Multilinear ICA [@multilinearica]

I also implemented CorrIndex [@corrindex], which is a performance index to evaluate ICA results.

# Example

ICA and plots in Figure \autoref{fig:ica} can be easily reproduced on any machine where R is pre-installed by using the following commands in R:

```r
# Install package required (one per computer)
install.packages("iTensor")

# Load required package (once per R instance)
library("iTensor")

# Load Toy data
data1 <- iTensor::toyModel("ICA_Type1")

# Perform ICA
set.seed(1234)
out.JADE <- ICA2(X=data1$X_observed, J=3, algorithm="JADE")

# Source Signal extracted by ICA (If it becomes an upright square,
# the calculation is successful)
pairs(data1$X_observed)
pairs(Re(out.JADE$S))

# CorrIndex (0.2211509, the closer to 0, the better the performance)
CorrIndex(cor(data1$S, Re(out.JADE$S)))
```

![ICA with time-independent sub-gaussian data\label{fig:ica}](figure.png){ width=100% }

# Related work

There are some packages to perform ICA for matrix, matrices, and tensor but such packages focus on only a few algorithms. \texttt{iTensor} is the most comprehensive and unified package to perform ICA-based matrix/tensor decomposition as follows.

| Name (function or package) | Language | ICA for matrix | ICA for matrices | ICA for tensor | Reference |
|:------ | :---- | ----: | ----: | ----: | :----: |
| \texttt{scikit-learn} | Python | 1 | - | - | @sklearn |
| \texttt{MNE} | Python | 1 | - | - | @mne |
| \texttt{rica} | MATLAB |1 | - | - | @rica |
| \texttt{fastICA} | R | 1 | - | - | @fastica |
| \texttt{fICA} | R | 1 | - | - | @fastica |
| \texttt{JADE} | R | 1 | - | - | @jade |
| \texttt{ProDenICA} | R | 1 | - | - | @prodenica |
| \texttt{ica} | R | 3 | - | - | @ica1; @ica2 |
| \texttt{groupICA} | R | - | 1 | - | @groupica2 |
| \texttt{coroICA} | R/Python/MATLAB | - | 2 | - | @coroica |
| \texttt{BrainVoyager} | MATLAB | 1 | - | - | @brainvoyager1; @brainvoyager2 |
| \texttt{FMRLAB} | MATLAB | 1 | - | - | @fmrlab |
| \texttt{GIFT} | MATLAB | 14 | 1 | - | @gift |
| \texttt{tensorBSS} | R | - | - | 6 | @tensorbss |
| \texttt{iTensor} | R | 12 | 2 | 1 | This paper |

: Existing ICA-related packages

For MICA [@mica] and Multilinear ICA [@multilinearica], there is no package without \texttt{iTensor} to perform them.

# References
