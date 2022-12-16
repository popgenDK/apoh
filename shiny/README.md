
# apoh - Admixture Pedigrees Of Hybrids

apoh is a software to infer, explore, rank and visualize recent admixture pedigrees of hybrids. As input it takes the estimates of paired ancestry proportions that can be produced by NGSremix -bothanc 1 option. This shiny provides a Graphical User Interface (GUI) to do so interactively.

## Usage

To upload the data 

Shiny to process paired ancestries and parental admixture proporitons and visualize most likely pedigrees to generate observed proporitons.

Input format (for N individuals, K ancestral populations):

- Parental admixture: file with N rows and 2*K columns, each row estimates for a sample, first K columns are admixture proporitons of parent 1 and K+1 to 2*K columns are admixture proprotions of parent 2.

- Paired ancestries (optional): file with N rows and (K*(K-1))/2 + K columns, with unordered ancestry proporitons estiamtes for a sample.

- Individual ID (optional): file with N rows and 1 column, with ID for each of the samples.

