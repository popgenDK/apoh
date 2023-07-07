
# apoh - Admixture Pedigrees Of Hybrids

apoh is a software to infer, explore, rank and visualize recent admixture pedigrees of hybrids. This shiny provides a Graphical User Interface (GUI) to do so interactively.

More details on apoh can be found in [the github page](https://github.com/popgenDK/apoh)

## GUI usage

### Input

The left panel allows you to enter the input:

- The main and only input apoh needs are the estimates of paired ancestry proportions that can be produced by NGSremix -bothanc 1 option (see [github to download and run NGSremix](https://github.com/popgenDK/apoh)). TheNGSremix output file needs to be uploaded as the model estiamtes in order to run apoh.

- Additionally you can upload a text file with indivdiual IDs you want apoh to use for each sample in the output. This must be a text file with a single column and one row per sample. If no file is provided, apoh will use the first columnt of the model estiamtes file as sample ID.

- Most of the information apoh reports are specific to a single sample, and thus require to select an indidivual to explore in the left panel.

- Finally, apoh deals with uncertainty by allowing to explore and compare multiple pedigrees that can lead to the observed estiamtes. In the left tab you can chose how many pedigrees you which to explore.

### Tabs

There are 6 tabs with different apoh outputs. The first 2 and the last one report information on all samples uploaded, while the other 3 middle ones to a single sample 

- Parental admixture plot: plot of estiamted parental admixutre proportions for all samples. Allow to visualize the parental admixutre proporitons estimated by NGSremix that apoh uses as input.

- Summary indices table: table with summary indices for all samples.

- Individual ordered paired ancestry plots: Plots of the ordered paired ancestry proportions for an individual. This is the main information on which apoh recent admixture pedigrees inference is based. It shows from left to right:
  - Model estimates: the estimates for that sample the parental admixture model by the NGSremix.
  - Expected under independent pedigree: the expectation under a model with independent ancestries, that corresponds to a case where the paired ancestries give no evidence to favor recent admixture.
  - Expected under pedigree x: expected under a certain recent admixutre pedigree; will show as many as selected top compatible pedigrees in the left tab.

- Individual pedigrees: visualization of the inferred pedigrees. It shows the independent ancestries pedigree, and the top compatible recent admixutre pedigrees selected.

- Individual unordered paired ancestries plots: Plots of the unordered paired ancestry proprotions. While less informative than their ordered version, the comparison betweent the two models of inference allows to evaluate if the model assumptions are met, and thus if the inference is reliable.

- Pedigrees bootstrap support table: for each pedigree (the independent pedigree and the number of recent admixutre pedigrees selected) shows the proporiton of bootstrap while each sample has that pedigree as the closest one to the estimates. Requires that the input file with model estimates has the bootstrap replicates (done with the flag -boot 1 when running NGSremix. Currently only works with beagle file input)