# BIRSBIO 2020 scNMTseq challenge using a PLS-based approach

pkdown report available at: https://ajabadi.github.io/BIRSBIO2020.scNMTseq.PLS

## How to take your vignette into a `pkdown` package deployed using GitHub Actions

:zero: Choose a name for your package. Recommended convention: `{hackathon_event}.{theme}.{method/topic}`. e.g. `BIRSBIO2020.scNMTseq.Benchmarking`, `BIRSBIO2020.scProteomics.LatentDirichlet`. Ensure that you also have a [Docker](https://hub.docker.com/) account for automatic containerization as well.

:one: Use https://github.com/seandavi/BuildABiocWorkshop2020 as template for your analysis package. You can simply click [here](https://github.com/seandavi/BuildABiocWorkshop2020/generate) to accomplish this. Make it public and include all branches (to keep `master` & `gh-pages`, you can delete the rest).

:two: Follow the steps outlined [here](https://github.com/seandavi/BuildABiocWorkshop2020/blob/master/vignettes/HOWTO_BUILD_WORKSHOP.Rmd) to set up your own workflow and create a package from your analyses. The notebooks should go in `./vignettes` folder and source files in `./R` (or simply include them in notebooks). Include all dependencies in `DESCRIPTION` and ensure it can be installed, Ensure `devtools::build_vignettes()` can successfully create the vignettes locally before testing using GitHub Actions.

   :bulb: If you use R and have python dependencies, [this setup](https://github.com/ajabadi/BIRSBIO2020.scNMTseq.PLS/commit/3155bbab63129e3734e155f9f245c3a386230627#diff-5822d7d51c0024ec80488aa8a41ba9caR5-R20) should help as an example.
      
   :bulb: The docker image name should only include lowercase letters, integers, `_` and `-`.
   
   :bulb: Follow [this article](https://docs.github.com/en/actions/configuring-and-managing-workflows/creating-and-storing-encrypted-secrets) to add Docker credentials to GitHub Secrets.
   
:three: Push and ensure the workflow deploys successfully using GitHub Actions (https://github.com/YOUR_GITHUB_USER/REPO/actions)

:four: Please update these steps through a PR if required
