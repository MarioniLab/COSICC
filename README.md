The COSICC R package can be installed using the devtools R package as follows.

if(!require(devtools)){
    install.packages("devtools")
}

if(!require(BiocManager)){
    install.packages("BiocManager")
}

options(repos = BiocManager::repositories())

BiocManager::install("scRNAseq")

devtools::install_github("https://github.com/theislab/destiny",dependencies = TRUE)

devtools::install_github("https://github.com/MarioniLab/COSICC",dependencies = TRUE, upgrade = TRUE,build_vignettes = TRUE)



There is a pfd vignette calledd Intro_to_COSICC_pdf.pdf.

Once the R package has been installed, a vignette can also be accessed by typing vignette("Intro_to_COSICC",package="COSICC").
