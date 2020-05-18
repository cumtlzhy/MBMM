# The MBMM user's guider
Multivariate Beta Mixture Model-Based Clustering Algorithm
## Function description
MBMM is used to find potential m<sup>6</sup>A co-methylation modules based on MeRIP-Seq data, where sites are hyper-methylated or hypo-methylated simultaneously across the specific conditions.
## Data discription
The data fed into MBMM should be **a numerical matrice (or dataframe), of which all the values should be in the range of 0 and 1**. As you known,  MeRIP-seq data profiles the m<sup>6</sup>A epi-transcriptome by *input* and *IP* samples, thus, you should preprocess the MeRIP-seq data which you obtained before using the MBMM.
## Demo
You can run the MBMM based on the simulated data we have supplied, then MBMM will give you the correct model and the related parameters.
```
   Data_lables_sample <- read.csv("dataset.csv", header = TRUE)
   Data <- Data_lables_sample[,1:4]
   model <- model_select(cmax = 6,x = Data,eps = 1e-5)
```
## Contact
Please contact the maintainer of MBMM if you have encountered any problems:

Zhaoyang Liu: <zy.liu@cumt.edu.cn>
