#' CTCF Demo data
#'
#' A transformed demo data set used for showing the workflow of BIMTR.
#'
#'
#' Xu B, Wang H, Wright S, Hyle J, Zhang Y, Shao Y, Niu M, Fan Y, Rosikiewicz W, Djekidel MN, Peng J, Lu R, Li C.
#' Acute depletion of CTCF rewires genome-wide chromatin accessibility. Genome Biol. 2021 Aug 24;22(1):244.
#' doi: 10.1186/s13059-021-02466-0. PMID: 34429148; PMCID: PMC8386078.
#'
#' @format ## 'CTCF_Demo'
#' A data frame with 13511 rows and 4 columns
#' \describe{
#'   \item{TF}{transcription factor HGNC symbol.}
#'   \item{GOOD}{number of good cases of demo input ATAC-seq peakset.}
#'   \item{TOTAL}{number of total informative cases of demo input ATAC-seq peakset.}
#'   \item{Ratio}{ratio of good cases over all informative cases.}
#' }
#' @source https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE153237
"CTCF_Demo"


#' ZBTB7A Demo data
#'
#' A transformed demo data set used for showing the workflow of BIMTR.
#'
#'
#' King AJ, Songdej D, Downes DJ, Beagrie RA, Liu S, Buckley M, Hua P, Suciu MC, Marieke Oudelaar A, Hanssen LLP, Jeziorska D,
#' Roberts N, Carpenter SJ, Francis H, Telenius J, Olijnik AA, Sharpe JA, Sloane-Stanley J, Eglinton J, Kassouf MT, Orkin SH,
#' Pennacchio LA, Davies JOJ, Hughes JR, Higgs DR, Babbs C.
#' Reactivation of a developmentally silenced embryonic globin gene. Nat Commun. 2021 Jul 21;12(1):4439.
#' doi: 10.1038/s41467-021-24402-3. PMID: 34290235; PMCID: PMC8295333.
#'
#' @format ## 'ZBTB7A_Demo'
#' A data frame with 13511 rows and 4 columns
#' \describe{
#'   \item{TF}{transcription factor HGNC symbol.}
#'   \item{GOOD}{number of good cases of demo input ATAC-seq peakset.}
#'   \item{TOTAL}{number of total informative cases of demo input ATAC-seq peakset.}
#'   \item{Ratio}{ratio of good cases over all informative cases.}
#' }
#' @source https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE173416
"ZBTB7A_Demo"

#' CREs_list
#'
#' A pre-processed list that contains indices of all cis-regulatory elements with 3 type of bin width.
#'
#'
#' ENCODE Project Consortium, Jill E. Moore, Michael J. Purcaro, Henry E. Pratt, Charles B. Epstein, Noam Shoresh, Jessika Adrian, et al. 2020.
#' Expanded Encyclopaedias of DNA Elements in the Human and Mouse Genomes. Nature 583 (7818): 699–710.
#' doi: 10.1038/s41586-020-2493-4. Epub 2020 Jul 29. Erratum in: Nature. 2022 May;605(7909):E3. PMID: 32728249; PMCID: PMC7410828.
#'
#' @format ## 'CREs_list'
#' A list with 3 vectors.
#' \describe{
#'   \item{CREs_1000}{cis-regulatory elements transformed to indices with 1000 bin width}
#'   \item{CREs_500}{cis-regulatory elements transformed to indices with 500 bin width}
#'   \item{CREs_100}{cis-regulatory elements transformed to indices with 100 bin width}
#' }
#' @source https://www.nature.com/articles/s41586-020-2493-4
"CREs_list"

#' PLS_list
#'
#' A pre-processed list that contains indices of all candidate promoters with 3 type of bin width.
#'
#'
#' ENCODE Project Consortium, Jill E. Moore, Michael J. Purcaro, Henry E. Pratt, Charles B. Epstein, Noam Shoresh, Jessika Adrian, et al. 2020.
#' Expanded Encyclopaedias of DNA Elements in the Human and Mouse Genomes. Nature 583 (7818): 699–710.
#' doi: 10.1038/s41586-020-2493-4. Epub 2020 Jul 29. Erratum in: Nature. 2022 May;605(7909):E3. PMID: 32728249; PMCID: PMC7410828.
#'
#' @format ## 'PLS_list'
#' A list with 3 vectors.
#' \describe{
#'   \item{PLS_1000}{candidate promoters transformed to indices with 1000 bin width}
#'   \item{PLS_500}{candidate promoters transformed to indices with 500 bin width}
#'   \item{PLS_100}{candidate promoters transformed to indices with 100 bin width}
#' }
#' @source https://www.nature.com/articles/s41586-020-2493-4
"PLS_list"

#' ELS_list
#'
#' A pre-processed list that contains indices of all candidate enhancers with 3 type of bin width.
#'
#'
#' ENCODE Project Consortium, Jill E. Moore, Michael J. Purcaro, Henry E. Pratt, Charles B. Epstein, Noam Shoresh, Jessika Adrian, et al. 2020.
#' Expanded Encyclopaedias of DNA Elements in the Human and Mouse Genomes. Nature 583 (7818): 699–710.
#' doi: 10.1038/s41586-020-2493-4. Epub 2020 Jul 29. Erratum in: Nature. 2022 May;605(7909):E3. PMID: 32728249; PMCID: PMC7410828.
#'
#' @format ## 'ELS_list'
#' A list with 3 vectors.
#' \describe{
#'   \item{ELS_1000}{candidate enhancers transformed to indices with 1000 bin width}
#'   \item{ELS_500}{candidate enhancers transformed to indices with 500 bin width}
#'   \item{ELS_100}{candidate enhancers transformed to indices with 100 bin width}
#' }
#' @source https://www.nature.com/articles/s41586-020-2493-4
"ELS_list"

