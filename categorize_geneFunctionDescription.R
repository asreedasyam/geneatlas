
source("geneFunctionDescription_utils.R")

## -- download annotation files for Gene Atlas species from Phytozome (https://phytozome-next.jgi.doe.gov/) to the wrk.dir
wrk.dir = "/Users/ga.wrk/Downloads"
ga.info = "GA.v2"
annotation.files <- unlist(strsplit("Athaliana_167_TAIR10.annotation_info.txt,Athaliana_447_Araport11.annotation_info.txt,Bdistachyon_314_v3.1.annotation_info.txt,Creinhardtii_281_v5.6.annotation_info.txt,Egrandis_297_v2.0.annotation_info.txt,Gmax_508_Wm82.a4.v1.annotation_info.txt,Kfedtschenkoi_382_v1.1.annotation_info.txt,Mtruncatula_285_Mt4.0v1.annotation_info.txt,Phallii_495_v3.1.annotation_info.txt,PhalliiHAL_496_v2.1.annotation_info.txt,Ppatens_318_v3.3.annotation_info.txt,Ptrichocarpa_533_v4.1.annotation_info.txt,Pvirgatum_516_v5.1.annotation_info.txt,Sbicolor_454_v3.1.1.annotation_info.txt,SbicolorRio_468_v2.1.annotation_info.txt,Sfallax_522_v1.1.annotation_info.txt,Sitalica_312_v2.2.annotation_info.txt,Sviridis_500_v2.1.annotation_info.txt,Lalbus_567_v1.annotation_info.txt",  split = ","))
names(annotation.files) <- unlist(strsplit("Athaliana.TAIR10,Athaliana.Araport11,Bdistachyon.v3.1,Creinhardtii.v5.6,Egrandis.v2.0,Gmax.Wm82.a4.v1,Kfedtschenkoi.v1.1,Mtruncatula.Mt4.0v1,Phallii.v3.1,PhalliiHAL.v2.1,Ppatens.v3.3,Ptrichocarpa.v4.1,Pvirgatum.v5.1,Sbicolor.v3.1.1,SbicolorRio.v2.1,Sfallax.v1.1,Sitalica.v2.2,Sviridis.v2.1,Lalbus.v1.1", split = ","))


## -- generate lexicon for gene function descriptions
generate.lexicon_for_geneDescriptions(annotation.files, 
                                      description.column = "Description",
                                      file.suffix = ga.info, 
                                      out.dir = wrk.dir)

polarity.file <- list.files(wrk.dir, glob2rx("functionDescriptions_sentiment_PolarityDT*.csv"), full.names = T)

sentiment_threshold = 0.1999

## -- Categorize gene function descriptions of each species based on the sentiment score
for(i in seq(annotation.files)){
  message(sprintf("*********** %s ***********", names(annotation.files)[i]))
  d <- categorize.functionDescriptions(annotation.file = annotation.files[i], 
                                       polarity.file = polarity.file, 
                                       file.suffix = names(annotation.files)[i],
                                       description.column = "Description",
                                       sentiment_threshold = sentiment_threshold,
                                       out.dir = wrk.dir,
                                       verbose = TRUE)
}

fd.files <- list.files(wrk.dir, glob2rx("FunctionDescription_categories_*.csv"), full = T)
names(fd.files) <- gsub("FunctionDescription_categories_|.csv", "", basename(fd.files))
fdList <- lapply(fd.files, function(x) {d <- fread(x); d$V1 = NULL; return(d)})
save(fdList, file = sprintf("%s/%s.species_FunctionDescriptionCategories.RData", out.dir, ga.info))

