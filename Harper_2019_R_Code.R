
# Code I wrote to conduct the analyses in Harper, Malone, and Iacono (2019; https://doi.org/10.1017/S0033291719003258). Includes data wrangling, data inspection, linear mixed modeling and generalized mixed effects regressions for longitudinal prediction, and data visualization 

library(plyr);library(reshape2);library(ggplot2);library(effects);library(lmerTest); library(pbkrtest); library(lme4);library(stats);library(car);

setwd('/LOCAL/');

# Load data and create variables from unstructured datasets

TD.data = read.table("./output_data/ES14_rho_stim_s25lLCD_erp-win-rs256-bsln_450250_p3.dat", header=TRUE, sep="\t", na="-999");
ROI.data = read.table("./output_data/ES14_rho_stim_s25lLCD_TP_wvltlog130-wintfd-rs256-t64-f25-fqAD-cnspc_450250_dltht2al.dat", header=TRUE, sep="\t", na="-999");

# Sort by individual
TD.data = TD.data[order(TD.data$subname),];
ROI.data = ROI.data[order(ROI.data$subname),];

    # ROI Total Power
TPROIthmFCZT = ROI.data$t2m[ROI.data$elecname=='FCZ'&ROI.data$catcodes=='TRUE'];
TPROIthmFCZF = ROI.data$t2m[ROI.data$elecname=='FCZ'&ROI.data$catcodes=='FALSE'];

TPROIthmCZT = ROI.data$t2m[ROI.data$elecname=='CZ '&ROI.data$catcodes=='TRUE'];
TPROIthmCZF = ROI.data$t2m[ROI.data$elecname=='CZ '&ROI.data$catcodes=='FALSE'];

TPROIthmFZT = ROI.data$t2m[ROI.data$elecname=='FZ '&ROI.data$catcodes=='TRUE'];
TPROIthmFZF = ROI.data$t2m[ROI.data$elecname=='FZ '&ROI.data$catcodes=='FALSE'];

TPROIthmFC1T = ROI.data$t2m[ROI.data$elecname=='FC1'&ROI.data$catcodes=='TRUE'];
TPROIthmFC1F = ROI.data$t2m[ROI.data$elecname=='FC1'&ROI.data$catcodes=='FALSE'];
TPROIthmFC2T = ROI.data$t2m[ROI.data$elecname=='FC2'&ROI.data$catcodes=='TRUE'];
TPROIthmFC2F = ROI.data$t2m[ROI.data$elecname=='FC2'&ROI.data$catcodes=='FALSE'];

TPROIthmFCZ12T = rowMeans(cbind(TPROIthmFCZT, TPROIthmFC1T, TPROIthmFC2T));
TPROIthmFCZ12F = rowMeans(cbind(TPROIthmFCZF, TPROIthmFC1F, TPROIthmFC2F));

TPROIthmFCZCZT = rowMeans(cbind(TPROIthmFCZT, TPROIthmCZT));
TPROIthmFCZCZF = rowMeans(cbind(TPROIthmFCZF, TPROIthmCZF));

TPROIthmFZFCZT = rowMeans(cbind(TPROIthmFCZT, TPROIthmFZT));
TPROIthmFZFCZF = rowMeans(cbind(TPROIthmFCZF, TPROIthmFZF));

TPROIthmFZFCZCZT = rowMeans(cbind(TPROIthmFZT, TPROIthmFCZT, TPROIthmCZT));
TPROIthmFZFCZCZF = rowMeans(cbind(TPROIthmFZF, TPROIthmFCZF, TPROIthmCZF));

    # TD

TDp3pPOZT = TD.data$p3p[TD.data$elecname=='POZ'&TD.data$catcodes=='TRUE'];
TDp3pPOZF = TD.data$p3p[TD.data$elecname=='POZ'&TD.data$catcodes=='FALSE'];

TDp3mPOZT = TD.data$p3m[TD.data$elecname=='POZ'&TD.data$catcodes=='TRUE'];
TDp3mPOZF = TD.data$p3m[TD.data$elecname=='POZ'&TD.data$catcodes=='FALSE'];

TDp3pPZT = TD.data$p3p[TD.data$elecname=='PZ '&TD.data$catcodes=='TRUE'];
TDp3pPZF = TD.data$p3p[TD.data$elecname=='PZ '&TD.data$catcodes=='FALSE'];

TDp3mPZT = TD.data$p3m[TD.data$elecname=='PZ '&TD.data$catcodes=='TRUE'];
TDp3mPZF = TD.data$p3m[TD.data$elecname=='PZ '&TD.data$catcodes=='FALSE'];

TDp3pP1T = TD.data$p3p[TD.data$elecname=='P1 '&TD.data$catcodes=='TRUE'];
TDp3pP1F = TD.data$p3p[TD.data$elecname=='P1 '&TD.data$catcodes=='FALSE'];

TDp3mP1T = TD.data$p3m[TD.data$elecname=='P1 '&TD.data$catcodes=='TRUE'];
TDp3mP1F = TD.data$p3m[TD.data$elecname=='P1 '&TD.data$catcodes=='FALSE'];

TDp3pP2T = TD.data$p3p[TD.data$elecname=='P2 '&TD.data$catcodes=='TRUE'];
TDp3pP2F = TD.data$p3p[TD.data$elecname=='P2 '&TD.data$catcodes=='FALSE'];

TDp3mP2T = TD.data$p3m[TD.data$elecname=='P2 '&TD.data$catcodes=='TRUE'];
TDp3mP2F = TD.data$p3m[TD.data$elecname=='P2 '&TD.data$catcodes=='FALSE'];

TDp3pPZPOZT = rowMeans(cbind(TDp3pPZT, TDp3pPOZT));
TDp3pPZPOZF = rowMeans(cbind(TDp3pPZF, TDp3pPOZF));

TDp3mPZPOZT = rowMeans(cbind(TDp3mPZT, TDp3mPOZT));
TDp3mPZPOZF = rowMeans(cbind(TDp3mPZF, TDp3mPOZF));

TDp3pPZ12T = rowMeans(cbind(TDp3pPZT, TDp3pP1T, TDp3pP2T));
TDp3pPZ12F = rowMeans(cbind(TDp3pPZF, TDp3pP1F, TDp3pP2F));

TDp3mPZ12T = rowMeans(cbind(TDp3mPZT, TDp3mP1T, TDp3mP2T));
TDp3mPZ12F = rowMeans(cbind(TDp3mPZF, TDp3mP1F, TDp3mP2F))

# Create/recode demographic variables
sex_org    = TD.data$sex[TD.data$elecname=='FCZ'&TD.data$catcodes=='TRUE'];
sex = sex_org; sex[sex==1] = -1; sex[sex==2] = 1; # effect code [-1, 1] for males and females
idyrfam    = TD.data$idyrfam[TD.data$elecname=='FCZ'&TD.data$catcodes=='TRUE'];
ID         = TD.data$id[TD.data$elecname=='FCZ'&TD.data$catcodes=='TRUE'];
idsc       = TD.data$idsc[TD.data$elecname=='FCZ'&TD.data$catcodes=='TRUE'];
subname    = TD.data$subname[TD.data$elecname=='FCZ'&TD.data$catcodes=='TRUE'];
zyg_org    = TD.data$zyg[TD.data$elecname=='FCZ'&TD.data$catcodes=='TRUE'];
zyg = zyg_org; zyg[zyg==1] = 1; zyg[zyg==2] = -1; # effect code [-1, 1] for DZ and MZ

# Put it all together into a data.frame
dataset = data.frame(subname, ID, idyrfam, idsc, sex, zyg, TDp3pPOZT, TDp3pPOZF, TDp3mPOZT, TDp3mPOZF, TDp3pPZT, TDp3pPZF, TDp3mPZT, TDp3mPZF, TDp3pPZ12T, TDp3pPZ12F, TDp3mPZ12T, TDp3mPZ12F, TDp3pPZPOZT, TDp3pPZPOZF, TDp3mPZPOZT, TDp3mPZPOZF, TPROIthmFCZT, TPROIthmFCZF, TPROIthmFCZ12T, TPROIthmFCZ12F, TPROIthmFCZCZT, TPROIthmFCZCZF, TPROIthmFZFCZCZT, TPROIthmFZFCZCZF, TPROIthmCZT,TPROIthmCZF, TPROIthmFZFCZT, TPROIthmFZFCZF);


# Data checking/cleaning/merging 
{
library(dplyr)

# Merge behavioral data with dataset (left join)
beh.data = read.table("./output_data/ES14_beg_behdata.dat", header=TRUE, sep="\t", na=c("NaN","-999")); 
dataset = merge(dataset,beh.data,by="ID",all.x=TRUE);
rm(list=setdiff(ls(), "dataset")); gc();

# Read in and merge external data
CSUdrink.data = read.table("../ES/data/ES111417_CSU_drink_index_TOB_CAN.dat", header=TRUE, sep="\t", na=c(-999,'NA'));
CSUdrink.data = CSUdrink.data[!duplicated(CSUdrink.data$ID),]; # remove duplicated rows
dataset = merge(dataset, CSUdrink.data, by="ID",all.x=TRUE);
rm(list=setdiff(ls(), "dataset")); gc();

# Recode 98 to 0 and 96 to NAs
recode_phenos <- c('TOB_EVER_USED_11', 'TOB_USED_MOST_11', 'ALC_EVER_DRUNK_11', 'ALC_EVER_DRUNK_12m_11', 'TOB_EVER_USED_14', 'TOB_USED_MOST_14', 'ALC_EVER_DRUNK_14', 'ALC_EVER_DRUNK_12m_14', 'TOB_EVER_USED_17', 'TOB_USED_MOST_17', 'ALC_EVER_DRUNK_17', 'ALC_EVER_DRUNK_12m_17')
dataset <- dataset %>% mutate_at(vars(recode_phenos), car::recode, "98=0;96=NA")

# code TOB_EVER_USED_11 to 0 if TOB_USED_MOST_11 is 0 (indicating no use)
dataset$TOB_EVER_USED_11[!complete.cases(dataset$TOB_EVER_USED_11)&dataset$TOB_USED_MOST_11 %in% c(0)] = 0

# code TOB_EVER_USED_14 to 0 if TOB_USED_MOST_14 is 0 (indicating no use)
dataset$TOB_EVER_USED_14[!complete.cases(dataset$TOB_EVER_USED_14)&dataset$TOB_USED_MOST_14 %in% c(0)] = 0

# code TOB_EVER_USED_17 to 0 if TOB_USED_MOST_17 is 0 (indicating no use)
dataset$TOB_EVER_USED_17[!complete.cases(dataset$TOB_EVER_USED_17)&dataset$TOB_USED_MOST_17 %in% c(0)] = 0

# Recode into age in years
unique(dataset$ALC_AGE_1st_DRUNK_11)
dataset$ALC_AGE_1st_DRUNK_11 = recode(dataset$ALC_AGE_1st_DRUNK_11,"1=8;2=9;3=10;4=11;5=12;6=13;7=14;8=15;9=16;  98=0;96=NA",as.numeric.result=FALSE)
unique(dataset$ALC_AGE_1st_DRUNK_14)
dataset$ALC_AGE_1st_DRUNK_14 = recode(dataset$ALC_AGE_1st_DRUNK_14,"1=8;2=9;3=10;4=11;5=12;6=13;7=14;8=15;9=16;98=0;96=NA",as.numeric.result=FALSE)
unique(dataset$ALC_AGE_1st_DRUNK_17)
dataset$ALC_AGE_1st_DRUNK_17 = recode(dataset$ALC_AGE_1st_DRUNK_17,"1=11;2=12;3=13;4=14;5=15;6=16;7=17;8=18;98=0;96=NA",as.numeric.result=FALSE)

# this ID is missing for Ever Drunk at 11, but reported No to ever drunk in past 12m, has a 0 drink score, and has a age of first intox as never (NA). Recode to 0
dataset$ALC_EVER_DRUNK_11[dataset$ID %in% c(11111111)] = 0

# this ID is missing for Ever Drunk at 17, but reported No to ever drunk in past 12m, has 0 drink score, and age of first intox as never (NA). recode to 0
dataset$ALC_EVER_DRUNK_17[dataset$ID %in% 22222222] = 0


# Create aggregate tobacco measure
dataset$TOB_FREQ_12m_r_11 = apply(data.frame(dataset$TOB_FREQ_12m_CIG_11_r,dataset$TOB_FREQ_12m_DIP_11_r,dataset$TOB_FREQ_12m_PIP_11_r),1, function(x) max(x))
dataset$TOB_AMT_12m_r_11 = apply(data.frame(dataset$TOB_AMT_12m_CIG_11_r,dataset$TOB_AMT_12m_DIP_11_r,dataset$TOB_AMT_12m_PIP_11_r),1, function(x) max(x))
dataset$TOB_FREQ_12m_r_14 = apply(data.frame(dataset$TOB_FREQ_12m_CIG_14_r,dataset$TOB_FREQ_12m_DIP_14_r,dataset$TOB_FREQ_12m_PIP_14_r),1, function(x) max(x))
dataset$TOB_AMT_12m_r_14 = apply(data.frame(dataset$TOB_AMT_12m_CIG_14_r,dataset$TOB_AMT_12m_DIP_14_r,dataset$TOB_AMT_12m_PIP_14_r),1, function(x) max(x))

# Load, process, and merge external data
drink.data = read.table("../ES/data/ESx_longitudinal_drink_index.dat", header=TRUE, sep="\t", na=(NA));
drink.data$ID = drink.data$id; drink.data$id = NULL
drink.data = drink.data[!duplicated(drink.data$ID),];
drink.data = subset(drink.data, select = c(ID, CSU_drink_index_IN, CSU_drink_index_FU1,  CSU_drink_index_FU2, SAM_drink_index_FU2, FREQ_SAM_FU2, AMT_SAM_FU2, MAX_SAM_FU2, INTOX_SAM_FU2, FREQ_CSU_FU2, AMT_CSU_FU2, MAX_CSU_FU2, INTOX_CSU_FU2))
dataset = merge(dataset, drink.data, by="ID",all.x=TRUE);
rm(list=setdiff(ls(), "dataset"));

binge.data = read.table("../ES17/data/ES17_SAM_BINGE.dat", header=TRUE, sep="\t", na=(NA));
binge.data = binge.data[!duplicated(binge.data $ID),];
dataset = merge(dataset, binge.data, by="ID",all.x=TRUE);
rm(list=setdiff(ls(), "dataset"));

# Recode variables according to handbook
unique(dataset$ALC_EVER_USED_SAM_17)
dataset$ALC_EVER_USED_SAM_17 = recode(dataset$ALC_EVER_USED_SAM_17,"1=0;2=0;5=1",as.numeric.result=FALSE)

unique(dataset$EVER_BINGE)
dataset$EVER_BINGE = recode(dataset$EVER_BINGE,"1=0;5=1",as.numeric.result=FALSE)

# code non-users binge from appropriately skipped to 0
dataset$EVER_BINGE[dataset$ALC_EVER_USED_SAM_17 %in% c(0)] = 0  
# code binge appropriately skipped to 0 
dataset$EVER_BINGE[dataset$ALC_EVER_USED_SAM_17 %in% c(1) & dataset$EVER_BINGE %in% c(98)] = 0 
dataset = rename(dataset,c("EVER_BINGE" = "EVER_BINGE_SAM_17"));

# Calculate alternate binge drinking item from the CSU max drinks item
dataset$CSU_BINGE_14 = -999
dataset$CSU_BINGE_14[(dataset$ALC_MAXDRINK_12m_14 %in% 1:6  & !dataset$ALC_MAXDRINK_12m_14 %in% 98) & dataset$sex %in% -1] = 1  # males ≥ 5 drinks
dataset$CSU_BINGE_14[(dataset$ALC_MAXDRINK_12m_14 %in% 1:7  & !dataset$ALC_MAXDRINK_12m_14 %in% 98) & dataset$sex %in%  1] = 1  # females ≥ 4 drinks
dataset$CSU_BINGE_14[(dataset$ALC_MAXDRINK_12m_14 %in% 7:11 |  dataset$ALC_MAXDRINK_12m_14 %in% 98 | dataset$ALC_MAXDRINK_12m_14 %in% 0) & dataset$sex %in% -1] = 0
dataset$CSU_BINGE_14[(dataset$ALC_MAXDRINK_12m_14 %in% 8:11 |  dataset$ALC_MAXDRINK_12m_14 %in% 98 | dataset$ALC_MAXDRINK_12m_14 %in% 0) & dataset$sex %in%  1] = 0

dataset$CSU_BINGE_17 = -999
dataset$CSU_BINGE_17[(dataset$ALC_MAXDRINK_12m_17 %in% 1:6  & !dataset$ALC_MAXDRINK_12m_17 %in% 98) & dataset$sex %in% -1] = 1  # males ≥ 5 drinks
dataset$CSU_BINGE_17[(dataset$ALC_MAXDRINK_12m_17 %in% 1:7  & !dataset$ALC_MAXDRINK_12m_17 %in% 98) & dataset$sex %in%  1] = 1  # females ≥ 4 drinks
dataset$CSU_BINGE_17[(dataset$ALC_MAXDRINK_12m_17 %in% 7:11 |  dataset$ALC_MAXDRINK_12m_17 %in% 98 | dataset$ALC_MAXDRINK_12m_17 %in% 0) & dataset$sex %in% -1] = 0
dataset$CSU_BINGE_17[(dataset$ALC_MAXDRINK_12m_17 %in% 8:11 |  dataset$ALC_MAXDRINK_12m_17 %in% 98 | dataset$ALC_MAXDRINK_12m_17 %in% 0) & dataset$sex %in%  1] = 0

# 9 IDs did not binge at 17 but did at 14, recode these instances to 1
table(dataset$CSU_BINGE_14,dataset$CSU_BINGE_17) 
dataset$CSU_BINGE_17[(dataset$CSU_BINGE_14 %in% c(1)&dataset$CSU_BINGE_17 %in% c(0))==TRUE] = 1

# Read in and merge demographic data
agesx = read.table("../ES23/gonogo/data/demodata_021518.txt", header=TRUE, sep=",", na="NA");
agesx = subset(agesx, select = c("ID","IDES","AGE_IN","AGE_FU1","AGE_FU2"));
dataset = merge(dataset,agesx,by="ID",all.x=TRUE);
rm(list=setdiff(ls(), "dataset"));

dataset$AGE_IN_scaled = scale(dataset$AGE_IN, scale = FALSE);
dataset$AGE_FU1_scaled = scale(dataset$AGE_FU1, scale = FALSE);
dataset$AGE_FU2_scaled = scale(dataset$AGE_FU2, scale = FALSE);

dataset$Gender = NULL;
dataset$Gender[dataset$sex ==  1] = 'F';
dataset$Gender[dataset$sex == -1] = 'M';
dataset$Gender = as.factor(dataset$Gender);

# calculate weighted effect coding across total sample
dataset$sex_wei[dataset$sex ==  1] =  1;
dataset$sex_wei[dataset$sex == -1] = -1*(sum(dataset$Gender=='F')/sum(dataset$Gender=='M'));

dataset$sex = as.factor(dataset$sex)

# Remove cases without complete CSU drinking data at 11 14 and 17
dataset = dataset[complete.cases(dataset$CSU_drink_index_11)&complete.cases(dataset$CSU_drink_index_14)&complete.cases(dataset$CSU_drink_index_17),];

# Calculate Age 11/14 composite measures
dataset$CSU_drink_index_1114 = integer(dim(dataset)[1])
dataset$CSU_drink_index_1114 = apply(data.frame(dataset$CSU_drink_index_11 ,dataset$CSU_drink_index_14),1, function(x) sum(x))

dataset$ALC_FREQ_12m_bin_1114 = integer(dim(dataset)[1])
dataset$ALC_FREQ_12m_bin_1114 = apply(data.frame(dataset$ALC_FREQ_12m_11_bin ,dataset$ALC_FREQ_12m_14_bin),1, function(x) sum(x))

dataset$ALC_AMT_12m_bin_1114 = integer(dim(dataset)[1])
dataset$ALC_AMT_12m_bin_1114 = apply(data.frame(dataset$ALC_AMT_12m_11_bin ,dataset$ALC_AMT_12m_14_bin),1, function(x) sum(x))

dataset$ALC_MAXDRINK_12m_bin_1114 = integer(dim(dataset)[1])
dataset$ALC_MAXDRINK_12m_bin_1114 = apply(data.frame(dataset$ALC_MAXDRINK_12m_11_bin ,dataset$ALC_MAXDRINK_12m_14_bin),1, function(x) sum(x))

dataset$ALC_INTOX_12m_bin_1114 = integer(dim(dataset)[1])
dataset$ALC_INTOX_12m_bin_1114 = apply(data.frame(dataset$ALC_INTOX_12m_11_bin ,dataset$ALC_INTOX_12m_14_bin),1, function(x) sum(x))

dataset$ALC_EVER_USED_NOPERM_14[(dataset$ALC_EVER_USED_NOPERM_1114 %in% c(1)&dataset$ALC_EVER_USED_NOPERM_14 %in% c(0))==TRUE] = 1

dataset$ALC_EVER_USED_NOPERM_1114[(dataset$ALC_EVER_USED_NOPERM_11 %in% c(0)&dataset$ALC_EVER_USED_NOPERM_14 %in% c(0))==TRUE] = 0
dataset$ALC_EVER_USED_NOPERM_1114[(dataset$ALC_EVER_USED_NOPERM_11 %in% c(0)&dataset$ALC_EVER_USED_NOPERM_14 %in% c(0))==FALSE] = 1

# out of the total 593, 17 reported + alc use by age-14 but - at 17; recode those to 1 at age-17
table(dataset$ALC_EVER_USED_NOPERM_1114,dataset$ALC_EVER_USED_NOPERM_17) 
dataset$ALC_EVER_USED_NOPERM_17[(dataset$ALC_EVER_USED_NOPERM_1114 %in% c(1)&dataset$ALC_EVER_USED_NOPERM_17 %in% c(0))==TRUE] = 1
table(dataset$ALC_EVER_USED_NOPERM_1114,dataset$ALC_EVER_USED_NOPERM_17)

dataset$ALC_EVER_DRUNK_1114[(dataset$ALC_EVER_DRUNK_11 %in% c(0)&dataset$ALC_EVER_DRUNK_14 %in% c(0))==TRUE] = 0
dataset$ALC_EVER_DRUNK_1114[(dataset$ALC_EVER_DRUNK_11 %in% c(0)&dataset$ALC_EVER_DRUNK_14 %in% c(0))==FALSE] = 1

# out of the total 593, 9 reported + alc use by age-14 but - at 17; recode those to 1 at age-17
table(dataset$ALC_EVER_DRUNK_1114,dataset$ALC_EVER_DRUNK_17) 
dataset$ALC_EVER_DRUNK_17[(dataset$ALC_EVER_DRUNK_1114 %in% c(1)&dataset$ALC_EVER_DRUNK_17 %in% c(0))==TRUE] = 1
table(dataset$ALC_EVER_DRUNK_1114,dataset$ALC_EVER_DRUNK_17)

dataset$TOB_AMT_12m_CIG_r_1114 = integer(dim(dataset)[1])
dataset$TOB_AMT_12m_CIG_r_1114 = apply(data.frame(dataset$TOB_AMT_12m_CIG_11_r ,dataset$TOB_AMT_12m_CIG_14_r),1, function(x) sum(x))

dataset$TOB_FREQ_12m_CIG_r_1114 = integer(dim(dataset)[1])
dataset$TOB_FREQ_12m_CIG_r_1114 = apply(data.frame(dataset$TOB_FREQ_12m_CIG_11_r ,dataset$TOB_FREQ_12m_CIG_14_r),1, function(x) sum(x))

dataset$TOB_FREQAMT_12m_CIG_r_1114 = integer(dim(dataset)[1])
dataset$TOB_FREQAMT_12m_CIG_r_1114 = apply(data.frame(dataset$TOB_AMT_12m_CIG_r_1114 ,dataset$TOB_FREQ_12m_CIG_r_1114),1, function(x) sum(x))

dataset$TOB_FREQAMT_12m_CIG_r_17 = integer(dim(dataset)[1])
dataset$TOB_FREQAMT_12m_CIG_r_17 = apply(data.frame(dataset$TOB_AMT_12m_CIG_17_r ,dataset$TOB_FREQ_12m_CIG_17_r),1, function(x) sum(x))

dataset$TOB_FREQ_12m_r_1114 = integer(dim(dataset)[1])
dataset$TOB_FREQ_12m_r_1114 = apply(data.frame(dataset$TOB_FREQ_12m_r_11 ,dataset$TOB_FREQ_12m_r_14),1, function(x) sum(x))

dataset$CAN_FREQ_12m_r_1114 = integer(dim(dataset)[1])
dataset$CAN_FREQ_12m_r_1114 = apply(data.frame(dataset$CAN_FREQ_12m_11_r ,dataset$CAN_FREQ_12m_14_r),1, function(x) sum(x))


# Code the tobacco measure to reflect ever trying more than a small amount
dataset$TOB_EVER_USEDFULL_11 = NA
dataset$TOB_EVER_USEDFULL_11[(dataset$TOB_EVER_USED_11 %in% c(0))==TRUE] = 0
dataset$TOB_EVER_USEDFULL_11[(dataset$TOB_EVER_USED_11 %in% c(1)&dataset$TOB_AGE_1st_USED_11 %in% c(0,1))==TRUE] = 0
dataset$TOB_EVER_USEDFULL_11[(dataset$TOB_EVER_USED_11 %in% c(1))==TRUE & (dataset$TOB_AGE_1st_USED_11 %in% c(0,1))==FALSE] = 1

dataset$TOB_EVER_USEDFULL_14 = NA
dataset$TOB_EVER_USEDFULL_14[(dataset$TOB_EVER_USED_14 %in% c(0))==TRUE] = 0
dataset$TOB_EVER_USEDFULL_14[(dataset$TOB_EVER_USED_14 %in% c(1)&dataset$TOB_AGE_1st_USED_14 %in% c(0,1))==TRUE] = 0
dataset$TOB_EVER_USEDFULL_14[(dataset$TOB_EVER_USED_14 %in% c(1))==TRUE & (dataset$TOB_AGE_1st_USED_14 %in% c(0,1))==FALSE] = 1

dataset$TOB_EVER_USEDFULL_17 = NA
dataset$TOB_EVER_USEDFULL_17[(dataset$TOB_EVER_USED_17 %in% c(0))==TRUE] = 0
dataset$TOB_EVER_USEDFULL_17[(dataset$TOB_EVER_USED_17 %in% c(1)&dataset$TOB_AGE_1st_USED_17 %in% c(0,1))==TRUE] = 0
dataset$TOB_EVER_USEDFULL_17[(dataset$TOB_EVER_USED_17 %in% c(1))==TRUE & (dataset$TOB_AGE_1st_USED_17 %in% c(0,1))==FALSE] = 1

dataset$TOB_EVER_USEDFULL_1114[(dataset$TOB_EVER_USEDFULL_11 %in% c(0)&dataset$TOB_EVER_USEDFULL_14 %in% c(0))==TRUE] = 0
dataset$TOB_EVER_USEDFULL_1114[(dataset$TOB_EVER_USEDFULL_11 %in% c(0)&dataset$TOB_EVER_USEDFULL_14 %in% c(0))==FALSE] = 1

# out of the total 594, 9 reported + tob use by age-14 but - at 17; recode those to 1 at age-17
table(dataset$TOB_EVER_USEDFULL_1114)
dataset$TOB_EVER_USEDFULL_17[(dataset$TOB_EVER_USEDFULL_1114 %in% c(1)&dataset$TOB_EVER_USEDFULL_17 %in% c(0))==TRUE] = 1

dataset$TOB_EVER_USED_1114[(dataset$TOB_EVER_USED_11 %in% c(0)&dataset$TOB_EVER_USED_14 %in% c(0))==TRUE] = 0
dataset$TOB_EVER_USED_1114[(dataset$TOB_EVER_USED_11 %in% c(0)&dataset$TOB_EVER_USED_14 %in% c(0))==FALSE] = 1

 # out of the total 594, 7 reported + tob use by age-14 but - at 17; recode those to 1 at age-17
table(dataset$TOB_EVER_USED_1114,dataset$TOB_EVER_USED_17)
dataset$TOB_EVER_USED_17[(dataset$TOB_EVER_USED_1114 %in% c(1)&dataset$TOB_EVER_USED_17 %in% c(0))==TRUE] = 1
table(dataset$TOB_EVER_USED_1114,dataset$TOB_EVER_USED_17)

dataset$Drug_EVER_USED_1114 = NULL
dataset$Drug_EVER_USED_1114[(dataset$Number_Drugs_Tried_11>0) |(dataset$Number_Drugs_Tried_14>0)] = 1
dataset$Drug_EVER_USED_1114[(dataset$Number_Drugs_Tried_11==0) &(dataset$Number_Drugs_Tried_14==0)] = 0

# Fix logical inconsistencies
dataset$CAN_EVER_USED_11[dataset$ID %in% 33333333] = 0 # age of can onset = 98, 0 frequency score, no use at 14
dataset$CAN_EVER_USED_14[dataset$ID %in% 44444444] = 0 # @14 age of can onset = 98, 0 frequency score; @17 age of initiation is 15
dataset$CAN_EVER_USED_11[dataset$ID %in% 55555555] = 0 # age of can onset = 98, 0 frequency score, no use at 14 or 17
dataset$CAN_EVER_USED_11[dataset$ID %in% 66666666] = 0 # age of can onset = 98, 0 frequency score, no use at 14 or 17

dataset$CAN_EVER_USED_1114[(dataset$CAN_EVER_USED_11 %in% c(0)&dataset$CAN_EVER_USED_14 %in% c(0))==TRUE] = 0
dataset$CAN_EVER_USED_1114[(dataset$CAN_EVER_USED_11 %in% c(0)&dataset$CAN_EVER_USED_14 %in% c(0))==FALSE] = 1
table(dataset$CAN_EVER_USED_1114,dataset$CAN_EVER_USED_17+10) # out of 594, 4 reported use by 14 but none at 17, recode those to 1 at 17

dataset$CAN_EVER_USED_17[(dataset$CAN_EVER_USED_1114 %in% c(1)&dataset$CAN_EVER_USED_17 %in% c(0))==TRUE] = 1
table(dataset$CAN_EVER_USED_1114,dataset$CAN_EVER_USED_17+10)

dataset$DrugCAN_EVER_USED_1114 = NULL
dataset$DrugCAN_EVER_USED_1114[(dataset$Number_Drugs_Tried_11>0) |(dataset$Number_Drugs_Tried_14>0) | dataset$CAN_EVER_USED_11 %in% 1 | dataset$CAN_EVER_USED_14 %in% 1 ] = 1
dataset$DrugCAN_EVER_USED_1114[(dataset$Number_Drugs_Tried_11 %in% 0) &(dataset$Number_Drugs_Tried_14 %in% 0)&dataset$CAN_EVER_USED_11 %in% 0 & dataset$CAN_EVER_USED_14 %in% 0] = 0


# Parental AUD status
ParentIN.data = read.table("../ES/data/PARENTS_IN_AUD_wide.dat", header=TRUE, sep="\t", na=c(-999));
dataset = merge(dataset, ParentIN.data,by="idyrfam",all.x=TRUE);
rm(list=setdiff(ls(), "dataset"));

ParentFU2.data = read.table("../ES/data/PARENTS_FU2_AUD_wide.dat", header=TRUE, sep="\t", na=c(-999));
dataset = merge(dataset, ParentFU2.data,by="idyrfam",all.x=TRUE);
rm(list=setdiff(ls(), "dataset"));

# Code as 1 if either parent met criteria for AUD at IN or FU2
dataset$ParINAUD = apply(data.frame(dataset$IN_AUD04_10,dataset$IN_AUD04_11,dataset$IN_AUD04_15, dataset$IN_AUD04_16,dataset$IN_AUD04_36), 1, function(x) as.numeric(length(which(x==1)) >=1))
dataset$ParFU2AUD = apply(data.frame(dataset$FU2_AUD64_10,dataset$FU2_AUD64_11,dataset$FU2_AUD64_15, dataset$FU2_AUD64_16,dataset$FU2_AUD64_36), 1, function(x) as.numeric(length(which(x==1)) >=1))

dataset$ParINFU2AUD = apply(data.frame(dataset$ParINAUD,dataset$ParFU2AUD), 1, function(x) as.numeric(length(which(x==1)) >=1))


# make new drug measure without pill question (questionable validity of the pills item, so do not use)
dataset$Number_Drugs_Tried_NoD7_11 = integer(dim(dataset)[1])-999
dataset$Number_Drugs_Tried_NoD7_11 = apply(data.frame(dataset$D1_11,dataset$D2_11,dataset$D3_11,dataset$D6_11,dataset$D8_11,dataset$D9_11,dataset$D10_11,dataset$D11_11),1, function(x) sum(x))

missidx = as.numeric(apply(data.frame(dataset$D1_11,dataset$D2_11,dataset$D3_11,dataset$D6_11,dataset$D8_11,dataset$D9_11,dataset$D10_11,dataset$D11_11),1, function(x){any(is.na(x))}))

NAIDs = dataset$ID[missidx == 1]

for(i in seq(along= NAIDs)) {
 	idx = match(NAIDs[i],dataset$ID); 	if(rowSums(data.frame(dataset$D1_11,dataset$D2_11,dataset$D3_11,dataset$D6_11,dataset$D8_11,dataset$D9_11,dataset$D10_11,dataset$D11_11)[idx,],na.rm=T) >= 1) {
 		dataset$Number_Drugs_Tried_NoD7_11[idx] = rowSums(data.frame(dataset$D1_11,dataset$D2_11,dataset$D3_11,dataset$D6_11,dataset$D8_11,dataset$D9_11,dataset$D10_11,dataset$D11_11)[idx,],na.rm=T)  }
                                      }

dataset$Number_Drugs_Tried_NoD7_14 = integer(dim(dataset)[1])-999
dataset$Number_Drugs_Tried_NoD7_14 = apply(data.frame(dataset$D1_14,dataset$D2_14,dataset$D3_14,dataset$D6_14,dataset$D8_14,dataset$D9_14,dataset$D10_14,dataset$D11_14),1, function(x) sum(x))

missidx = as.numeric(apply(data.frame(dataset$D1_14,dataset$D2_14,dataset$D3_14,dataset$D6_14,dataset$D8_14,dataset$D9_14,dataset$D10_14,dataset$D11_14),1, function(x){any(is.na(x))}))

NAIDs = dataset$ID[missidx == 1]

for(i in seq(along= NAIDs)) {
 	idx = match(NAIDs[i],dataset$ID); 	if(rowSums(data.frame(dataset$D1_14,dataset$D2_14,dataset$D3_14,dataset$D6_14,dataset$D8_14,dataset$D9_14,dataset$D10_14,dataset$D11_14)[idx,],na.rm=T) >= 1) {
 		dataset$Number_Drugs_Tried_NoD7_14[idx] = rowSums(data.frame(dataset$D1_14,dataset$D2_14,dataset$D3_14,dataset$D6_14,dataset$D8_14,dataset$D9_14,dataset$D10_14,dataset$D11_14)[idx,],na.rm=T)  }
                                      }
dataset$Number_Drugs_Tried_NoD7_17 = integer(dim(dataset)[1])-999
dataset$Number_Drugs_Tried_NoD7_17 = apply(data.frame(dataset$D1_17,dataset$D2_17,dataset$D3_17,dataset$D6_17,dataset$D8_17,dataset$D9_17,dataset$D10_17,dataset$D11_17),1, function(x) sum(x))

missidx = as.numeric(apply(data.frame(dataset$D1_17,dataset$D2_17,dataset$D3_17,dataset$D6_17,dataset$D8_17,dataset$D9_17,dataset$D10_17,dataset$D11_17),1, function(x){any(is.na(x))}))

NAIDs = dataset$ID[missidx == 1]

for(i in seq(along= NAIDs)) {
 	idx = match(NAIDs[i],dataset$ID); 	if(rowSums(data.frame(dataset$D1_17,dataset$D2_17,dataset$D3_17,dataset$D6_17,dataset$D8_17,dataset$D9_17,dataset$D10_17,dataset$D11_17)[idx,],na.rm=T) >= 1) {
 		dataset$Number_Drugs_Tried_NoD7_17[idx] = rowSums(data.frame(dataset$D1_17,dataset$D2_17,dataset$D3_17,dataset$D6_17,dataset$D8_17,dataset$D9_17,dataset$D10_17,dataset$D11_17)[idx,],na.rm=T)  }
                                      }
                                      
dataset$DrugNoD7_EVER_USED_11 = NULL
dataset$DrugNoD7_EVER_USED_11[(dataset$Number_Drugs_Tried_NoD7_11>0)] = 1
dataset$DrugNoD7_EVER_USED_11[(dataset$Number_Drugs_Tried_NoD7_11==0)] = 0

dataset$DrugNoD7_EVER_USED_14 = NULL
dataset$DrugNoD7_EVER_USED_14[(dataset$Number_Drugs_Tried_NoD7_14>0)] = 1
dataset$DrugNoD7_EVER_USED_14[(dataset$Number_Drugs_Tried_NoD7_14==0)] = 0

dataset$DrugNoD7_EVER_USED_17 = NULL
dataset$DrugNoD7_EVER_USED_17[(dataset$Number_Drugs_Tried_NoD7_17>0)] = 1
dataset$DrugNoD7_EVER_USED_17[(dataset$Number_Drugs_Tried_NoD7_17==0)] = 0

dataset$DrugNoD7_EVER_USED_1114 = NULL
dataset$DrugNoD7_EVER_USED_1114[(dataset$Number_Drugs_Tried_NoD7_11>0) |(dataset$Number_Drugs_Tried_NoD7_14>0)] = 1
dataset$DrugNoD7_EVER_USED_1114[(dataset$Number_Drugs_Tried_NoD7_11==0) &(dataset$Number_Drugs_Tried_NoD7_14==0)] = 0

dataset$DrugNoD7CAN_EVER_USED_1114 = NULL
dataset$DrugNoD7CAN_EVER_USED_1114[(dataset$Number_Drugs_Tried_NoD7_11>0) |(dataset$Number_Drugs_Tried_NoD7_14>0) | dataset$CAN_EVER_USED_11 %in% 1 | dataset$CAN_EVER_USED_14 %in% 1 ] = 1
dataset$DrugNoD7CAN_EVER_USED_1114[(dataset$Number_Drugs_Tried_NoD7_11 %in% 0) &(dataset$Number_Drugs_Tried_NoD7_14 %in% 0)&dataset$CAN_EVER_USED_11 %in% 0 & dataset$CAN_EVER_USED_14 %in% 0] = 0

# calculate within-cluster mean and difference scores for TPROIthmFCZT
fubar = data.frame(dataset$subname, dataset$idyrfam, dataset$TPROIthmFCZT);
names(fubar)[names(fubar)=="dataset.subname"] = "subname";
names(fubar)[names(fubar)=="dataset.idyrfam"] = "idyrfam";
names(fubar)[names(fubar)=="dataset.TPROIthmFCZT"] = "TPROIthmFCZT";
twin.means = aggregate(fubar$TPROIthmFCZT, by=list(idyrfam=fubar$idyrfam), FUN=mean, na.rm=TRUE);
fubar = merge(fubar, twin.means, by="idyrfam");
fubar$TPROIthmFCZTdiff = fubar$TPROIthmFCZT - fubar$x;
names(fubar)[names(fubar)=="x"] = "TPROIthmFCZTmean";
fubar$TPROIthmFCZT = NULL;
dataset = merge(dataset, fubar, by = c("subname", "idyrfam")); 

# Remove IDs with ACC <50% in either condition
dataset$ID[((dataset$target_correct_per)<50|(dataset$freq_correct_per)<50)];
dataset = dataset[!dataset$ID %in% dataset$ID[((dataset$target_correct_per)<50|(dataset$freq_correct_per)<50)],]
dataset$ID[((dataset$target_correct_per)<50|(dataset$freq_correct_per)<50)];

}

# End data cleaning


# add misc functions
# R2 analog for linear mixed models (lmer objects)
Snijders_Bosker <-
function(mod,NULL_mod){
	sigma_1 = mod$vcov[2]
	tau_1   = mod$vcov[1]
	sigma_0 = NULL_mod$vcov[2]
	tau_0   = NULL_mod$vcov[1]
	R_squared = 1 - (sigma_1 + tau_1)/(sigma_0 + tau_0)
	print(cbind(R_squared, sigma_1, tau_1, sigma_0, tau_0));
					}

# Transform unstandardized beta to standardized beta in linear mixed models (lmer)
lm.beta.lmer <- function(mod) {
   b <- fixef(mod)[-1]
   sd.x <- apply(cbind(getME(mod,"X")[,-1]),2,sd) # added cbind to account for 1-IV models 27/02/18
  #sd.x <- apply(getME(mod,"X")[,-1],2,sd) # 
   sd.y <- sd(getME(mod,"y"))
   b*sd.x/sd.y
}


summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,conf.interval=.95, .drop=TRUE) {
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
 
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
    .fun = function(xx, col) {
      c(N    = length2(xx[[col]], na.rm=na.rm),
        mean = mean   (xx[[col]], na.rm=na.rm),
        sd   = sd     (xx[[col]], na.rm=na.rm),
        min  = min    (xx[[col]], na.rm=na.rm),
        max  = max    (xx[[col]], na.rm=na.rm))},
    measurevar)
 
  # Rename the "mean" column    
  datac <- plyr::rename(datac, c("mean" = measurevar))
 
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
 
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
 
  return(datac)
}

# Calculate contrast matrix for backward difference coding (e.g., 2 vs. 1, 3 vs. 2, etc)
backward_difference_coding <- function(lvls){
  conts <- lvls - 1
  contrast_mat <- matrix(data = NA, nrow = lvls, ncol = conts)
  for(i in seq(conts)) {
    cont_vec <- c(rep( -(lvls - i)/lvls , i) , rep(i/lvls , lvls-i))
    contrast_mat[,i] <- cont_vec
  }
  print(contrast_mat)
  return(contrast_mat)
}

cbPalette <- c("#0072B2", "#009E73","#D55E00", "#999999", "#E69F00", "#56B4E9",  "#F0E442", "#CC79A7")


## Begin Analysis ##

# Calculate Cronbach's alpha to assess internal consistency of the composite measure at ages 11/14 and 17
psych::alpha(as.data.frame(cbind(dataset$ALC_INTOX_12m_bin_1114, dataset$ALC_MAXDRINK_12m_bin_1114, dataset$ALC_AMT_12m_bin_1114, dataset$ALC_FREQ_12m_bin_1114)))
cormatrix = cor(as.data.frame(cbind(dataset$ALC_INTOX_12m_bin_1114, dataset$ALC_MAXDRINK_12m_bin_1114, dataset$ALC_AMT_12m_bin_1114, dataset$ALC_FREQ_12m_bin_1114)),use='complete.obs')
cormatrix
mean(cormatrix[lower.tri(cor(cormatrix))])
range(cormatrix[lower.tri(cor(cormatrix))])

psych::alpha(as.data.frame(cbind(dataset$ALC_INTOX_12m_17_bin, dataset$ALC_MAXDRINK_12m_17_bin, dataset$ALC_AMT_12m_17_bin, dataset$ALC_FREQ_12m_17_bin)))
cormatrix = cor(as.data.frame(cbind(dataset$ALC_INTOX_12m_17_bin, dataset$ALC_MAXDRINK_12m_17_bin, dataset$ALC_AMT_12m_17_bin, dataset$ALC_FREQ_12m_17_bin)),use='complete.obs')
cormatrix
mean(cormatrix[lower.tri(cor(cormatrix))])
range(cormatrix[lower.tri(cor(cormatrix))])


# Descriptive statistics and change over time

# Restructure data into wide format to visualize and test differences over time
melted = NULL;
melted = melt(dataset, id.vars=c("ID", "idyrfam","sex"), measure.vars = c("CSU_drink_index_11", "CSU_drink_index_14","CSU_drink_index_17"));
melted$age = as.factor(c(rep('11', length(melted$variable[melted$variable=="CSU_drink_index_11"])),rep('14', length(melted$variable[melted$variable=="CSU_drink_index_14"])), rep('17', length(melted$variable[melted$variable=="CSU_drink_index_17"]))));
melted = plyr::rename(melted,c("value" = "CSU_drink_index"));

# Set up backward difference contrast matrix
contrasts(melted$age) = backward_difference_coding(length(unique(melted$age)))

# Descriptives for full sample
tgc <- summarySE(melted, measurevar="CSU_drink_index", groupvars=c("age"))

# Descriptives for only Age-17 initiators
tgc_drink17 <- summarySE(melted[melted$ID %in% dataset[dataset$ALC_EVER_USED_NOPERM_17==1,'ID'],], measurevar="CSU_drink_index", groupvars=c("age"))

tgc$Group = 'All'
tgc_drink17$Group = 'FU2-Initiators'
tgc = rbind(tgc,tgc_drink17)
tgc

# Plot mean and 95% CI for each assessment, split by group, and plot jittered data points
ggplot(tgc, aes(age, CSU_drink_index,group=Group))+  geom_jitter(data = melted, aes(x=age, y = CSU_drink_index, group = 1),width=.25,height=.15,alpha=.5) +geom_errorbar(aes(ymin=CSU_drink_index-ci, ymax=CSU_drink_index+ci), width=.075) + geom_line(aes(color=Group)) + geom_point(aes(color=Group),size=3) + theme_classic() + scale_color_manual(values= cbPalette) + labs(x = "Assessment Age", y = "Score", title = "Drink Index by Age and Group") + scale_y_continuous(breaks=c(0:17))


# Test for overall differences in scores as a function of assessment (include a random effect at the family-level to account for clustered data)

CSU_LMM <- lmer(CSU_drink_index ~ as.factor(age) + (1 | idyrfam), data = melted)

anova(CSU_LMM, ddf = "Kenward-Roger")
summary(CSU_LMM, ddf = "Kenward-Roger") # provides sequential contrasts to test growth differences

visreg::visreg(CSU_LMM,"age", overlay=TRUE, type = 'conditional', line=list(col=c("red3")),fill=list(col=c("#000000")),points=list(col=c("#000000"), pch=c(1),size=1.5), xlab="Assessment Age", ylab="Drink Index Score", gg=TRUE) + ggplot2::theme_classic() + ggplot2::ggtitle( 'LMM Partial Residual Plot')


####

# Linear mixed models predicting Age-17 drinking (continuous) with Age-14 measures                        
{

# Create covariate-only model
Drink1114_ParAUD_sex_LMM.model = lmer(CSU_drink_index_17 ~ CSU_drink_index_1114 + sex + ParINFU2AUD + (1 | idyrfam), data = dataset);

# Create full model
P3Theta_LMM.model = lmer(CSU_drink_index_17 ~ TDp3mPZPOZT + TPROIthmFCZT + sex + CSU_drink_index_1114 + ParINFU2AUD +(1 | idyrfam), data = dataset);

# Compare model fit of covariate-only and full model
anova(Drink1114_ParAUD_sex_LMM.model, P3Theta_LMM.model);

summary(P3Theta_LMM.model, ddf = "Kenward-Roger");
confint(P3Theta_LMM.model) # profile-based confidence intervals
lm.beta.lmer(P3Theta_LMM.model); # standardized betas

# Calculate change in R2 following Snijders and Bosker (2012) method
# Quantifies increase in R2 for the full model relative to the covariate model
Snijders_Bosker(as.data.frame(VarCorr(P3Theta_LMM.model)),as.data.frame(VarCorr(Drink1114_ParAUD_sex_LMM.model)));

}


# Invert measures to aid interpretability of odds ratios
dataset$TDp3mPZPOZT_inv = dataset$TDp3mPZPOZT*-1
dataset$TPROIthmFCZT_inv = dataset$TPROIthmFCZT*-1

# logistic generalized linear mixed regression model predicting new development of drinking by Age-17
# Covariate-only logistic glmer
gmCov = glmer(ALC_EVER_USED_NOPERM_17 ~  sex + ParINFU2AUD + (1 | idyrfam), data = dataset[dataset$ALC_EVER_USED_NOPERM_1114 %in% c(0),], family = binomial);

# Full model logistic glmer 
gm1 = glmer(ALC_EVER_USED_NOPERM_17 ~ scale(TDp3mPZPOZT_inv) + scale(TPROIthmFCZT_inv) + sex + ParINFU2AUD + (1 | idyrfam), data = dataset[dataset$ALC_EVER_USED_NOPERM_1114 %in% c(0),], family = binomial) # control = glmerControl(optCtrl = list(maxfun = 300000), optimizer='bobyqa'), nAGQ=2)

# Compare fit
anova(gmCov,gm1)

summary(gm1)

# Get CIs for the betas, then convert to Odds Ratios
cc <- confint(gm1,parm="beta_") #,devtol=Inf);
ctab <- cbind(est=fixef(gm1),cc)
rtab <- exp(ctab)
round(rtab,digits=2)

# logistic generalized linear mixed regression model predicting Initiation of Intoxication, excluding those use before 14
gmCov = glmer(ALC_EVER_DRUNK_17 ~  sex + ParINFU2AUD + (1 | idyrfam), data = dataset[dataset$ALC_EVER_USED_NOPERM_1114 %in% c(0),], family = binomial);

gm1 = glmer(ALC_EVER_DRUNK_17 ~ scale(TDp3mPZPOZT_inv) + scale(TPROIthmFCZT_inv) + sex + ParINFU2AUD + (1 | idyrfam), data = dataset[dataset$ALC_EVER_USED_NOPERM_1114 %in% c(0),], family = binomial);

# Compare fit
anova(gmCov,gm1)

summary(gm1)

# Get CIs for the betas, then convert to Odds Ratios
cc <- confint(gm1,parm="beta_");
ctab <- cbind(est=fixef(gm1),cc)
rtab <- exp(ctab)
round(rtab,digits=2)
data.frame(r2beta(gm1, partial=TRUE, method='nsj', data = dataset)); # if using transforms (e.g., log1p) in model, need to specify data

# logistic generalized linear mixed regression model predicting Initiation of Binge (≥5 drinks), excluding those use before 14
gmCov = glmer(CSU_BINGE_17 ~  sex + ParINFU2AUD + (1 | idyrfam), data = dataset[dataset$ALC_EVER_USED_NOPERM_1114 %in% c(0),], family = binomial);

gm1 = glmer(CSU_BINGE_17 ~ scale(TDp3mPZPOZT_inv) + scale(TPROIthmFCZT_inv) + sex + ParINFU2AUD + (1 | idyrfam), data = dataset[dataset$ALC_EVER_USED_NOPERM_1114 %in% c(0),], family = binomial);

# Compare fit
anova(gmCov,gm1)

summary(gm1)

# Get CIs for the betas, then convert to Odds Ratios
cc <- confint(gm1,parm="beta_");
ctab <- cbind(est=fixef(gm1),cc)
rtab <- exp(ctab)
round(rtab,digits=2)
data.frame(r2beta(gm1, partial=TRUE, method='nsj', data = dataset)); # if using transforms (e.g., log1p) in model, need to specify data


# Plot the ORs and CIs for the logistic glmer models
target=NULL
target$Model = as.character(c("Ever used alcohol?",
"Ever used alcohol?", 
"Ever been drunk/intoxicated?",
"Ever been drunk/intoxicated?",
"Ever binged? (≥ 5 drinks on one occasion)",
"Ever binged? (≥ 5 drinks on one occasion)"))
target$Predictor = as.character(c("P3", "Theta","P3", "Theta","P3", "Theta"))
target$Est  = as.numeric(c(1.49,1.48,1.68,2.82,2.73,2.50))
target$low  = as.numeric(c(1.09,1.08,0.81,1.32,1.26,1.20))
target$high = as.numeric(c(2.16,2.16,3.80,7.11,6.67,5.88))
target = data.frame(target)

adj=.2
ggplot(data = target, aes(y= Predictor, x=Est, color = Model)) +
  geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed") +
  geom_errorbarh(data = dplyr::filter(target, Model== target$Model[1]), aes(xmax = high, xmin = low), size = .5, height = .1, color = "gray50", position = position_nudge(y = adj)) +
  geom_point(data = dplyr::filter(target, Model == target$Model[1]), size = 4, position = position_nudge(y = adj)) +
  geom_errorbarh(data = dplyr::filter(target, Model == target$Model[3]), aes(xmax = high, xmin = low), size = .5, height = .1, color = "gray50") +
  geom_point(data = dplyr::filter(target, Model == target$Model[3]), size = 4) +
  geom_errorbarh(data = dplyr::filter(target, Model == target$Model[5]), aes(xmax = high, xmin = low), size = .5, height = .1, color = "gray50", position = position_nudge(y = - adj)) +
  geom_point(data = dplyr::filter(target, Model == target$Model[5]), size = 4, position = position_nudge(y = - adj)) + theme_bw() + scale_x_log10("Odds Ratio",breaks=c(1.00, 1.38, 1.91, 2.65, 3.66, 5.06, 7.00), limits=c(.75,7.5))  + theme(panel.grid.minor = element_blank()) 


# Supplemental analyses are below

# DRINK INDEX - ONLY NO DRINK BY 14

Cov_LMM.model = lmer(CSU_drink_index_17 ~ sex + ParINFU2AUD + (1 | idyrfam), data = dataset[dataset$ALC_EVER_USED_NOPERM_1114 %in% c(0),]);

P3Theta_LMM.model = lmer(CSU_drink_index_17 ~ TDp3mPZPOZT + TPROIthmFCZT + sex + ParINFU2AUD + (1 | idyrfam), data = dataset[dataset$ALC_EVER_USED_NOPERM_1114 %in% c(0),]);

anova(Cov_LMM.model, P3Theta_LMM.model)

summary(P3Theta_LMM.model, ddf = "Kenward-Roger");
lm.beta.lmer(P3Theta_LMM.model);

Snijders_Bosker(as.data.frame(VarCorr(P3Theta_LMM.model)),as.data.frame(VarCorr(Drink1114_ParAUD_sex_LMM.model)));


# Re-run the models to test if effects hold after adjusting for potentially important covariate (other drug use)

dataset$TOBCANDrugNoD7_1114 = NULL
dataset$TOBCANDrugNoD7_1114[dataset$DrugNoD7CAN_EVER_USED_1114 %in% c(0) & dataset$TOB_EVER_USEDFULL_1114 %in% c(0)] = 0
dataset$TOBCANDrugNoD7_1114[dataset$DrugNoD7CAN_EVER_USED_1114 %in% c(1) | dataset$TOB_EVER_USEDFULL_1114 %in% c(1)] = 1

datasetCov = dataset[complete.cases(dataset$DrugNoD7CAN_EVER_USED_1114)&complete.cases(dataset$TOB_EVER_USEDFULL_1114),]
table(datasetCov$DrugNoD7CAN_EVER_USED_1114)
table(datasetCov$TOB_EVER_USEDFULL_1114)
 

P3Theta_LMM.model = lmer(CSU_drink_index_17 ~ TDp3mPZPOZT + TPROIthmFCZT + sex + CSU_drink_index_1114 + ParINFU2AUD + TOB_EVER_USEDFULL_1114 + DrugNoD7CAN_EVER_USED_1114 + (1 | idyrfam), data = datasetCov);
summary(P3Theta_LMM.model, ddf = "Kenward-Roger");
lm.beta.lmer(P3Theta_LMM.model)

table(datasetCov$DrugNoD7CAN_EVER_USED_1114[datasetCov$ALC_EVER_USED_NOPERM_1114 %in% c(0)])
table(datasetCov$TOB_EVER_USEDFULL_1114[datasetCov$ALC_EVER_USED_NOPERM_1114 %in% c(0)])
table(datasetCov$TOBCANDrugNoD7_1114[datasetCov$ALC_EVER_USED_NOPERM_1114 %in% c(0)])

  # Predicting Initiation of Alcohol Use, excluding those use before 14
gm1 = glmer(ALC_EVER_USED_NOPERM_17 ~ scale(TDp3mPZPOZT_inv) + scale(TPROIthmFCZT_inv) + sex + ParINFU2AUD + TOBCANDrugNoD7_1114 + (1 | idyrfam), data = datasetCov[datasetCov$ALC_EVER_USED_NOPERM_1114 %in% c(0),], family = binomial);
summary(gm1)
cc <- confint(gm1,parm="beta_");
ctab <- cbind(est=fixef(gm1),cc)
rtab <- exp(ctab)
round(rtab,digits=2)

  # Predicting Initiation of Intoxication, excluding those use before 14
gm1 = glmer(ALC_EVER_DRUNK_17 ~ scale(TDp3mPZPOZT_inv) + scale(TPROIthmFCZT_inv) + sex + ParINFU2AUD + TOBCANDrugNoD7_1114  + (1 | idyrfam), data = datasetCov[datasetCov$ALC_EVER_USED_NOPERM_1114 %in% c(0),], family = binomial);
summary(gm1)
cc <- confint(gm1,parm="beta_");
ctab <- cbind(est=fixef(gm1),cc)
rtab <- exp(ctab)
round(rtab,digits=2)

  # Predicting Initiation of Binge (≥5 drinks), excluding those use before 14
gm1 = glmer(CSU_BINGE_17 ~ scale(TDp3mPZPOZT_inv) + scale(TPROIthmFCZT_inv) + sex + ParINFU2AUD + TOBCANDrugNoD7_1114  + (1 | idyrfam), data = datasetCov[datasetCov$ALC_EVER_USED_NOPERM_1114 %in% c(0),], family = binomial);
summary(gm1)
cc <- confint(gm1,parm="beta_");
ctab <- cbind(est=fixef(gm1),cc)
rtab <- exp(ctab)
round(rtab,digits=2)


P3Theta_LMM.model = lmer(CSU_drink_index_17 ~ TDp3mPZPOZT + TPROIthmFCZT + sex  + ParINFU2AUD + TOBCANDrugNoD7_1114 + (1 | idyrfam), data = datasetCov[datasetCov$ALC_EVER_USED_NOPERM_1114 %in% c(0),]);
summary(P3Theta_LMM.model, ddf = "Kenward-Roger");
lm.beta.lmer(P3Theta_LMM.model)


P3Theta_LMM.model = lmer(CSU_drink_index_17 ~ TDp3mPZPOZT + TPROIthmFCZT + sex  + ParINFU2AUD  + (1 | idyrfam), data = datasetCov[datasetCov$ALC_EVER_USED_NOPERM_1114 %in% c(0)&datasetCov$TOBCANDrugNoD7_1114 %in% c(0),]);
summary(P3Theta_LMM.model, ddf = "Kenward-Roger");



# Behavior

dataset$target_correct_per    = dataset$target_correct_per/100
dataset$freq_error_per        = dataset$freq_error_per/100
dataset$target_correct_meanRT = dataset$target_correct_meanRT/100
dataset$target_correct_sdRT   = dataset$target_correct_sdRT/100


BEH_LMM.model = lmer(CSU_drink_index_17 ~ asin(sqrt(target_correct_per)) + target_correct_meanRT + target_correct_sdRT + sex + CSU_drink_index_1114 + ParINFU2AUD+ (1 | idyrfam), data = dataset);
summary(BEH_LMM.model, ddf = "Kenward-Roger");

gm1 = glmer(ALC_EVER_USED_NOPERM_17 ~ asin(sqrt(target_correct_per)) + target_correct_meanRT + target_correct_sdRT + sex + ParINFU2AUD+(1 | idyrfam), data = dataset[dataset$ALC_EVER_USED_NOPERM_1114 %in% c(0),], family = binomial);
summary(gm1)
cc <- confint(gm1,parm="beta_");
ctab <- cbind(est=fixef(gm1),cc)
rtab <- exp(ctab)
round(rtab,digits=2)

  # Predicting Initiation of Intoxication, excluding those use before 14
gm1 = glmer(ALC_EVER_DRUNK_17 ~ asin(sqrt(target_correct_per)) + target_correct_meanRT + target_correct_sdRT + sex +ParINFU2AUD+ (1 | idyrfam), data = dataset[dataset$ALC_EVER_USED_NOPERM_1114 %in% c(0),], family = binomial);
summary(gm1)
cc <- confint(gm1,parm="beta_");
ctab <- cbind(est=fixef(gm1),cc)
rtab <- exp(ctab)
round(rtab,digits=2)

  # Predicting Initiation of Binging Past 12m, excluding those use before 14
gm1 = glmer(CSU_BINGE_17 ~ asin(sqrt(target_correct_per)) + target_correct_meanRT + target_correct_sdRT + sex + ParINFU2AUD+(1 | idyrfam), data = dataset[dataset$ALC_EVER_USED_NOPERM_1114 %in% c(0),], family = binomial);
summary(gm1)
cc <- confint(gm1,parm="beta_");
ctab <- cbind(est=fixef(gm1),cc)
rtab <- exp(ctab)
round(rtab,digits=2)


# Rerun all models testing for an interaction between main predictors and sex
{

# DRINK INDEX
P3Theta_LMM.model = lmer(CSU_drink_index_17 ~ TDp3mPZPOZT*sex + TPROIthmFCZT*sex + CSU_drink_index_1114 + ParINFU2AUD +(1 | idyrfam), data = dataset);
summary(P3Theta_LMM.model, ddf = "Kenward-Roger");
lm.beta.lmer(P3Theta_LMM.model);

P3Theta_LMM_NOINT.model = lmer(CSU_drink_index_17 ~ TDp3mPZPOZT + TPROIthmFCZT+sex + CSU_drink_index_1114+ParINFU2AUD + (1 | idyrfam), data = dataset);

# Test difference in model with with or without interaction
anova(P3Theta_LMM_NOINT.model, P3Theta_LMM.model);

}

  # Predicting Initiation of Alcohol Use, excluding those use before 14
gm1 = glmer(ALC_EVER_USED_NOPERM_17 ~ scale(TDp3mPZPOZT_inv)*sex + scale(TPROIthmFCZT_inv) *sex + ParINFU2AUD + (1 | idyrfam), data = dataset[dataset$ALC_EVER_USED_NOPERM_1114 %in% c(0),], family = binomial);
summary(gm1)
cc <- confint(gm1,parm="beta_");
ctab <- cbind(est=fixef(gm1),cc)
rtab <- exp(ctab)
round(rtab,digits=2)

gm1NOINT = glmer(ALC_EVER_USED_NOPERM_17 ~ scale(TDp3mPZPOZT_inv) + scale(TPROIthmFCZT_inv) + sex + ParINFU2AUD + (1 | idyrfam), data = dataset[dataset$ALC_EVER_USED_NOPERM_1114 %in% c(0),], family = binomial);
anova(gm1, gm1NOINT)

  # Predicting Initiation of Intoxication, excluding those use before 14
gm1 = glmer(ALC_EVER_DRUNK_17 ~ scale(TDp3mPZPOZT_inv)*sex + scale(TPROIthmFCZT_inv)*sex + ParINFU2AUD + (1 | idyrfam), data = dataset[dataset$ALC_EVER_USED_NOPERM_1114 %in% c(0),], family = binomial);
summary(gm1)
cc <- confint(gm1,parm="beta_");
ctab <- cbind(est=fixef(gm1),cc)
rtab <- exp(ctab)
round(rtab,digits=2)

gm1NOINT = glmer(ALC_EVER_DRUNK_17 ~ scale(TDp3mPZPOZT_inv) + scale(TPROIthmFCZT_inv) + sex + ParINFU2AUD + (1 | idyrfam), data = dataset[dataset$ALC_EVER_USED_NOPERM_1114 %in% c(0),], family = binomial);
anova(gm1, gm1NOINT)


  # Predicting Initiation of Binge (≥5 drinks), excluding those use before 14
gm1 = glmer(CSU_BINGE_17 ~ scale(TDp3mPZPOZT_inv)*sex + scale(TPROIthmFCZT_inv)*sex + ParINFU2AUD + (1 | idyrfam), data = dataset[dataset$ALC_EVER_USED_NOPERM_1114 %in% c(0),], family = binomial);
summary(gm1)
cc <- confint(gm1,parm="beta_");
ctab <- cbind(est=fixef(gm1),cc)
rtab <- exp(ctab)
round(rtab,digits=2)

gm1NOINT = glmer(CSU_BINGE_17 ~ scale(TDp3mPZPOZT_inv) + scale(TPROIthmFCZT_inv) + sex + ParINFU2AUD + (1 | idyrfam), data = dataset[dataset$ALC_EVER_USED_NOPERM_1114 %in% c(0),], family = binomial);

anova(gm1,gm1NOINT)


P3Theta_LMM.model = lmer(CSU_drink_index_17 ~ TDp3mPZPOZT*sex + TPROIthmFCZT*sex + ParINFU2AUD + (1 | idyrfam), data = dataset[dataset$ALC_EVER_USED_NOPERM_1114 %in% c(0),]);
summary(P3Theta_LMM.model, ddf = "Kenward-Roger");
P3Theta_LMM_NOINT.model = lmer(CSU_drink_index_17 ~ TDp3mPZPOZT + TPROIthmFCZT+sex + ParINFU2AUD + (1 | idyrfam), data = dataset[dataset$ALC_EVER_USED_NOPERM_1114 %in% c(0),]);
anova(P3Theta_LMM_NOINT.model, P3Theta_LMM.model);

