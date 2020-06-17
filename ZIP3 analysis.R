library(data.table)
library(lubridate)
library(readxl)
library(ggplot2)
library(foreign)



######
# Read in 2010 census by ZIP3
#######
zip_pop <- read.dbf("data/ZCTA_2010Census_DP1.dbf")
zip_pop <- data.table(zip_pop)
zip_pop <- zip_pop[ , c("ZCTA5CE10", "DP0010001")]
colnames(zip_pop) <- c("zip5", "pop2010")
zip_pop[ , zip3 := substr(zip5, 1, 3)]
zip_pop <- zip_pop[ , sum(pop2010), by = zip3]
colnames(zip_pop) <- c("zip3", "pop2010")
dt <- expand.grid(zip3 = zip_pop$zip3, season = c("H", "L"), year = c(2017, 2018, 2019), TTI = c("WNV", "ZIKV"))
dt <- zip_pop[dt, on = "zip3"]


###
# Read in AABB NAT testing yield data
########
zikv_yield <- data.table(read_excel("data/zikv_yield_AABB.xls"))
zikv_yield[ , "collection_date" := parse_date_time(collection_date, orders = c("mdy"))]
zikv_yield[ , interpretation_simp := ifelse(interpretation %in% c("Confirmed", "False Positive"), interpretation, "Cannot Conclude")]
zikv_yield[ , TTI := "ZIKV"]
# table(zikv_yield$interpretation_simp, year(zikv_yield$collection_date))


wnv_yield <- fread("data/wnv_yield_AABB.csv")
wnv_yield <- wnv_yield[nchar(zip5) == 5]
wnv_yield[, collection_date := parse_date_time(collection_date, orders = c("mdy"))]
wnv_yield[ , interpretation := ifelse(interpretation == "Confirmed Positive", "Confirmed", interpretation)]
wnv_yield[ , interpretation_simp := ifelse(interpretation %in% c("Confirmed", "False Positive"), interpretation, "Cannot Conclude")]
wnv_yield[ , TTI := "WNV"]
# table(wnv_yield$interpretation_simp, year(wnv_yield$collection_date))

tti_yield <- rbind(zikv_yield[ , c("TTI", "collection_date", "zip5", "interpretation_simp")],
                   wnv_yield[ , c("TTI", "collection_date", "zip5", "interpretation_simp")])
tti_yield[ , year := year(collection_date)]
tti_yield[ , month := month(collection_date)]
# table(tti_yield[interpretation_simp == "Confirmed"]$month)
tti_yield[ , season := ifelse(month >= 6 & month <= 11, "H", "L")]
# table(tti_yield$season, tti_yield$interpretation_simp)

PPV_by_year <- tti_yield[ , .N, by = .(TTI, interpretation_simp, year(collection_date))]
PPV_by_year <- dcast(data = PPV_by_year, TTI + year ~ interpretation_simp, value.var = "N")
PPV_by_year[ , PPV := ifelse(is.na(Confirmed/(Confirmed + `False Positive`)), 0, Confirmed/(Confirmed + `False Positive`))]

#setkey(tti_yield, TTI, year)
#setkey(PPV_by_year, TTI, year)
tti_yield <- tti_yield[PPV_by_year[ , c("TTI", "year", "PPV")], on = c("TTI", "year")]
tti_yield[ , case := ifelse(interpretation_simp == "Confirmed", 1, ifelse(interpretation_simp == "False Positive", 0, PPV))]
tti_yield[ , zip3 := substr(zip5, 1, 3)]


tti_yield <- tti_yield[ , sum(case), by = .(TTI, zip3, season, year)][V1 != 0]
colnames(tti_yield)[5] <- "n_cases"

rm(PPV_by_year); rm(wnv_yield); rm(zikv_yield)
dt <- tti_yield[dt, on = c("TTI", "zip3", "season", "year")]


######
# Read in states name and id by zip
#######
zip_lookup <- fread("data/zip_lookup.csv")
zip_lookup[ , zip5 := unlist(lapply(as.character(zip5), function(x) paste0(c(rep(0, times = 5 - nchar(x)), x), collapse = '')  ))]
zip_lookup[ , zip3 := substr(zip5, 1, 3)]
zip_lookup <- unique(zip_lookup[ , c("state_id", "state_name", "zip3")])
# 063 is in CT, not NY. 834 is in Idaho, not WY
zip_lookup <- zip_lookup[ !(state_id == "NY" & zip3 == "063") & !(state_id == "WY" & zip3 == "834")]

dt <- zip_lookup[dt, on = "zip3"]



######
# Read in states pop change from 2010 census
#######
state_pop_change <- data.table(read_excel("data/state_pop_change.xlsx", sheet = "Use"))
state_pop_change[ , mult2019 := pop2019/Census]
state_pop_change[ , mult2018 := pop2018/Census]
state_pop_change[ , mult2017 := pop2017/Census]
state_pop_change <- state_pop_change[ , c("Area", "mult2017", "mult2018", "mult2019")]
setnames(state_pop_change, "Area", "state_name")
state_pop_change <- melt(state_pop_change, id.vars = "state_name", value.name = "mult", variable.name = "year")
state_pop_change[ , year := as.numeric(substr(year, 5, 8))]

dt <- state_pop_change[dt, on = c("state_name", "year")]

dt[ , pop := pop2010*mult]

n_donations <- 10327046

dt[ , donations := n_donations*pop/sum(0.5*pop), by = year]
dt[ , n_cases := ifelse(is.na(n_cases), 0, n_cases)]
dt[ , prev := n_cases/donations]
dt <- dt[!is.na(prev)]
# dt[ , c("pop2010", "n_cases", "mult") := NULL]
dt <- dcast(dt, year + zip3 + state_name + state_id + season + donations ~ TTI, value.var = "prev")

dt[ , donorGroup := paste(state_id, zip3, year, season, sep = "-")]
setnames(dt, "WNV", "pWNV")
setnames(dt, "ZIKV", "pZIKV")
setnames(dt, "donations", "numDonor")
fwrite(dt, "data/donor_groups_zip3.csv")
