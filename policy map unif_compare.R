library(ggplot2)
library(data.table)
# library(mapproj)
library(scales)
library(sf)
library(rmapshaper)
library(maptools)
library(grid)
library(gtable)

library(RColorBrewer)
theme_set(theme_minimal())

# writePolyShape(test, fn = "data/dt_zip3.shp")

zip3_shapes <- readShapePoly(fn = "data/dt_zip3.shp")
zip3_sf <- st_as_sf(zip3_shapes)



dt_zip3 <- data.table(
  zip3_sf
)

setnames(dt_zip3, "SP_ID", "zip3")
dt_zip3$dummy <- NULL




# Read in PSA percents
psa_percents_zip <- fread("results/PSA_zip3_opt_by_group_combined.csv")[ , lapply(.SD, sum), .SD = 3:8, by = c("group", "year")]
psa_percents_unif <- fread("results/PSA_unif_opt_by_group_combined.csv")[ , lapply(.SD, sum), .SD = 3:8, by = c("group", "year")]

psa_percents_combined <- rbind(
  cbind(psa_percents_zip, scenario = "Main PSA scenario"),
  cbind(psa_percents_unif, scenario = "Uniform distributions")
)

psa_percents_combined[ , z := (10000 - z)/10000]
psa_percents_combined[, c("a1", "a2", "a3", "a4", "m1") :=lapply(.SD, function(x) x/10000), .SDcols=4:8]
psa_percents_combined[ , STUSPS := substr(group, 1, 2)]
psa_percents_combined[ , season := substr(group, 13, 13)]
psa_percents_combined[ , zip3 := substr(group, 4, 6)]

dt <- rbind(
  psa_percents_combined[year == "2017" & season == "H", c("scenario", "year", "season", "a1", "a2", "a3", "a4", "zip3", "STUSPS")],
  psa_percents_combined[year == "2018" & season == "H", c("scenario", "year", "season", "a1", "a2", "a3", "a4", "zip3", "STUSPS")],
  psa_percents_combined[year == "2019" & season == "H", c("scenario", "year", "season", "a1", "a2", "a3", "a4", "zip3", "STUSPS")],
  psa_percents_combined[year == "2017" & season == "L", c("scenario", "year", "season", "a1", "a2", "a3", "a4", "zip3", "STUSPS")],
  psa_percents_combined[year == "2018" & season == "L", c("scenario", "year", "season", "a1", "a2", "a3", "a4", "zip3", "STUSPS")],
  psa_percents_combined[year == "2019" & season == "L", c("scenario", "year", "season", "a1", "a2", "a3", "a4", "zip3", "STUSPS")]
)


dt <- dcast(dt, zip3 + STUSPS + year + scenario ~ season, value.var = c("a1", "a2", "a3", "a4"))
dt[, WNVH := a3_H+a4_H]
dt[, WNVL := a3_L+a4_L]
dt[, ZIKVH := a1_H+a2_H]
dt[, ZIKVL := a1_L+a2_L]

dt <- dt_zip3[dt, on = "zip3"]




zip_to_state <- unique(dt[ , c("zip3", "STUSPS")])
# dt_states <- dt_zip3[, .(geom = st_union(geometry)),by=STUSPS]
# State outlines separated
# dt_outline_48 <- dt_states[ !(STUSPS %in% c("HI", "AK", "PR"))]
# dt_outline_HI <- dt_states[STUSPS == "HI"]
# dt_outline_AK <- dt_states[STUSPS == "AK"]
# dt_outline_PR <- dt_states[STUSPS == "PR"]
# Zip3 data tables separated
dt_48 <- dt[ !(STUSPS %in% c("HI", "AK", "PR"))]
# dt_HI <- dt[STUSPS == "HI"]
# dt_AK <- dt[STUSPS == "AK"] 
# dt_PR <- dt[STUSPS == "PR"] 

# # dt_outline_48_mini <- dt_outline_48[ STUSPS %in% c("WA", "FL", "ME")]
# dt_48_mini <- dt_48[ STUSPS %in% c("WA", "FL", "ME")]



# 
us_outline <- st_read('data/gz_2010_us_outline_500k.shp')
us_outline <- data.table(us_outline)
# STATEFP is FIPS codes. PR is 72; AK is 02; HI is 15
# us_outline <- ms_simplify(us_outline, keep = 0.01)
# us_outline_HI <- us_outline[R_STATEFP == "15" | L_STATEFP == "15" ]
# us_outline_AK <- us_outline[R_STATEFP == "02" | L_STATEFP == "02" ]
# us_outline_PR <- us_outline[R_STATEFP == "72" | L_STATEFP == "72" ]
us_outline_48 <- us_outline[!(R_STATEFP %in% c("15", "02", "72")) & !(L_STATEFP %in% c("15", "02", "72")) ]


#######.
# Plot the data
######

#Plotting functions ##############
makeplots <- function(dt_zip, dt_outline, fill, max = 1){
  return(
    ggplot()+
      geom_sf(data = dt_outline, aes(geometry = geometry), fill = NA)+
      geom_sf(data = dt_zip, aes(geometry = geometry, fill = fill), color = NA) + 
      facet_grid(cols = vars(scenario), rows = vars(year)) +
      theme(panel.grid=element_blank(),
            axis.ticks=element_blank(),
            axis.text=element_blank(),
            axis.title = element_blank(),
            legend.position="bottom",
            panel.border=element_blank(),
            strip.background=element_blank(),
            strip.text=element_text(face="bold", hjust=0, size=10),
            panel.spacing = unit(0, "lines"))+
      scale_fill_gradientn(colours=c("#ffffff", brewer.pal(n=9, name="YlOrRd")),
                           na.value="#ffffff", name="% of iterations",
                           labels = percent, limits = c(0,max))
    
  )
  
}

annotation_custom2 <- 
  function (grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, data){
    layer(data = data, stat = StatIdentity, position = PositionIdentity, 
          geom = ggplot2:::GeomCustomAnn,
            inherit.aes = TRUE, params = list(grob = grob, 
                                              xmin = xmin, xmax = xmax, 
                                              ymin = ymin, ymax = ymax))
    }
  
gen_annotations <- function(row){
  return(annotation_custom2(
    grob = row$grobs[[1]],
    xmin = row$xmin[[1]],
    xmax = row$xmax[[1]],
    ymin = row$ymin[[1]],
    ymax = row$ymax[[1]],
    data = data.frame(year = row$year[[1]])))
}

##########

#
# WNVH
#


# gg_HI_WNVH <- makeplots(dt_HI, us_outline_HI, fill = dt_HI$WNVH) + coord_sf(xlim = c(-161, -154), ylim = c(19, 23))
# #ggsave("results/WNVH_HI_zip3.png", plot = gg_HI_WNVH, width = 6.5, height = 5)
# 
# gg_AK_WNVH <- makeplots(dt_AK, us_outline_AK, fill = dt_AK$WNVH) + coord_sf(xlim = c(-180, -129))
# #ggsave("Results/WNVH_AK_zip3.png", plot = gg_AK_WNVH, width = 6.5, height = 5)
# 
# gg_PR_WNVH <- makeplots(dt_PR, us_outline_PR, fill = dt_PR$WNVH)
# #ggsave("Results/WNVH_PR_zip3.png", plot = gg_PR_WNVH, width = 6.5, height = 5)

gg_48_WNVH <- makeplots(dt_48, us_outline_48, fill = dt_48$WNVH)
ggsave("WNVH_48_unif_compare.png", plot = gg_48_WNVH, width = 6.5, height = 6)